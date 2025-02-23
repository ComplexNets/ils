import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import streamlit as st
import plotly.graph_objects as go
from core.combine_fragments import combine_fragments, get_filtered_fragments
from core.heat_capacity import calculate_ionic_liquid_heat_capacity
from core.density import calculate_density, validate_density
from core.toxicity import calculate_ionic_liquid_toxicity
from core.solubility import calculate_solubility, validate_solubility
from core.pareto_optimizer import ParetoOptimizer
from frontend.property_input import PropertyRanges, PropertyCriteria
import pandas as pd
import numpy as np
import math
from multiprocessing import Pool, cpu_count
from functools import partial
from models.shortList_frag import fragments as short_fragments
from models.mediumList_frag import fragments as medium_fragments
from models.longList_frag import fragments as long_fragments
from utils.utils import generate_il_name, is_in_il_thermo
from utils.validation_rules import MolecularValidator
import random
import io

# Initialize session state
if 'property_ranges' not in st.session_state:
    st.session_state.property_ranges = PropertyRanges()
if 'optimizer' not in st.session_state:
    st.session_state.optimizer = ParetoOptimizer()
if 'fragment_list' not in st.session_state:
    st.session_state.fragment_list = 'short'
if 'validation_params' not in st.session_state:
    st.session_state.validation_params = {
        'max_total_groups': 6,
        'min_total_groups': 0,
        'max_groups_per_type': 6,
        'min_groups_per_type': 0,
        'max_alkyl_chain_length': 12,
        'max_alkyl_chains': 6,
        'max_groups_per_chain': 6,
        'min_groups_per_chain': 0,
        'max_group_occurrences': 6,
        'min_group_occurrences': 0,
        'max_total_group_occurrences': 6,
        'min_total_group_occurrences': 0
    }

def get_fragment_list_choice():
    """Get user's choice of fragment list"""
    st.sidebar.header("Fragment List Selection")
    
    # Radio button for fragment list selection
    list_choice = st.sidebar.radio(
        "Choose Fragment List:",
        ('Short List', 'Medium List', 'Long List'),
        index=0 if st.session_state.fragment_list == 'short' else 
              1 if st.session_state.fragment_list == 'medium' else 2,
        help="""
        Short List: Basic set of fragments
        Medium List: Intermediate set of fragments
        Long List: Extended set of fragments
        """
    )
    
    # Update session state
    if list_choice == 'Short List':
        st.session_state.fragment_list = 'short'
    elif list_choice == 'Medium List':
        st.session_state.fragment_list = 'medium'
    else:
        st.session_state.fragment_list = 'long'
    
    # Show fragment counts
    if st.session_state.fragment_list == 'short':
        fragments = short_fragments
    elif st.session_state.fragment_list == 'medium':
        fragments = medium_fragments
    else:
        fragments = long_fragments
        
    cation_count = len([f for f in fragments if f['fragment_type'] == 'cation'])
    anion_count = len([f for f in fragments if f['fragment_type'] == 'anion'])
    alkyl_count = len([f for f in fragments if f['fragment_type'] == 'alkyl_chain'])
    functional_group_count = len([f for f in fragments if f['fragment_type'] == 'functional_group'])
    
    st.sidebar.markdown(f"""
    **Fragment Counts:**
    - Cations: {cation_count}
    - Anions: {anion_count}
    - Alkyl Chains: {alkyl_count}
    - Functional Groups: {functional_group_count}
    - Total: {len(fragments)}
    """)
    
    return st.session_state.fragment_list

def get_user_ranges():
    """Get user-defined property ranges and update optimizer"""
    prop_ranges = st.session_state.property_ranges
    optimizer = st.session_state.optimizer
    
    # Add Validation Settings first
    st.sidebar.header("Structure Validation Settings")
    with st.sidebar.expander("Validation Rules", expanded=False):
        # Total groups constraints
        st.markdown("#### Total Groups Constraints")
        col1, col2 = st.columns(2)
        with col1:
            st.session_state.validation_params['min_total_groups'] = st.number_input(
                "Min Total Groups",
                min_value=0,
                max_value=10,
                value=st.session_state.validation_params['min_total_groups'],
                help="Minimum total number of functional groups allowed"
            )
        with col2:
            st.session_state.validation_params['max_total_groups'] = st.number_input(
                "Max Total Groups",
                min_value=1,
                max_value=10,
                value=st.session_state.validation_params['max_total_groups'],
                help="Maximum total number of functional groups allowed"
            )
            
        # Groups per type constraints
        st.markdown("#### Groups per Type Constraints")
        col1, col2 = st.columns(2)
        with col1:
            st.session_state.validation_params['min_groups_per_type'] = st.number_input(
                "Min Groups per Type",
                min_value=0,
                max_value=10,
                value=st.session_state.validation_params['min_groups_per_type'],
                help="Minimum number of groups of a specific type"
            )
        with col2:
            st.session_state.validation_params['max_groups_per_type'] = st.number_input(
                "Max Groups per Type",
                min_value=1,
                max_value=10,
                value=st.session_state.validation_params['max_groups_per_type'],
                help="Maximum number of groups of a specific type"
            )
            
        # Groups per chain constraints
        st.markdown("#### Groups per Chain Constraints")
        col1, col2 = st.columns(2)
        with col1:
            st.session_state.validation_params['min_groups_per_chain'] = st.number_input(
                "Min Groups per Chain",
                min_value=0,
                max_value=10,
                value=st.session_state.validation_params['min_groups_per_chain'],
                help="Minimum number of groups per chain"
            )
        with col2:
            st.session_state.validation_params['max_groups_per_chain'] = st.number_input(
                "Max Groups per Chain",
                min_value=1,
                max_value=10,
                value=st.session_state.validation_params['max_groups_per_chain'],
                help="Maximum number of groups per chain"
            )
            
        # Group occurrence constraints
        st.markdown("#### Group Occurrence Constraints")
        col1, col2 = st.columns(2)
        with col1:
            st.session_state.validation_params['min_group_occurrences'] = st.number_input(
                "Min Group Occurrences",
                min_value=0,
                max_value=10,
                value=st.session_state.validation_params['min_group_occurrences'],
                help="Minimum occurrences of any group"
            )
        with col2:
            st.session_state.validation_params['max_group_occurrences'] = st.number_input(
                "Max Group Occurrences",
                min_value=1,
                max_value=10,
                value=st.session_state.validation_params['max_group_occurrences'],
                help="Maximum occurrences of any group"
            )
            
        # Total group occurrence constraints
        st.markdown("#### Total Group Occurrence Constraints")
        col1, col2 = st.columns(2)
        with col1:
            st.session_state.validation_params['min_total_group_occurrences'] = st.number_input(
                "Min Total Occurrences",
                min_value=0,
                max_value=10,
                value=st.session_state.validation_params['min_total_group_occurrences'],
                help="Minimum total occurrences of all groups"
            )
        with col2:
            st.session_state.validation_params['max_total_group_occurrences'] = st.number_input(
                "Max Total Occurrences",
                min_value=1,
                max_value=10,
                value=st.session_state.validation_params['max_total_group_occurrences'],
                help="Maximum total occurrences of all groups"
            )
            
        # Alkyl chain constraints
        st.markdown("#### Alkyl Chain Constraints")
        col1, col2 = st.columns(2)
        with col1:
            st.session_state.validation_params['max_alkyl_chain_length'] = st.number_input(
                "Max Alkyl Chain Length",
                min_value=1,
                max_value=20,
                value=st.session_state.validation_params['max_alkyl_chain_length'],
                help="Maximum length of alkyl chains"
            )
        with col2:
            st.session_state.validation_params['max_alkyl_chains'] = st.number_input(
                "Max Alkyl Chains",
                min_value=1,
                max_value=10,
                value=st.session_state.validation_params['max_alkyl_chains'],
                help="Maximum number of alkyl chains"
            )
            
        # Chain length constraints
        col1, col2 = st.columns(2)
        with col1:
            min_length = st.number_input(
                "Min Chain Length",
                min_value=1,
                max_value=20,
                value=prop_ranges.validation.min_chain_length,
                help="Minimum number of carbon atoms in chain"
            )
            prop_ranges.validation.min_chain_length = min_length
        with col2:
            max_length = st.number_input(
                "Max Chain Length",
                min_value=1,
                max_value=20,
                value=prop_ranges.validation.max_chain_length,
                help="Maximum number of carbon atoms in chain"
            )
            prop_ranges.validation.max_chain_length = max_length
            
        # Group occurrence constraints
        st.markdown("#### Group Occurrence Constraints")
        col1, col2 = st.columns(2)
        with col1:
            upper_limit = st.number_input(
                "Upper Limit (t1)",
                min_value=0,
                max_value=10,
                value=prop_ranges.validation.group_occurrence_upper,
                help="Maximum number of occurrences of any group (Eq. 19)"
            )
            prop_ranges.validation.group_occurrence_upper = upper_limit
            
            lower_limit = st.number_input(
                "Lower Limit (t2)",
                min_value=0,
                max_value=10,
                value=prop_ranges.validation.group_occurrence_lower,
                help="Minimum number of occurrences of any group (Eq. 20)"
            )
            prop_ranges.validation.group_occurrence_lower = lower_limit
        
        with col2:
            exact_count = st.number_input(
                "Exact Count (t3)",
                min_value=-1,
                max_value=10,
                value=prop_ranges.validation.group_occurrence_exact,
                help="Exact number of occurrences required (-1 to disable) (Eq. 21)"
            )
            prop_ranges.validation.group_occurrence_exact = exact_count
    
    # Now add Property Ranges section
    st.sidebar.header("Property Ranges")
    
    # Density range (kg/m³)
    st.sidebar.subheader("Density (kg/m³)")
    density_col1, density_col2, density_col3 = st.sidebar.columns(3)
    with density_col1:
        density_min = st.number_input(
            "Minimum",
            value=float(prop_ranges.properties['density'].range[0]),
            step=50.0,
            format="%.1f",
            key="density_min"
        )
    with density_col2:
        density_max = st.number_input(
            "Maximum",
            value=float(prop_ranges.properties['density'].range[1]),
            step=50.0,
            format="%.1f",
            key="density_max"
        )
    with density_col3:
        density_importance = st.slider(
            "Importance",
            min_value=1,
            max_value=5,
            value=prop_ranges.properties['density'].importance,
            key="density_importance"
        )
    density_optimize_higher = st.sidebar.radio(
        "Density Optimization",
        ["Higher is better", "Lower is better"],
        index=0 if prop_ranges.properties['density'].optimize_higher else 1,
        horizontal=True,
        key="density_optimize"
    )

    # Heat capacity range (J/mol·K)
    st.sidebar.subheader("Heat Capacity (J/mol·K)")
    cp_col1, cp_col2, cp_col3 = st.sidebar.columns(3)
    with cp_col1:
        cp_min = st.number_input(
            "Minimum",
            value=float(prop_ranges.properties['heat_capacity'].range[0]),
            step=10.0,
            format="%.1f",
            key="cp_min"
        )
    with cp_col2:
        cp_max = st.number_input(
            "Maximum",
            value=float(prop_ranges.properties['heat_capacity'].range[1]),
            step=10.0,
            format="%.1f",
            key="cp_max"
        )
    with cp_col3:
        cp_importance = st.slider(
            "Importance",
            min_value=1,
            max_value=5,
            value=prop_ranges.properties['heat_capacity'].importance,
            key="cp_importance"
        )
    cp_optimize_higher = st.sidebar.radio(
        "Heat Capacity Optimization",
        ["Higher is better", "Lower is better"],
        index=0 if prop_ranges.properties['heat_capacity'].optimize_higher else 1,
        horizontal=True,
        key="cp_optimize"
    )

    # Toxicity range (IC50 in mM)
    st.sidebar.subheader("Toxicity (IC50 in mM)")
    toxicity_col1, toxicity_col2, toxicity_col3 = st.sidebar.columns(3)
    with toxicity_col1:
        toxicity_min = st.number_input(
            "Minimum IC50",
            value=float(prop_ranges.properties.get('toxicity', PropertyCriteria(range=(0.1, 100.0), importance=3, unit="mM")).range[0]),
            step=0.1,
            format="%.1f",
            key="toxicity_min"
        )
    with toxicity_col2:
        toxicity_max = st.number_input(
            "Maximum IC50",
            value=float(prop_ranges.properties.get('toxicity', PropertyCriteria(range=(0.1, 100.0), importance=3, unit="mM")).range[1]),
            step=0.1,
            format="%.1f",
            key="toxicity_max"
        )
    with toxicity_col3:
        toxicity_importance = st.slider(
            "Importance",
            min_value=1,
            max_value=5,
            value=prop_ranges.properties.get('toxicity', PropertyCriteria(range=(0.1, 100.0), importance=3, unit="mM")).importance,
            key="toxicity_importance"
        )
    toxicity_optimize_higher = st.sidebar.radio(
        "Toxicity Optimization",
        ["Higher is better (less toxic)", "Lower is better (more toxic)"],
        index=0 if prop_ranges.properties.get('toxicity').optimize_higher else 1,
        horizontal=True,
        key="toxicity_optimize"
    )

    # Solubility range (g/L)
    st.sidebar.subheader("Solubility (g/L)")
    solubility_col1, solubility_col2, solubility_col3 = st.sidebar.columns(3)
    with solubility_col1:
        solubility_min = st.number_input(
            "Minimum",
            value=float(prop_ranges.properties.get('solubility', PropertyCriteria(range=(0.1, 1000.0), importance=3, unit="g/L")).range[0]),
            step=0.1,
            format="%.1f",
            key="solubility_min"
        )
    with solubility_col2:
        solubility_max = st.number_input(
            "Maximum",
            value=float(prop_ranges.properties.get('solubility', PropertyCriteria(range=(0.1, 1000.0), importance=3, unit="g/L")).range[1]),
            step=1.0,
            format="%.1f",
            key="solubility_max"
        )
    with solubility_col3:
        solubility_importance = st.slider(
            "Importance",
            min_value=1,
            max_value=5,
            value=prop_ranges.properties.get('solubility', PropertyCriteria(range=(0.1, 1000.0), importance=3, unit="g/L")).importance,
            key="solubility_importance"
        )
    solubility_optimize_higher = st.sidebar.radio(
        "Solubility Optimization",
        ["Higher is better", "Lower is better"],
        index=0 if prop_ranges.properties.get('solubility', PropertyCriteria(range=(0.1, 1000.0), importance=3, unit="g/L")).optimize_higher else 1,
        horizontal=True,
        key="solubility_optimize"
    )

    # Update property ranges and optimizer constraints
    prop_ranges.update_property(
        'density',
        (density_min, density_max),
        weight=density_importance/5.0,
        optimize_higher=density_optimize_higher == "Higher is better"
    )
    optimizer.set_constraint(
        'density', 
        density_min, 
        density_max, 
        density_importance/5.0,
        optimize_higher=density_optimize_higher == "Higher is better"
    )
    
    prop_ranges.update_property(
        'heat_capacity',
        (cp_min, cp_max),
        weight=cp_importance/5.0,
        optimize_higher=cp_optimize_higher == "Higher is better"
    )
    optimizer.set_constraint(
        'heat_capacity', 
        cp_min, 
        cp_max, 
        cp_importance/5.0,
        optimize_higher=cp_optimize_higher == "Higher is better"
    )
    
    prop_ranges.update_property(
        'toxicity',
        (toxicity_min, toxicity_max),
        weight=toxicity_importance/5.0,
        optimize_higher=toxicity_optimize_higher == "Higher is better (less toxic)"
    )
    optimizer.set_constraint(
        'toxicity', 
        toxicity_min, 
        toxicity_max, 
        toxicity_importance/5.0,
        optimize_higher=toxicity_optimize_higher == "Higher is better (less toxic)"
    )
    
    prop_ranges.update_property(
        'solubility',
        (solubility_min, solubility_max),
        weight=solubility_importance/5.0,
        optimize_higher=solubility_optimize_higher == "Higher is better"
    )
    optimizer.set_constraint(
        'solubility', 
        solubility_min, 
        solubility_max, 
        solubility_importance/5.0,
        optimize_higher=solubility_optimize_higher == "Higher is better"
    )
    
    return prop_ranges

def calculate_properties():
    """Calculate properties for all valid combinations"""
    try:
        print("\n=== Starting Property Calculation ===")
        print(f"Using fragment list: {st.session_state.fragment_list}")
        
        # Clear previous results
        if 'combinations' in st.session_state:
            del st.session_state.combinations
        if 'pareto_front' in st.session_state:
            del st.session_state.pareto_front
        
        # Create progress containers
        progress_container = st.empty()
        with progress_container.container():
            # Progress indicators
            validation_status = st.empty()
            validation_progress = st.progress(0.0)
            
            # Get fragment list choice
            list_choice = st.session_state.fragment_list
            
            # Get validation parameters from session state
            validator_params = st.session_state.validation_params
            
            fragment_list = st.session_state.fragment_list
            print(f"\nStarting combine_fragments with list: {fragment_list}")
            
            valid_combinations = combine_fragments(
                status_text=validation_status,
                progress_bar=validation_progress,
                validator_params=validator_params,
                list_name=fragment_list
            )
            
            print(f"\nReturned from combine_fragments with {len(valid_combinations) if valid_combinations else 0} combinations")
            
            if not valid_combinations:
                print("No valid combinations found!")
                st.error("No valid ionic liquid combinations found. Try adjusting your criteria.")
                st.stop()
            
            total_feasible = len(valid_combinations)
            
            # Step 2: Calculate properties for valid combinations
            calculation_status = st.empty()
            calculation_progress = st.progress(0.0)
            
            combinations = []
            total = len(valid_combinations)
            
            for i, combo in enumerate(valid_combinations):
                progress = (i + 1) / total
                calculation_status.write(f"Calculating properties for combination {i+1}/{total}")
                calculation_progress.progress(progress)
                
                try:
                    # Calculate properties for this combination
                    heat_capacity = calculate_ionic_liquid_heat_capacity(combo)
                    if heat_capacity is None:
                        continue
                    
                    density = calculate_density(combo['cation'], combo['anion'], combo['alkyl_chain'])
                    if density is None:
                        continue
                    
                    toxicity_result = calculate_ionic_liquid_toxicity(combo)
                    if toxicity_result is None:
                        continue
                    
                    toxicity = toxicity_result.get('ic50_mm', 0.0)
                    
                    # Calculate solubility
                    print(f"\nCalculating solubility for {combo['name']}")
                    max_solubility = st.session_state.property_ranges.properties['solubility'].range[1]
                    solubility = calculate_solubility(combo['cation'], combo['anion'], combo['alkyl_chain'], 
                                                    max_solubility=max_solubility)
                    print(f"Solubility result: {solubility}")
                    if solubility is None:
                        print("Skipping combination due to None solubility")
                        continue
                    
                    # Check if properties are within user-defined ranges
                    solubility_range = st.session_state.property_ranges.properties['solubility'].range
                    print(f"Checking solubility {solubility} against range {solubility_range}")
                    
                    if not (st.session_state.property_ranges.properties['density'].range[0] <= density <= st.session_state.property_ranges.properties['density'].range[1]):
                        print("Skipping due to density range")
                        continue
                    
                    if not (st.session_state.property_ranges.properties['heat_capacity'].range[0] <= heat_capacity <= st.session_state.property_ranges.properties['heat_capacity'].range[1]):
                        print("Skipping due to heat capacity range")
                        continue
                    
                    if not (st.session_state.property_ranges.properties['toxicity'].range[0] <= toxicity <= st.session_state.property_ranges.properties['toxicity'].range[1]):
                        print("Skipping due to toxicity range")
                        continue
                        
                    if not (st.session_state.property_ranges.properties['solubility'].range[0] <= solubility <= st.session_state.property_ranges.properties['solubility'].range[1]):
                        print(f"Skipping due to solubility range: {solubility} not in {st.session_state.property_ranges.properties['solubility'].range}")
                        continue
                    
                    # Add to combinations list if all properties were calculated and within ranges
                    combinations.append({
                        'name': combo['name'],
                        'cation': combo['cation']['name'],
                        'anion': combo['anion']['name'],
                        'alkyl_chain': combo['alkyl_chain']['name'],
                        'molecular_weight': combo['molecular_weight'],
                        'heat_capacity': heat_capacity,
                        'density': density,
                        'toxicity': toxicity,
                        'solubility': solubility,
                        'in_ilthermo': combo.get('in_ilthermo', False)
                    })
                except Exception as e:
                    print(f"Error calculating properties for {combo['name']}: {str(e)}")
                    continue
        
        # Store the results in session state
        st.session_state.combinations = combinations
        st.session_state.total_feasible = total_feasible
        
        if not combinations:
            st.warning("No combinations found within the specified property ranges.")
            st.stop()
            
        # Step 3: Get Pareto front
        with st.spinner("Calculating Pareto front..."):
            pareto_front = st.session_state.optimizer.get_pareto_front(combinations)
            
            # Calculate Pareto scores for all combinations
            for combo in combinations:
                pareto_score = 0.0
                for prop_name, constraint in st.session_state.optimizer.properties.items():
                    val = combo.get(prop_name, 0)
                    norm_val = st.session_state.optimizer._normalize_property(val, constraint)
                    pareto_score += norm_val * constraint.weight
                combo['pareto_score'] = pareto_score / len(st.session_state.optimizer.properties)
            
            # Add scores to Pareto front solutions
            for solution in pareto_front:
                matching_combo = next(c for c in combinations if c['name'] == solution['name'])
                solution['pareto_score'] = matching_combo['pareto_score']
        
        # Set final status message
        calculation_status.write(f"✅ Found {len(combinations)} combinations within specified ranges")
        
        return combinations, pareto_front
        
    except Exception as e:
        st.error(f"Error calculating properties: {str(e)}")
        raise e

# Main UI layout
st.set_page_config(page_title="Ionic Liquid Optimizer", layout="wide")

# Title and description
st.title("Ionic Liquid Property Optimizer")
st.markdown("""
This tool helps design ionic liquids by combining molecular fragments and optimizing for desired properties.
Choose your fragment list and set property ranges in the sidebar to begin.
""")

# Get fragment list choice (before property ranges)
fragment_list = get_fragment_list_choice()

# Initialize property ranges in sidebar
get_user_ranges()

# Add calculate button to sidebar
if st.sidebar.button("Find Optimal Ionic Liquids", key="calculate_button"):
    combinations, pareto_front = calculate_properties()
    if not combinations:
        st.warning("No valid combinations found.")
        st.stop()
    
    # Create tabs for different visualizations and analysis
    tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs([
        "Parallel Coordinates", "Property Correlation", 
        "Property Distributions", "Lead Ionic Liquids", "Statistics", "Fragment Lists"
    ])

    with tab1:
        # Helper function for normalization
        def get_normalized_range(values):
            min_val = min(values) if values else 0
            max_val = max(values) if values else 1
            # If all values are the same, add a small range
            if min_val == max_val:
                min_val = min_val - 0.5
                max_val = max_val + 0.5
            return min_val, max_val

        # Parallel coordinates plot
        fig_parallel = go.Figure()
        
        # Get global min/max for each property across all combinations
        property_ranges = {}
        for prop_name in st.session_state.property_ranges.properties:
            values = [c[prop_name] for c in combinations]
            property_ranges[prop_name] = get_normalized_range(values)
        
        # Add non-Pareto solutions
        non_pareto = [c for c in combinations if c not in pareto_front]
        if non_pareto:
            dims = []
            for prop_name, prop in st.session_state.property_ranges.properties.items():
                values = [c[prop_name] for c in non_pareto]
                min_val, max_val = property_ranges[prop_name]
                
                if prop_name == 'toxicity':
                    # For toxicity (IC50), use log scale since values can span orders of magnitude
                    # Higher IC50 = lower toxicity = better
                    log_values = [math.log10(max(v, 0.1)) for v in values]  # Use 0.1 mM as minimum to avoid log(0)
                    log_min = math.log10(0.1)  # 0.1 mM minimum
                    log_max = math.log10(100)  # 100 mM maximum
                    dims.append(dict(range=[log_min, log_max],
                               label=f"{prop_name} (IC50, mM)",
                               values=log_values,
                               ticktext=[f"{10**x:.1f}" for x in range(int(log_min), int(log_max)+1)],
                               tickvals=list(range(int(log_min), int(log_max)+1))))
                else:
                    # Scale other properties normally
                    scaled_values = [(v - min_val) / (max_val - min_val) for v in values]
                    dims.append(dict(range=[0, 1],
                               label=f"{prop_name} ({prop.unit})",
                               values=scaled_values))
            
            fig_parallel.add_trace(go.Parcoords(
                line=dict(color='rgba(128,128,128,0.3)',
                         colorscale=[[0, 'rgba(128,128,128,0.3)'], 
                                   [1, 'rgba(128,128,128,0.3)']]),
                dimensions=dims
            ))
        
        # Add Pareto solutions
        if pareto_front:
            dims = []
            for prop_name, prop in st.session_state.property_ranges.properties.items():
                values = [c[prop_name] for c in pareto_front]
                min_val, max_val = property_ranges[prop_name]
                
                if prop_name == 'toxicity':
                    # For toxicity (IC50), use log scale normalization
                    log_values = [math.log10(max(v, 0.1)) for v in values]
                    log_min = math.log10(0.1)
                    log_max = math.log10(100)
                    dims.append(dict(range=[log_min, log_max],
                               label=f"{prop_name} (IC50, mM)",
                               values=log_values,
                               ticktext=[f"{10**x:.1f}" for x in range(int(log_min), int(log_max)+1)],
                               tickvals=list(range(int(log_min), int(log_max)+1))))
                else:
                    # Scale other properties normally
                    scaled_values = [(v - min_val) / (max_val - min_val) for v in values]
                    dims.append(dict(range=[0, 1],
                               label=f"{prop_name} ({prop.unit})",
                               values=scaled_values))
            
            fig_parallel.add_trace(go.Parcoords(
                line=dict(color='rgba(0,128,255,1)',
                         colorscale=[[0, 'rgba(0,128,255,1)'], 
                                   [1, 'rgba(0,128,255,1)']]),
                dimensions=dims,
                name='Pareto Solutions'
            ))
        
        fig_parallel.update_layout(
            title="Parallel Coordinates Plot of Properties",
            height=600,
            showlegend=True
        )
        st.plotly_chart(fig_parallel, use_container_width=True)

    with tab2:
        st.subheader("Property Correlation Analysis")
        # Create correlation matrix
        props = list(combinations[0].keys())
        props = [p for p in props if p in ['heat_capacity', 'density', 'toxicity', 'solubility']]
        
        corr_data = []
        for p1 in props:
            row = []
            for p2 in props:
                x = [c[p1] for c in combinations]
                y = [c[p2] for c in combinations]
                corr = np.corrcoef(x, y)[0, 1]
                row.append(corr)
            corr_data.append(row)
        
        fig_corr = go.Figure(data=go.Heatmap(
            z=corr_data,
            x=props,
            y=props,
            text=[[f"{val:.2f}" for val in row] for row in corr_data],
            texttemplate="%{text}",
            textfont={"size": 10},
            hoverongaps=False,
            colorscale="RdBu"
        ))
        
        fig_corr.update_layout(
            title="Property Correlation Matrix",
            height=400
        )
        
        # Create SPLOM
        fig_splom = go.Figure(data=go.Splom(
            dimensions=[dict(label=p, values=[c[p] for c in combinations]) for p in props],
            showupperhalf=False,
            diagonal_visible=False,
            text=[c.get('name', f'IL {i+1}') for i, c in enumerate(combinations)],
            hovertemplate='<b>%{text}</b><br>' +
                         '%{xaxis.title.text}: %{x:.2f}<br>' +
                         '%{yaxis.title.text}: %{y:.2f}<br>',
            marker=dict(
                color='white',
                size=6,
                line=dict(
                    color='darkgray',
                    width=1
                )
            )
        ))
        
        fig_splom.update_layout(
            title="Scatter Plot Matrix",
            height=600
        )
        
        # Display both plots
        st.plotly_chart(fig_corr, use_container_width=True)
        st.plotly_chart(fig_splom, use_container_width=True)

    with tab3:
        st.subheader("Property Distributions")
        # Create figure with secondary y-axis
        fig = go.Figure()
        
        # Heat Capacity Distribution
        fig.add_trace(
            go.Histogram(
                x=[sol['heat_capacity'] for sol in combinations],
                name="Heat Capacity",
                nbinsx=30,
                marker_color='blue',
                opacity=0.6
            )
        )
        
        # Density Distribution
        fig.add_trace(
            go.Histogram(
                x=[sol['density'] for sol in combinations],
                name="Density",
                nbinsx=30,
                marker_color='red',
                opacity=0.6
            )
        )
        
        # Toxicity Distribution
        fig.add_trace(
            go.Histogram(
                x=[sol['toxicity'] for sol in combinations],
                name="Toxicity",
                nbinsx=30,
                marker_color='green',
                opacity=0.6
            )
        )
        
        # Solubility Distribution
        fig.add_trace(
            go.Histogram(
                x=[sol['solubility'] for sol in combinations],
                name="Solubility",
                nbinsx=30,
                marker_color='yellow',
                opacity=0.6
            )
        )
        
        fig.update_layout(
            title="Property Distributions",
            xaxis_title="Value",
            yaxis_title="Count",
            barmode='overlay',
            showlegend=True,
            legend=dict(
                yanchor="top",
                y=0.99,
                xanchor="right",
                x=0.99
            )
        )
        st.plotly_chart(fig, use_container_width=True)

    with tab4:
        st.subheader("Lead Ionic Liquids")
        
        # Add ILThermo validation filter
        show_validated = st.checkbox("Show only ILThermo validated", key="show_validated")
        
        # Display the table of top solutions
        if filtered_solutions := [
            sol for sol in pareto_front 
            if (not show_validated or sol.get('in_ilthermo', False))
        ]:
            # Prepare data for display
            solution_data = []
            for sol in filtered_solutions[:10]:
                solution_data.append({
                    'Name': sol['name'],
                    'Heat Capacity (J/mol·K)': f"{sol['heat_capacity']:.1f}",
                    'Density (kg/m³)': f"{sol['density']:.1f}",
                    'Toxicity (IC50 in mM)': f"{sol['toxicity']:.1f}",
                    'Solubility (g/L)': f"{sol['solubility']:.1f}",
                    'Pareto Score': f"{sol.get('pareto_score', 0):.3f}",
                    'In ILThermo': 'Yes' if sol.get('in_ilthermo', False) else 'No'
                })
            
            # Display the table
            st.dataframe(
                pd.DataFrame(solution_data),
                use_container_width=True,
                hide_index=True
            )
        else:
            st.warning("No solutions match the current filters.")
            
        # Add export section
        st.divider()
        st.subheader("Export All Ionic Liquids")
        
        # Create the export data
        def create_excel_download():
            # Prepare data
            all_data = []
            for combo in combinations:
                all_data.append({
                    'Name': combo['name'],
                    'Heat Capacity (J/mol·K)': combo['heat_capacity'],
                    'Density (kg/m³)': combo['density'],
                    'Toxicity (IC50 in mM)': combo['toxicity'],
                    'Solubility (g/L)': combo['solubility'],
                    'Pareto Score': combo.get('pareto_score', 0),
                    'In ILThermo': 'Yes' if combo.get('in_ilthermo', False) else 'No'
                })
            
            # Create DataFrame
            df = pd.DataFrame(all_data)
            
            # Create Excel file in memory
            output = io.BytesIO()
            with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
                df.to_excel(writer, sheet_name='All Ionic Liquids', index=False)
                
                # Get the worksheet
                worksheet = writer.sheets['All Ionic Liquids']
                
                # Auto-adjust columns
                for idx, col in enumerate(df.columns):
                    series = df[col]
                    max_len = max(
                        series.astype(str).apply(len).max(),
                        len(str(col))
                    ) + 1
                    worksheet.set_column(idx, idx, max_len)
            
            return output.getvalue()
        
        # Create download button
        if len(combinations) > 0:
            try:
                excel_data = create_excel_download()
                st.download_button(
                    label=f"Download All Ionic Liquids ({len(combinations)} compounds)",
                    data=excel_data,
                    file_name="all_ionic_liquids.xlsx",
                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                )
                st.caption(f"File contains {len(combinations)} ionic liquids with their properties.")
            except Exception as e:
                st.error(f"Error creating Excel file: {str(e)}")
        else:
            st.warning("No ionic liquids available for export.")

    with tab5:
        st.subheader("Optimization Statistics")
        
        # Create columns for statistics
        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.metric(
                "Total Feasible ILs",
                f"{st.session_state.total_feasible}",
                help="Total number of chemically feasible ionic liquids"
            )
            st.metric(
                "Valid ILs (In Range)",
                f"{len(combinations)}",
                help="Number of ionic liquids within specified property ranges"
            )
            st.metric(
                "ILThermo Validated",
                f"{sum(1 for c in combinations if c['in_ilthermo'])}",
                help="Number of ionic liquids found in ILThermo database"
            )
        
        with col2:
            # Calculate property ranges in the results
            density_values = [s['density'] for s in combinations if 'density' in s]
            cp_values = [s['heat_capacity'] for s in combinations if 'heat_capacity' in s]
            toxicity_values = [s['toxicity'] for s in combinations if 'toxicity' in s]
            solubility_values = [s['solubility'] for s in combinations if 'solubility' in s]
            
            if density_values:
                st.metric("Density Range", 
                         f"{min(density_values):.1f} - {max(density_values):.1f} kg/m³")
            if cp_values:
                st.metric("Heat Capacity Range", 
                         f"{min(cp_values):.1f} - {max(cp_values):.1f} J/mol·K")
            if toxicity_values:
                st.metric("Toxicity Range", 
                         f"{min(toxicity_values):.1f} - {max(toxicity_values):.1f} mM")
            if solubility_values:
                st.metric("Solubility Range", 
                         f"{min(solubility_values):.1f} - {max(solubility_values):.1f} g/L")
                         
        with col3:
            # Show average scores
            if combinations:
                avg_density = sum(s.get('density', 0) for s in combinations) / len(combinations)
                avg_cp = sum(s.get('heat_capacity', 0) for s in combinations) / len(combinations)
                avg_toxicity = sum(s.get('toxicity', 0) for s in combinations) / len(combinations)
                avg_solubility = sum(s.get('solubility', 0) for s in combinations) / len(combinations)
                
                st.metric("Average Density", f"{avg_density:.1f} kg/m³")
                st.metric("Average Heat Capacity", f"{avg_cp:.1f} J/mol·K")
                st.metric("Average Toxicity", f"{avg_toxicity:.1f} mM")
                st.metric("Average Solubility", f"{avg_solubility:.1f} g/L")
        
        # Summary statistics in sidebar
        st.sidebar.subheader("Quick Summary")
        st.sidebar.write(f"Total combinations: {len(combinations)}")
        st.sidebar.write(f"Pareto-optimal solutions: {len(pareto_front)}")
        st.sidebar.write(f"ILThermo validated: {sum(1 for s in combinations if s.get('in_ilthermo', False))}")
        
        # Export results button
        if st.sidebar.button("Export Results", key="export_results"):
            csv_data = []
            for sol in pareto_front:
                csv_data.append({
                    'Name': sol['name'],
                    'Heat_Capacity': sol['heat_capacity'],
                    'Density': sol['density'],
                    'Toxicity': sol['toxicity'],
                    'Solubility': sol['solubility'],
                    'Pareto_Score': sol.get('pareto_score', 0),
                    'In_ILThermo': sol.get('in_ilthermo', False)
                })
            df = pd.DataFrame(csv_data)
            csv = df.to_csv(index=False)
            st.sidebar.download_button(
                "Download CSV",
                csv,
                "ionic_liquids.csv",
                "text/csv",
                key='download-csv'
            )

    with tab6:
        st.subheader("Fragment Lists")
        
        # Get the current fragment list
        current_list = st.session_state.get('fragment_list', [])
        
        # Create a DataFrame for better display
        fragments_data = []
        if st.session_state.fragment_list == 'short':
            fragments = short_fragments
        elif st.session_state.fragment_list == 'medium':
            fragments = medium_fragments
        else:
            fragments = long_fragments
        
        for frag in fragments:
            fragments_data.append({
                'Name': frag['name'],
                'SMILES': frag['smiles'],
                'Type': frag['fragment_type'].replace('_', ' ').title()
            })
        
        if fragments_data:
            # Convert to DataFrame
            import pandas as pd
            df = pd.DataFrame(fragments_data)
            
            # Group by fragment type
            for frag_type in sorted(df['Type'].unique()):
                st.write(f"### {frag_type}")
                type_df = df[df['Type'] == frag_type][['Name', 'SMILES']]
                st.dataframe(type_df, hide_index=True, use_container_width=True)
        else:
            st.warning("No fragments loaded. Please select a fragment list from the sidebar.")