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
from core.hydrophobicity import calculate_ionic_liquid_hydrophobicity, validate_hydrophobicity
from core.viscosity import calculate_viscosity, validate_viscosity
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
    prop_ranges = PropertyRanges()
    prop_ranges.update_property('heat_capacity', (100.0, 2200.0), 0.5, True)
    prop_ranges.update_property('density', (100.0, 4000.0), 0.5, False)
    prop_ranges.update_property('toxicity', (0.01, 200.0), 0.5, True)
    prop_ranges.update_property('solubility', (0.001, 4000.0), 0.5, True)
    prop_ranges.update_property('hydrophobicity', (-20.0, 20.0), 0.5, True)
    prop_ranges.update_property('viscosity', (0.0005, 2000.0), 0.5, False)
    st.session_state.property_ranges = prop_ranges
else:
    prop_ranges = st.session_state.property_ranges

if 'optimizer' not in st.session_state:
    optimizer = ParetoOptimizer()
    st.session_state.optimizer = optimizer

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
            "Min",
            value=100.0,
            step=10.0,
            key="density_min"
        )
    with density_col2:
        density_max = st.number_input(
            "Max",
            value=4000.0,
            step=10.0,
            key="density_max"
        )
    with density_col3:
        density_importance = st.slider(
            "Importance",
            min_value=0,
            max_value=10,
            value=5,
            help="0 = ignore this property, 10 = highest importance",
            key="density_importance_slider"
        ) / 10.0
    
    density_optimize_higher = st.sidebar.radio(
        "Density Optimization",
        ["Higher is better", "Lower is better"],
        index=1,
        key="density_optimize"
    )
    
    # Heat Capacity
    st.sidebar.subheader("Heat Capacity (J/K·mol)") # Updated unit label
    cp_col1, cp_col2, cp_col3 = st.sidebar.columns(3)
    
    with cp_col1:
        cp_min = st.number_input(
            "Min",
            value=100.0,
            step=10.0,
            key="cp_min"
        )
    with cp_col2:
        cp_max = st.number_input(
            "Max",
            value=2200.0,
            step=10.0,
            key="cp_max"
        )
    with cp_col3:
        heat_capacity_importance = st.slider(
            "Importance",
            min_value=0,
            max_value=10,
            value=5,
            help="0 = ignore this property, 10 = highest importance",
            key="heat_capacity_importance_slider"
        ) / 10.0
        
    cp_optimize_higher = st.sidebar.radio(
        "Heat Capacity Optimization",
        ["Higher is better", "Lower is better"],
        index=0,
        key="cp_optimize"
    )

    # Toxicity range (IC50 in mM)
    st.sidebar.subheader("Toxicity (IC50 in mM)")
    toxicity_col1, toxicity_col2, toxicity_col3 = st.sidebar.columns(3)
    with toxicity_col1:
        toxicity_min = st.number_input(
            "Min IC50",
            value=0.01,
            step=0.1,
            format="%.1f",
            key="toxicity_min"
        )
    with toxicity_col2:
        toxicity_max = st.number_input(
            "Max IC50",
            value=200.0,
            step=0.1,
            format="%.1f",
            key="toxicity_max"
        )
    with toxicity_col3:
        toxicity_importance = st.slider(
            "Importance",
            min_value=0,
            max_value=10,
            value=5,
            help="0 = ignore this property, 10 = highest importance",
            key="toxicity_importance_slider"
        ) / 10.0
        
    toxicity_optimize_higher = st.sidebar.radio(
        "Toxicity Optimization",
        ["Higher is better (less toxic)", "Lower is better (more toxic)"],
        index=0,
        horizontal=True,
        key="toxicity_optimize"
    )

    # Solubility range (g/L)
    st.sidebar.subheader("Solubility (g/L)")
    solubility_col1, solubility_col2, solubility_col3 = st.sidebar.columns(3)
    with solubility_col1:
        solubility_min = st.number_input(
            "Min",
            value=0.001,
            step=0.1,
            format="%.1f",
            key="solubility_min"
        )
    with solubility_col2:
        solubility_max = st.number_input(
            "Max",
            value=4000.0,
            step=1.0,
            format="%.1f",
            key="solubility_max"
        )
    with solubility_col3:
        solubility_importance = st.slider(
            "Importance",
            min_value=0,
            max_value=10,
            value=5,
            help="0 = ignore this property, 10 = highest importance",
            key="solubility_importance_slider"
        ) / 10.0
        
    solubility_optimize_higher = st.sidebar.radio(
        "Solubility Optimization",
        ["Higher is better", "Lower is better"],
        index=0,
        horizontal=True,
        key="solubility_optimize"
    )

    # Hydrophobicity range (logP)
    st.sidebar.subheader("Hydrophobicity (logP)")
    hydrophobicity_col1, hydrophobicity_col2, hydrophobicity_col3 = st.sidebar.columns(3)
    with hydrophobicity_col1:
        hydrophobicity_min = st.number_input(
            "Min",
            value=-20.0,
            step=0.1,
            format="%.1f",
            key="hydrophobicity_min"
        )
    with hydrophobicity_col2:
        hydrophobicity_max = st.number_input(
            "Max",
            value=20.0,
            step=0.1,
            format="%.1f",
            key="hydrophobicity_max"
        )
    with hydrophobicity_col3:
        hydrophobicity_importance = st.slider(
            "Importance",
            min_value=0,
            max_value=10,
            value=5,
            help="0 = ignore this property, 10 = highest importance",
            key="hydrophobicity_importance_slider"
        ) / 10.0
        
    hydrophobicity_optimize_higher = st.sidebar.radio(
        "Hydrophobicity Optimization",
        ["Higher is better", "Lower is better"],
        index=0,
        horizontal=True,
        key="hydrophobicity_optimize"
    )

    # Viscosity range (Pa·s)
    st.sidebar.subheader("Viscosity (Pa·s)")
    viscosity_col1, viscosity_col2, viscosity_col3 = st.sidebar.columns(3)
    with viscosity_col1:
        viscosity_min = st.number_input(
            "Min",
            value=0.0005,
            step=0.0001,
            format="%.4f",
            key="viscosity_min"
        )
    with viscosity_col2:
        viscosity_max = st.number_input(
            "Max",
            value=2000.0,
            step=1.0,
            format="%.1f",
            key="viscosity_max"
        )
    with viscosity_col3:
        viscosity_importance = st.slider(
            "Importance",
            min_value=0,
            max_value=10,
            value=5,
            help="0 = ignore this property, 10 = highest importance",
            key="viscosity_importance_slider"
        ) / 10.0
        
    viscosity_optimize_higher = st.sidebar.radio(
        "Viscosity Optimization",
        ["Higher is better", "Lower is better"],
        index=1,
        horizontal=True,
        key="viscosity_optimize"
    )

    # Update property ranges and optimizer constraints
    prop_ranges.update_property(
        'heat_capacity',
        (cp_min, cp_max),
        weight=heat_capacity_importance,
        optimize_higher=cp_optimize_higher == "Higher is better"
    )
    optimizer.set_constraint(
        'heat_capacity', 
        cp_min, 
        cp_max, 
        heat_capacity_importance,
        optimize_higher=cp_optimize_higher == "Higher is better"
    )
    
    prop_ranges.update_property(
        'density',
        (density_min, density_max),
        weight=density_importance,
        optimize_higher=density_optimize_higher == "Higher is better"
    )
    optimizer.set_constraint(
        'density', 
        density_min, 
        density_max, 
        density_importance,
        optimize_higher=density_optimize_higher == "Higher is better"
    )
    
    prop_ranges.update_property(
        'toxicity',
        (toxicity_min, toxicity_max),
        weight=toxicity_importance,
        optimize_higher=toxicity_optimize_higher == "Higher is better"
    )
    optimizer.set_constraint(
        'toxicity', 
        toxicity_min, 
        toxicity_max, 
        toxicity_importance,
        optimize_higher=toxicity_optimize_higher == "Higher is better"
    )
    
    prop_ranges.update_property(
        'solubility',
        (solubility_min, solubility_max),
        weight=solubility_importance,
        optimize_higher=solubility_optimize_higher == "Higher is better"
    )
    optimizer.set_constraint(
        'solubility', 
        solubility_min, 
        solubility_max, 
        solubility_importance,
        optimize_higher=solubility_optimize_higher == "Higher is better"
    )
    
    prop_ranges.update_property(
        'hydrophobicity',
        (hydrophobicity_min, hydrophobicity_max),
        weight=hydrophobicity_importance,
        optimize_higher=hydrophobicity_optimize_higher == "Higher is better"
    )
    optimizer.set_constraint(
        'hydrophobicity', 
        hydrophobicity_min, 
        hydrophobicity_max, 
        hydrophobicity_importance,
        optimize_higher=hydrophobicity_optimize_higher == "Higher is better"
    )
    
    prop_ranges.update_property(
        'viscosity',
        (viscosity_min, viscosity_max),
        weight=viscosity_importance,
        optimize_higher=viscosity_optimize_higher == "Higher is better"
    )
    optimizer.set_constraint(
        'viscosity', 
        viscosity_min, 
        viscosity_max, 
        viscosity_importance,
        optimize_higher=viscosity_optimize_higher == "Higher is better"
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
                    heat_capacity_results = calculate_ionic_liquid_heat_capacity(combo)
                    if heat_capacity_results is None:
                        print(f"Skipping {combo.get('name', 'Unknown')} due to None heat capacity result.")
                        continue
                    heat_capacity_molar = heat_capacity_results.get('molar')
                    # heat_capacity_gravimetric = heat_capacity_results.get('gravimetric') # Available if needed later

                    if heat_capacity_molar is None:
                         print(f"Skipping {combo.get('name', 'Unknown')} due to missing molar heat capacity.")
                         continue

                    density = calculate_density(
                        cation=combo['cation'], 
                        anion=combo['anion'], 
                        alkyl_chain=combo['alkyl_chain'],
                        functional_group=combo.get('functional_group')
                    )
                    if density is None:
                        continue
                    
                    toxicity_result = calculate_ionic_liquid_toxicity(combo)
                    if toxicity_result is None:
                        continue
                    
                    toxicity = toxicity_result.get('ic50_mm', 0.0)
                    
                    # Calculate solubility
                    print(f"\nCalculating solubility for {combo['name']}")
                    max_solubility = st.session_state.property_ranges.properties['solubility'].range[1]
                    solubility = calculate_solubility(
                        cation=combo['cation'], 
                        anion=combo['anion'], 
                        alkyl_chain=combo['alkyl_chain'],
                        functional_group=combo.get('functional_group'),
                        max_solubility=max_solubility
                    )
                    print(f"Solubility result: {solubility}")
                    if solubility is None:
                        print("Skipping combination due to None solubility")
                        continue
                    
                    hydrophobicity = calculate_ionic_liquid_hydrophobicity(combo)
                    if hydrophobicity is None:
                        continue
                    
                    # Extract individual components for viscosity calculation
                    viscosity = calculate_viscosity(
                        cation=combo['cation'],
                        anion=combo['anion'],
                        alkyl_chain=combo['alkyl_chain'],
                        functional_group=combo.get('functional_group')
                    )
                    if viscosity is None:
                        continue

                    # Check if properties are within user-defined ranges
                    solubility_range = st.session_state.property_ranges.properties['solubility'].range
                    print(f"Checking solubility {solubility} against range {solubility_range}")

                    if not (st.session_state.property_ranges.properties['density'].range[0] <= density <= st.session_state.property_ranges.properties['density'].range[1]):
                        print(f"Skipping {combo.get('name', 'Unknown')} due to density range ({density:.2f})")
                        continue

                    if not (st.session_state.property_ranges.properties['heat_capacity'].range[0] <= heat_capacity_molar <= st.session_state.property_ranges.properties['heat_capacity'].range[1]):
                        print(f"Skipping {combo.get('name', 'Unknown')} due to heat capacity range ({heat_capacity_molar:.2f})")
                        continue

                    if not (st.session_state.property_ranges.properties['toxicity'].range[0] <= toxicity <= st.session_state.property_ranges.properties['toxicity'].range[1]):
                        print("Skipping due to toxicity range")
                        continue
                        
                    if not (st.session_state.property_ranges.properties['solubility'].range[0] <= solubility <= st.session_state.property_ranges.properties['solubility'].range[1]):
                        print(f"Skipping due to solubility range: {solubility} not in {st.session_state.property_ranges.properties['solubility'].range}")
                        continue
                    
                    if not (st.session_state.property_ranges.properties['hydrophobicity'].range[0] <= hydrophobicity <= st.session_state.property_ranges.properties['hydrophobicity'].range[1]):
                        print(f"Skipping due to hydrophobicity range: {hydrophobicity} not in {st.session_state.property_ranges.properties['hydrophobicity'].range}")
                        continue
                    
                    if not (st.session_state.property_ranges.properties['viscosity'].range[0] <= viscosity <= st.session_state.property_ranges.properties['viscosity'].range[1]):
                        print(f"Skipping due to viscosity range: {viscosity} not in {st.session_state.property_ranges.properties['viscosity'].range}")
                        continue
                    
                    # Add to combinations list if all properties were calculated and within ranges
                    combo_dict = {
                        'name': combo['name'],
                        'cation': combo['cation'],
                        'anion': combo['anion'],
                        'alkyl_chain': combo['alkyl_chain'],
                        'heat_capacity': heat_capacity_molar, # Use molar value
                        # 'heat_capacity_gravimetric': heat_capacity_gravimetric, # Store if needed
                        'density': density,
                        'toxicity': toxicity,
                        'solubility': solubility,
                        'hydrophobicity': hydrophobicity,
                        'viscosity': viscosity,
                        'in_ilthermo': combo.get('in_ilthermo', False),  # Copy the in_ilthermo property
                        # --- Explicitly copy canonical SMILES keys if they exist ---
                        'cation_canonical_smiles': combo.get('cation_canonical_smiles'),
                        'anion_canonical_smiles': combo.get('anion_canonical_smiles'),
                        'alkyl_chain_canonical_smiles': combo.get('alkyl_chain_canonical_smiles')
                    }
                    if 'functional_group_canonical_smiles' in combo:
                         combo_dict['functional_group_canonical_smiles'] = combo.get('functional_group_canonical_smiles')
                    # --- End copying SMILES ---

                    # Add functional group if present
                    if 'functional_group' in combo and combo['functional_group']:
                        combo_dict['functional_group'] = combo['functional_group']
                        
                    combinations.append(combo_dict)
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
            pareto_front = st.session_state.optimizer.find_pareto_front(combinations)
            
            # Calculate Pareto scores for all combinations
            for combo in combinations:
                total_score = 0
                total_weight = 0
                
                # Calculate weighted score for each property
                for prop_name, constraint in st.session_state.optimizer.constraints.items():
                    if prop_name not in combo:
                        continue
                        
                    value = combo[prop_name]
                    weight = constraint['weight']
                    
                    # Skip if weight is 0 (property is ignored)
                    if weight == 0:
                        continue
                    
                    # Normalize the value
                    norm_value = st.session_state.optimizer.normalize_value(
                        value, 
                        constraint['min'], 
                        constraint['max']
                    )
                    
                    # Invert if we want to minimize
                    if not constraint['optimize_higher']:
                        norm_value = 1.0 - norm_value
                    
                    # Add to weighted sum
                    total_score += norm_value * weight
                    total_weight += weight
                
                # Calculate final score
                combo['pareto_score'] = total_score / total_weight if total_weight > 0 else 0.0
        
            # Add scores to Pareto front solutions
            for solution in pareto_front:
                matching_combo = next(c for c in combinations if c['name'] == solution['name'])
                solution['pareto_score'] = matching_combo['pareto_score']
            
            # Store results in session state
            st.session_state.pareto_front = pareto_front
            
            # Print summary
            print(f"\nOptimization complete!")
            print(f"Total combinations: {len(combinations)}")
            print(f"Pareto front size: {len(pareto_front)}")
            for prop_name, constraint in st.session_state.optimizer.constraints.items():
                print(f"\n{prop_name}:")
                print(f"  Range: {constraint['min']} to {constraint['max']}")
                print(f"  Weight: {constraint['weight']}")
                print(f"  Optimize higher: {constraint['optimize_higher']}")
            
            return combinations, pareto_front
        
    except Exception as e:
        st.error(f"Error calculating properties: {str(e)}")
        raise e

# Main UI layout
st.set_page_config(page_title="Ionic Liquid Optimizer", layout="wide")

# Set sidebar width
st.markdown(
    """
    <style>
    [data-testid="stSidebar"][aria-expanded="true"]{
        min-width: 450px;
        max-width: 450px;
    }
    </style>
    """,
    unsafe_allow_html=True,
)

# Title and description
st.title("IL-SCOPE")
st.markdown("""
The IL Screening, Computational Optimization, and Property Estimation tool helps design ionic liquids by combining molecular fragments and optimizing for desired properties. Choose your fragment list and set property ranges in the sidebar to begin.
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
            if not values:
                return 0, 1
            min_val = min(values)
            max_val = max(values)
            # If all values are the same or very close
            if abs(max_val - min_val) < 1e-10:
                # Use property range as fallback
                prop_name = next((name for name, prop in st.session_state.property_ranges.properties.items() 
                                if any(abs(v - values[0]) < 1e-10 for v in [prop.range[0], prop.range[1]]))
                              , None)
                if prop_name:
                    min_val, max_val = st.session_state.property_ranges.properties[prop_name].range
                else:
                    # If we can't find the property, add a small range around the value
                    min_val = values[0] * 0.9 if values[0] != 0 else -0.1
                    max_val = values[0] * 1.1 if values[0] != 0 else 0.1
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
                elif prop_name == 'hydrophobicity':
                    # For hydrophobicity (logP), use log scale since values can span orders of magnitude
                    log_values = [math.log10(max(v, 0.1)) for v in values]  # Use 0.1 as minimum to avoid log(0)
                    log_min = math.log10(0.1)  # 0.1 minimum
                    log_max = math.log10(10)  # 10 maximum
                    dims.append(dict(range=[log_min, log_max],
                               label=f"{prop_name} (logP)",
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
                elif prop_name == 'hydrophobicity':
                    # For hydrophobicity (logP), use log scale normalization
                    log_values = [math.log10(max(v, 0.1)) for v in values]
                    log_min = math.log10(0.1)
                    log_max = math.log10(10)
                    dims.append(dict(range=[log_min, log_max],
                               label=f"{prop_name} (logP)",
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
        props = [p for p in props if p in ['heat_capacity', 'density', 'toxicity', 'solubility', 'hydrophobicity', 'viscosity']]
        
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
        
        # Hydrophobicity Distribution
        fig.add_trace(
            go.Histogram(
                x=[sol['hydrophobicity'] for sol in combinations],
                name="Hydrophobicity",
                nbinsx=30,
                marker_color='orange',
                opacity=0.6
            )
        )
        
        # Viscosity Distribution
        fig.add_trace(
            go.Histogram(
                x=[sol['viscosity'] for sol in combinations],
                name="Viscosity",
                nbinsx=30,
                marker_color='purple',
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
        
        # Get all solutions that match property ranges
        valid_solutions = []
        for sol in combinations:  
            # Check if solution is within property ranges
            is_valid = True
            for prop_name, prop in st.session_state.property_ranges.properties.items():
                if prop_name not in sol:
                    continue
                value = sol[prop_name]
                min_val, max_val = prop.range
                
                # Use small tolerance for zero values
                if min_val == 0 and max_val == 0:
                    if abs(value) > 1e-10:
                        is_valid = False
                        break
                else:
                    # Add small tolerance to bounds
                    if not (min_val - 1e-10 <= value <= max_val + 1e-10):
                        is_valid = False
                        break
            
            if is_valid and (not show_validated or sol.get('in_ilthermo', False)):
                valid_solutions.append(sol)
        
        # Sort solutions by Pareto score
        valid_solutions.sort(key=lambda x: x.get('pareto_score', 0), reverse=True)
        
        # Display the table of solutions
        if valid_solutions:
            st.write(f"Found {len(valid_solutions)} ionic liquids matching the criteria")
            
            # Prepare data for display
            solution_data = []
            for sol in valid_solutions[:10]:  
                solution_data.append({
                    'Name': sol['name'],
                    'Heat Capacity (J/mol·K)': f"{sol['heat_capacity']:.1f}",
                    'Density (kg/m³)': f"{sol['density']:.1f}",
                    'Toxicity (IC50 in mM)': f"{sol['toxicity']:.1f}",
                    'Solubility (g/L)': f"{sol['solubility']:.1f}",
                    'Hydrophobicity (logP)': f"{sol['hydrophobicity']:.1f}",
                    'Viscosity (Pa·s)': f"{sol['viscosity']:.4f}",
                    'Pareto Score': f"{sol.get('pareto_score', 0):.3f}",
                    'In ILThermo': 'Yes' if sol.get('in_ilthermo', False) else 'No'
                })
            
            # Display the table
            st.dataframe(
                pd.DataFrame(solution_data),
                use_container_width=True,
                hide_index=True
            )

            # --- Add table for Component Canonical SMILES ---
            st.subheader("Component Canonical SMILES for Top 10")
            smiles_data = []
            for sol in valid_solutions[:10]:
                 smiles_entry = {'Name': sol['name']}
                 # Add component smiles if they exist in the sol dictionary
                 for comp_key in ['cation_canonical_smiles', 'anion_canonical_smiles', 'alkyl_chain_canonical_smiles', 'functional_group_canonical_smiles']:
                      if comp_key in sol:
                           # Create a more readable column name
                           col_name = comp_key.replace('_canonical_smiles', '').replace('_', ' ').title() + ' SMILES'
                           smiles_entry[col_name] = sol[comp_key]
                 smiles_data.append(smiles_entry)

            if smiles_data:
                 st.dataframe(
                      pd.DataFrame(smiles_data),
                      use_container_width=True,
                      hide_index=True
                 )
            # --- End SMILES table ---

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
                    'Cation': combo['cation']['name'] if 'cation' in combo else None,
                    'Anion': combo['anion']['name'] if 'anion' in combo else None,
                    'Alkyl Chain': combo['alkyl_chain']['name'] if 'alkyl_chain' in combo else None,
                    'Functional Group': combo['functional_group']['name'] if 'functional_group' in combo and combo['functional_group'] else None,
                    'Heat Capacity (J/mol·K)': combo['heat_capacity'],
                    'Density (kg/m³)': combo['density'],
                    'Toxicity (IC50 in mM)': combo['toxicity'],
                    'Solubility (g/L)': combo['solubility'],
                    'Hydrophobicity (logP)': combo['hydrophobicity'],
                    'Viscosity (Pa·s)': combo['viscosity'],
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
                "Total Combinations Generated",
                f"{len(combinations)}",
                help="Total number of valid ionic liquid combinations found"
            )
            
            # Calculate percentage of combinations in ILThermo database
            total_in_ilthermo = sum(1 for c in combinations if c.get('in_ilthermo', False))
            percent_in_ilthermo = (total_in_ilthermo / len(combinations) * 100) if combinations else 0
            
            st.metric(
                "Combinations in ILThermo",
                f"{total_in_ilthermo}",
                f"{percent_in_ilthermo:.1f}% of total",
                help="Number of combinations that exist in the ILThermo database"
            )
        
        with col2:
            # Calculate property ranges in the results
            density_values = [s['density'] for s in combinations if 'density' in s]
            cp_values = [s['heat_capacity'] for s in combinations if 'heat_capacity' in s]
            toxicity_values = [s['toxicity'] for s in combinations if 'toxicity' in s]
            solubility_values = [s['solubility'] for s in combinations if 'solubility' in s]
            hydrophobicity_values = [s['hydrophobicity'] for s in combinations if 'hydrophobicity' in s]
            viscosity_values = [s['viscosity'] for s in combinations if 'viscosity' in s]
            
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
            if hydrophobicity_values:
                st.metric("Hydrophobicity Range", 
                         f"{min(hydrophobicity_values):.1f} - {max(hydrophobicity_values):.1f} logP")
            if viscosity_values:
                st.metric("Viscosity Range", 
                         f"{min(viscosity_values):.4f} - {max(viscosity_values):.4f} Pa·s")
                         
        with col3:
            # Show average scores
            if combinations:
                avg_density = sum(s.get('density', 0) for s in combinations) / len(combinations)
                avg_cp = sum(s.get('heat_capacity', 0) for s in combinations) / len(combinations)
                avg_toxicity = sum(s.get('toxicity', 0) for s in combinations) / len(combinations)
                avg_solubility = sum(s.get('solubility', 0) for s in combinations) / len(combinations)
                avg_hydrophobicity = sum(s.get('hydrophobicity', 0) for s in combinations) / len(combinations)
                avg_viscosity = sum(s.get('viscosity', 0) for s in combinations) / len(combinations)
                
                st.metric("Average Density", f"{avg_density:.1f} kg/m³")
                st.metric("Average Heat Capacity", f"{avg_cp:.1f} J/mol·K")
                st.metric("Average Toxicity", f"{avg_toxicity:.1f} mM")
                st.metric("Average Solubility", f"{avg_solubility:.1f} g/L")
                st.metric("Average Hydrophobicity", f"{avg_hydrophobicity:.1f} logP")
                st.metric("Average Viscosity", f"{avg_viscosity:.4f} Pa·s")
        
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
                    'Hydrophobicity': sol['hydrophobicity'],
                    'Viscosity': sol['viscosity'],
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
