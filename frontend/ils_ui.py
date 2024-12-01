import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import streamlit as st
import plotly.graph_objects as go
from core.combine_fragments import combine_fragments, get_filtered_fragments
from core.heat_capacity import calculate_ionic_liquid_heat_capacity
from core.density import calculate_density, validate_density
from core.toxicity import calculate_ionic_liquid_toxicity
from core.pareto_optimizer import ParetoOptimizer
from frontend.property_input import PropertyRanges, PropertyCriteria
import pandas as pd
import numpy as np
import math
from multiprocessing import Pool, cpu_count
from functools import partial
from models.shortList_frag import fragments
from utils.utils import generate_il_name, is_in_il_thermo
from utils.validation_rules import MolecularValidator
from utils.utils import generate_il_name

# Initialize session state
if 'property_ranges' not in st.session_state:
    st.session_state.property_ranges = PropertyRanges()
if 'optimizer' not in st.session_state:
    st.session_state.optimizer = ParetoOptimizer()

def get_user_ranges():
    """Get user-defined property ranges and update optimizer"""
    prop_ranges = st.session_state.property_ranges
    optimizer = st.session_state.optimizer
    
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
    
    return prop_ranges

def calculate_properties():
    """Calculate properties for all valid combinations"""
    try:
        # Get user-defined ranges
        prop_ranges = st.session_state.property_ranges
        optimizer = st.session_state.optimizer
        
        # Step 1: Get fragments and generate combinations
        validation_status = st.empty()
        validation_progress = st.progress(0.0)
        
        valid_combinations = combine_fragments(
            status_text=validation_status,
            progress_bar=validation_progress
        )
        
        if not valid_combinations:
            st.warning("No valid ionic liquid combinations found.")
            validation_status.empty()
            validation_progress.empty()
            return [], []
        
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
                
                # Check if properties are within user-defined ranges
                if not (prop_ranges.properties['density'].range[0] <= density <= prop_ranges.properties['density'].range[1]):
                    continue
                
                if not (prop_ranges.properties['heat_capacity'].range[0] <= heat_capacity <= prop_ranges.properties['heat_capacity'].range[1]):
                    continue
                
                if not (prop_ranges.properties['toxicity'].range[0] <= toxicity <= prop_ranges.properties['toxicity'].range[1]):
                    continue
                
                # Add to combinations list if all properties were calculated and within ranges
                combinations.append({
                    'name': combo['name'],
                    'cation': combo['cation']['name'],
                    'anion': combo['anion']['name'],
                    'alkyl_chain': combo['alkyl_chain']['name'],
                    'in_ilthermo': is_in_il_thermo(combo['name']),
                    'heat_capacity': float(heat_capacity),
                    'density': float(density),
                    'toxicity': float(toxicity)
                })
            except Exception as e:
                st.error(f"Error calculating properties for {combo['name']}: {str(e)}")
                continue
        
        if not combinations:
            st.warning("No combinations found within the specified property ranges.")
            validation_status.empty()
            validation_progress.empty()
            calculation_status.empty()
            calculation_progress.empty()
            return [], []
        
        calculation_status.write(f"✅ Found {len(combinations)} combinations within specified ranges")
        
        # Step 3: Get Pareto front
        with st.spinner("Calculating Pareto front..."):
            pareto_front = st.session_state.optimizer.get_pareto_front(combinations)
            
            # Calculate Pareto scores for all combinations
            for combo in combinations:
                pareto_score = 0.0
                for prop_name, constraint in optimizer.properties.items():
                    val = combo.get(prop_name, 0)
                    norm_val = optimizer._normalize_property(val, constraint)
                    pareto_score += norm_val * constraint.weight
                combo['pareto_score'] = pareto_score / len(optimizer.properties)
            
            # Add scores to Pareto front solutions
            for solution in pareto_front:
                matching_combo = next(c for c in combinations if c['name'] == solution['name'])
                solution['pareto_score'] = matching_combo['pareto_score']
        
        # Clean up progress indicators
        validation_status.empty()
        validation_progress.empty()
        calculation_status.empty()
        calculation_progress.empty()
        
        return combinations, pareto_front
        
    except Exception as e:
        st.error(f"Error calculating properties: {str(e)}")
        raise e

def plot_pareto_front(combinations, pareto_front, prop_ranges):
    """Create Pareto front visualizations"""
    st.subheader("Multi-Property Visualization")
    
    # Helper function for normalization
    def get_normalized_range(values):
        min_val = min(values) if values else 0
        max_val = max(values) if values else 1
        # If all values are the same, add a small range
        if min_val == max_val:
            min_val = min_val - 0.5
            max_val = max_val + 0.5
        return min_val, max_val
    
    # Create tabs for different visualization types
    viz_tab1, viz_tab2 = st.tabs(["Parallel Coordinates", "Radar Plot"])
    
    with viz_tab1:
        # Parallel coordinates plot
        fig_parallel = go.Figure()
        
        # Get global min/max for each property across all combinations
        property_ranges = {}
        for prop_name in prop_ranges.properties:
            values = [c[prop_name] for c in combinations]
            property_ranges[prop_name] = get_normalized_range(values)
        
        # Add non-Pareto solutions
        non_pareto = [c for c in combinations if c not in pareto_front]
        if non_pareto:
            dims = []
            for prop_name, prop in prop_ranges.properties.items():
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
            for prop_name, prop in prop_ranges.properties.items():
                values = [c[prop_name] for c in pareto_front]
                min_val, max_val = property_ranges[prop_name]
                
                if prop_name == 'toxicity':
                    # For toxicity (IC50), use log scale
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
                line=dict(color='rgba(255,0,0,1)',
                         colorscale=[[0, 'rgba(255,0,0,1)'], 
                                   [1, 'rgba(255,0,0,1)']]),
                dimensions=dims
            ))
        
        fig_parallel.update_layout(
            title="Parallel Coordinates Plot of Properties",
            height=600,
            showlegend=True
        )
        st.plotly_chart(fig_parallel, use_container_width=True)
    
    with viz_tab2:
        if not pareto_front:
            st.write("No Pareto solutions available for visualization")
            return
            
        num_solutions = min(5, len(pareto_front))
        selected_indices = st.multiselect(
            "Select solutions to compare (max 5)",
            range(len(pareto_front)),
            default=range(min(3, num_solutions)),
            format_func=lambda x: f"Solution {x+1}: {pareto_front[x].get('name', f'Combination {x+1}')}",
            key="pareto_solutions"
        )
        
        # Create figure container
        fig_radar = go.Figure()
        
        # Get labels for the radar plot (do this once, outside the loop)
        labels = []
        for prop_name, prop in prop_ranges.properties.items():
            if prop_name == 'toxicity':
                labels.append(f"{prop_name}\n(IC50, mM)")
            else:
                labels.append(f"{prop_name}\n({prop.unit})")
        
        if selected_indices:
            # Get global min/max for each property across all combinations
            property_ranges = {}
            for prop_name in prop_ranges.properties:
                values = [c[prop_name] for c in combinations]
                property_ranges[prop_name] = get_normalized_range(values)
            
            for idx in selected_indices:
                solution = pareto_front[idx]
                values = []
                
                for prop_name, prop in prop_ranges.properties.items():
                    val = solution[prop_name]
                    min_val, max_val = property_ranges[prop_name]
                    
                    if prop_name == 'toxicity':
                        # For toxicity (IC50), use log scale normalization
                        log_val = math.log10(max(val, 0.1))
                        log_min = math.log10(0.1)
                        log_max = math.log10(100)
                        norm_val = (log_val - log_min) / (log_max - log_min)
                    else:
                        if min_val == max_val:
                            norm_val = 1.0
                        else:
                            norm_val = (val - min_val) / (max_val - min_val)
                    
                    values.append(norm_val)
                
                fig_radar.add_trace(go.Scatterpolar(
                    r=values,
                    theta=labels,
                    name=solution.get('name', f'Solution {idx+1}'),
                    fill='toself'
                ))
        else:
            # Add an empty trace to keep the plot structure
            fig_radar.add_trace(go.Scatterpolar(
                r=[0] * len(labels),
                theta=labels,
                showlegend=False
            ))
        
        # Always update layout
        fig_radar.update_layout(
            polar=dict(
                radialaxis=dict(visible=True, range=[0, 1], 
                              ticktext=['0%', '25%', '50%', '75%', '100%'],
                              tickvals=[0, 0.25, 0.5, 0.75, 1])
            ),
            showlegend=True,
            title="Radar Plot of Selected Solutions"
        )
        st.plotly_chart(fig_radar, use_container_width=True)

def plot_property_correlation(combinations):
    """Create correlation matrix and SPLOM for properties"""
    if not combinations:
        return
        
    # Create correlation matrix
    props = list(combinations[0].keys())
    props = [p for p in props if p in ['heat_capacity', 'density', 'toxicity']]
    
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
        diagonal_visible=False
    ))
    
    fig_splom.update_layout(
        title="Scatter Plot Matrix",
        height=600
    )
    
    # Display both plots
    st.plotly_chart(fig_corr, use_container_width=True)
    st.plotly_chart(fig_splom, use_container_width=True)

def plot_property_distribution(combinations):
    """Create distribution plots for properties"""
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
    
    return fig

def display_results():
    """Display optimization results"""
    st.title("Ionic Liquid Property Optimization")
    
    # Calculate properties and get results
    combinations, pareto_front = calculate_properties()
    if not combinations:
        return
    
    # Create tabs for different visualizations
    tab1, tab2, tab3, tab4, tab5 = st.tabs([
        "Pareto Front", 
        "Property Correlation",
        "Property Distributions",
        "Top Solutions",
        "Statistics"
    ])
    
    with tab1:
        st.subheader("Property Distribution and Pareto Front")
        fig_pareto = plot_pareto_front(combinations, pareto_front, st.session_state.property_ranges)
    
    with tab2:
        st.subheader("Property Correlation Analysis")
        plot_property_correlation(combinations)
    
    with tab3:
        st.subheader("Property Distributions")
        fig_dist = plot_property_distribution(combinations)
        st.plotly_chart(fig_dist, use_container_width=True)
    
    with tab4:
        st.subheader("Top Ionic Liquid Combinations")
        
        # Add filters
        col1, col2 = st.columns(2)
        with col1:
            min_score = st.slider(
                "Minimum Pareto Score",
                min_value=0.0,
                max_value=1.0,
                value=0.0,
                step=0.1,
                key="min_pareto_score"
            )
        with col2:
            show_ilthermo = st.checkbox("Show only ILThermo validated", value=False, key="show_ilthermo")
        
        # Filter solutions
        filtered_solutions = [
            sol for sol in pareto_front 
            if sol.get('pareto_score', 0) >= min_score
            and (not show_ilthermo or sol.get('in_ilthermo', False))
        ]
        
        # Display solutions in a table
        if filtered_solutions:
            solution_data = []
            for sol in filtered_solutions[:10]:
                solution_data.append({
                    'Name': sol['name'],
                    'Heat Capacity (J/mol·K)': f"{sol['heat_capacity']:.1f}",
                    'Density (kg/m³)': f"{sol['density']:.1f}",
                    'Toxicity (IC50 in mM)': f"{sol['toxicity']:.1f}",
                    'Pareto Score': f"{sol.get('pareto_score', 0):.3f}",
                    'In ILThermo': '✓' if sol.get('in_ilthermo', False) else '✗'
                })
            
            st.dataframe(
                solution_data,
                use_container_width=True,
                hide_index=True
            )
        else:
            st.warning("No solutions match the current filters.")

    with tab5:
        st.subheader("Optimization Statistics")
        
        # Create three columns for statistics
        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.metric("Total Combinations", len(combinations))
            st.metric("Pareto-optimal Solutions", len(pareto_front))
            st.metric("ILThermo Validated", sum(1 for s in combinations if s.get('in_ilthermo', False)))
        
        with col2:
            # Calculate property ranges in the results
            density_values = [s['density'] for s in combinations if 'density' in s]
            cp_values = [s['heat_capacity'] for s in combinations if 'heat_capacity' in s]
            toxicity_values = [s['toxicity'] for s in combinations if 'toxicity' in s]
            
            if density_values:
                st.metric("Density Range", 
                         f"{min(density_values):.1f} - {max(density_values):.1f} kg/m³")
            if cp_values:
                st.metric("Heat Capacity Range", 
                         f"{min(cp_values):.1f} - {max(cp_values):.1f} J/mol·K")
            if toxicity_values:
                st.metric("Toxicity Range", 
                         f"{min(toxicity_values):.1f} - {max(toxicity_values):.1f} mM")
                         
        with col3:
            # Show average scores
            if combinations:
                avg_density = sum(s.get('density', 0) for s in combinations) / len(combinations)
                avg_cp = sum(s.get('heat_capacity', 0) for s in combinations) / len(combinations)
                avg_toxicity = sum(s.get('toxicity', 0) for s in combinations) / len(combinations)
                
                st.metric("Average Density", f"{avg_density:.1f} kg/m³")
                st.metric("Average Heat Capacity", f"{avg_cp:.1f} J/mol·K")
                st.metric("Average Toxicity", f"{avg_toxicity:.1f} mM")
        
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
        
# Main UI layout
st.set_page_config(page_title="Ionic Liquid Optimizer", layout="wide")

# Initialize property ranges in sidebar
get_user_ranges()

# Add calculate button to sidebar
if st.sidebar.button("Find Optimal Ionic Liquids", key="calculate_button"):
    display_results()