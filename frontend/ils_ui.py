import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import streamlit as st
import plotly.graph_objects as go
from core.combine_fragments import combine_fragments
from core.heat_capacity import calculate_ionic_liquid_heat_capacity
from core.density import calculate_ionic_liquid_density
from core.toxicity import calculate_ionic_liquid_toxicity
from core.pareto_optimizer import ParetoOptimizer
from frontend.property_input import PropertyRanges, PropertyCriteria
import pandas as pd
import numpy as np
import math

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
    density_min = st.sidebar.number_input(
        "Minimum Density",
        value=prop_ranges.properties['density'].range[0],
        step=100
    )
    density_max = st.sidebar.number_input(
        "Maximum Density",
        value=prop_ranges.properties['density'].range[1],
        step=100
    )
    density_importance = st.sidebar.slider(
        "Density Importance",
        min_value=1,
        max_value=5,
        value=prop_ranges.properties['density'].importance
    )
    density_strict = st.sidebar.checkbox("Strict Density Constraint", value=False)
    
    # Heat capacity range (J/mol·K)
    st.sidebar.subheader("Heat Capacity (J/mol·K)")
    cp_min = st.sidebar.number_input(
        "Minimum Heat Capacity",
        value=prop_ranges.properties['heat_capacity'].range[0],
        step=50
    )
    cp_max = st.sidebar.number_input(
        "Maximum Heat Capacity",
        value=prop_ranges.properties['heat_capacity'].range[1],
        step=50
    )
    cp_importance = st.sidebar.slider(
        "Heat Capacity Importance",
        min_value=1,
        max_value=5,
        value=prop_ranges.properties['heat_capacity'].importance
    )
    cp_strict = st.sidebar.checkbox("Strict Heat Capacity Constraint", value=False)
    
    # Toxicity range (IC50 in mM)
    st.sidebar.subheader("Toxicity (IC50 in mM)")
    toxicity_min = st.sidebar.number_input(
        "Minimum IC50",
        value=prop_ranges.properties.get('toxicity', PropertyCriteria(range=(0.1, 100.0), importance=3, unit="mM")).range[0],
        step=0.1,
        format="%.1f"
    )
    toxicity_max = st.sidebar.number_input(
        "Maximum IC50",
        value=prop_ranges.properties.get('toxicity', PropertyCriteria(range=(0.1, 100.0), importance=3, unit="mM")).range[1],
        step=0.1,
        format="%.1f"
    )
    toxicity_importance = st.sidebar.slider(
        "Toxicity Importance",
        min_value=1,
        max_value=5,
        value=prop_ranges.properties.get('toxicity', PropertyCriteria(range=(0.1, 100.0), importance=3, unit="mM")).importance
    )
    toxicity_strict = st.sidebar.checkbox("Strict Toxicity Constraint", value=False)
    
    # Update property ranges
    prop_ranges.update_property(
        'density',
        (density_min, density_max),
        weight=density_importance/5.0,
        is_strict=density_strict
    )
    
    prop_ranges.update_property(
        'heat_capacity',
        (cp_min, cp_max),
        weight=cp_importance/5.0,
        is_strict=cp_strict
    )
    
    prop_ranges.update_property(
        'toxicity',
        (toxicity_min, toxicity_max),
        weight=toxicity_importance/5.0,
        is_strict=toxicity_strict
    )
    
    # Update optimizer constraints
    optimizer.set_constraint('density', density_min, density_max, 
                           density_importance/5.0, density_strict)
    optimizer.set_constraint('heat_capacity', cp_min, cp_max,
                           cp_importance/5.0, cp_strict)
    optimizer.set_constraint('toxicity', toxicity_min, toxicity_max,
                           toxicity_importance/5.0, toxicity_strict, inverse=True)
    
    return prop_ranges

def calculate_properties():
    """Calculate properties for all valid combinations"""
    try:
        combinations = combine_fragments()
        if not combinations:
            st.error("No valid combinations found!")
            return
            
        # Calculate properties for each combination
        for combo in combinations:
            combo['heat_capacity'] = calculate_ionic_liquid_heat_capacity(combo)
            combo['density'] = calculate_ionic_liquid_density(combo)
            toxicity_result = calculate_ionic_liquid_toxicity(combo)
            if toxicity_result:
                combo['toxicity'] = toxicity_result.get('ic50_mm', 0.0)
            else:
                combo['toxicity'] = 0.0
                
        # Get Pareto front
        optimizer = st.session_state.optimizer
        
        # Set up property constraints
        prop_ranges = st.session_state.property_ranges
        
        # Update optimizer with current property ranges and weights
        for prop_name, criteria in prop_ranges.properties.items():
            is_inverse = prop_name == 'toxicity'  # Higher IC50 is better for toxicity
            optimizer.set_constraint(
                prop_name,
                criteria.range[0],
                criteria.range[1],
                criteria.importance / 5.0,  # Convert 1-5 scale to 0-1
                is_strict=False,
                inverse=is_inverse
            )
        
        pareto_front = optimizer.get_pareto_front(combinations)
        
        # Calculate Pareto scores for all solutions
        ranked_solutions = optimizer.rank_solutions(combinations)
        for solution in combinations:
            matching_ranked = next((s for s in ranked_solutions if s['name'] == solution['name']), None)
            if matching_ranked:
                solution['pareto_score'] = matching_ranked.get('pareto_score', 0.0)
        
        return combinations, pareto_front
        
    except Exception as e:
        st.error(f"Error calculating properties: {str(e)}")
        return None, None

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
            format_func=lambda x: f"Solution {x+1}: {pareto_front[x].get('name', f'Combination {x+1}')}"
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
    st.title("Ionic Liquid Property Optimization")
    
    with st.spinner('Optimizing ionic liquids...'):
        try:
            # Calculate properties
            combinations, pareto_front = calculate_properties()
            if not combinations:
                st.warning("No valid ionic liquid combinations found.")
                return
            
            # Create tabs for different visualizations
            tab1, tab2, tab3, tab4 = st.tabs([
                "Pareto Front", 
                "Property Correlation",
                "Property Distributions",
                "Top Solutions"
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
                        step=0.1
                    )
                with col2:
                    show_ilthermo = st.checkbox("Show only ILThermo validated", value=False)
                
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
            
            # Summary statistics in sidebar
            st.sidebar.subheader("Optimization Summary")
            st.sidebar.write(f"Total combinations: {len(combinations)}")
            st.sidebar.write(f"Pareto-optimal solutions: {len(pareto_front)}")
            st.sidebar.write(f"ILThermo validated: {sum(1 for s in combinations if s.get('in_ilthermo', False))}")
            
            # Export results button
            if st.sidebar.button("Export Results"):
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
            
        except Exception as e:
            st.error(f"Error calculating properties: {str(e)}")
            raise e

# Main UI layout
st.set_page_config(page_title="Ionic Liquid Optimizer", layout="wide")

# Initialize property ranges in sidebar
get_user_ranges()

# Add calculate button to sidebar
if st.sidebar.button("Find Optimal Ionic Liquids"):
    display_results()