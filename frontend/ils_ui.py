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
    optimizer.set_constraint(
        "density", density_min, density_max,
        weight=density_importance/5.0,
        is_strict=density_strict
    )
    optimizer.set_constraint(
        "heat_capacity", cp_min, cp_max,
        weight=cp_importance/5.0,
        is_strict=cp_strict
    )
    optimizer.set_constraint(
        "toxicity", toxicity_min, toxicity_max,
        weight=toxicity_importance/5.0,
        is_strict=toxicity_strict
    )
    
    return prop_ranges

def calculate_properties():
    """Get valid combinations and calculate their properties"""
    # Get user-defined ranges
    prop_ranges = st.session_state.property_ranges
    
    # First get valid combinations from validation rules
    combinations = combine_fragments()
    if not combinations:
        return []
        
    # Calculate properties for each combination
    valid_combinations = []
    for combo in combinations:
        # Calculate density, heat capacity, and toxicity using detailed methods
        density = calculate_ionic_liquid_density(combo)
        heat_capacity = calculate_ionic_liquid_heat_capacity(combo)
        toxicity_result = calculate_ionic_liquid_toxicity(combo)
        
        # Check if properties are within user-defined ranges
        if density is not None and heat_capacity is not None and toxicity_result is not None:
            density_range = prop_ranges.properties['density'].range
            cp_range = prop_ranges.properties['heat_capacity'].range
            toxicity_range = prop_ranges.properties['toxicity'].range
            
            toxicity_ic50 = toxicity_result['ic50_mm']
            
            if (density_range[0] <= density <= density_range[1] and
                cp_range[0] <= heat_capacity <= cp_range[1] and
                toxicity_range[0] <= toxicity_ic50 <= toxicity_range[1]):
                
                combo.update({
                    'density': density,
                    'heat_capacity': heat_capacity,
                    'toxicity': toxicity_ic50,
                    'toxicity_components': toxicity_result['components']
                })
                valid_combinations.append(combo)
    
    return valid_combinations

def plot_pareto_front(combinations, pareto_front, prop_ranges):
    """Create Pareto front visualization"""
    fig = go.Figure()
    
    # All solutions
    fig.add_trace(go.Scatter(
        x=[sol['heat_capacity'] for sol in combinations],
        y=[sol['density'] for sol in combinations],
        mode='markers',
        name='All Solutions',
        marker=dict(color='blue', size=8, opacity=0.5)
    ))
    
    # Pareto front
    fig.add_trace(go.Scatter(
        x=[sol['heat_capacity'] for sol in pareto_front],
        y=[sol['density'] for sol in pareto_front],
        mode='markers+lines',
        name='Pareto Front',
        marker=dict(color='red', size=10)
    ))
    
    # Add property ranges
    fig.add_shape(
        type="rect",
        x0=prop_ranges.properties['heat_capacity'].range[0],
        x1=prop_ranges.properties['heat_capacity'].range[1],
        y0=prop_ranges.properties['density'].range[0],
        y1=prop_ranges.properties['density'].range[1],
        line=dict(color="green", width=2),
        fillcolor="green",
        opacity=0.1
    )
    
    fig.update_layout(
        title="Pareto Front of Ionic Liquid Properties",
        xaxis_title="Heat Capacity (J/mol·K)",
        yaxis_title="Density (kg/m³)",
        showlegend=True
    )
    
    return fig

def plot_property_correlation(combinations):
    """Create correlation plot between properties"""
    fig = go.Figure()
    
    # Create scatter plot with density coloring
    fig.add_trace(go.Scatter(
        x=[sol['heat_capacity'] for sol in combinations],
        y=[sol['density'] for sol in combinations],
        mode='markers',
        marker=dict(
            size=8,
            color=[sol.get('pareto_score', 0) for sol in combinations],
            colorscale='Viridis',
            showscale=True,
            colorbar=dict(title="Pareto Score")
        ),
        text=[sol['name'] for sol in combinations],
        hovertemplate="<b>%{text}</b><br>" +
                     "Heat Capacity: %{x:.1f} J/mol·K<br>" +
                     "Density: %{y:.1f} kg/m³<br>" +
                     "<extra></extra>"
    ))
    
    fig.update_layout(
        title="Property Correlation with Pareto Scores",
        xaxis_title="Heat Capacity (J/mol·K)",
        yaxis_title="Density (kg/m³)",
        showlegend=False
    )
    
    return fig

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
            combinations = calculate_properties()
            if not combinations:
                st.warning("No valid ionic liquid combinations found.")
                return
            
            # Get Pareto front and rankings
            optimizer = st.session_state.optimizer
            pareto_front, ranked_solutions = optimizer.optimize_combinations(combinations)
            
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
                st.plotly_chart(fig_pareto, use_container_width=True)
            
            with tab2:
                st.subheader("Property Correlation Analysis")
                fig_corr = plot_property_correlation(combinations)
                st.plotly_chart(fig_corr, use_container_width=True)
            
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
                    sol for sol in ranked_solutions 
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
                for sol in ranked_solutions:
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