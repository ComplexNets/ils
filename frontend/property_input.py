from dataclasses import dataclass
from typing import Tuple, Dict
import numpy as np
import plotly.graph_objects as go
from core.pareto_optimizer import ParetoOptimizer

@dataclass
class PropertyCriteria:
    range: Tuple[float, float]
    importance: int  # 1-5 scale
    unit: str

class PropertyRanges:
    def __init__(self):
        self.properties = {
            'heat_capacity': PropertyCriteria(
                range=(0, 8000),
                importance=3,
                unit="J/mol·K"
            ),
            'density': PropertyCriteria(
                range=(0, 8000),
                importance=3,
                unit="kg/m³"
            )
        }

    def get_user_input(self):
        """Get property ranges and importance weights from user input"""
        try:
            print("\nEnter desired property ranges and importance (1-5 scale):")
            
            # Heat Capacity Input
            print("\nHeat Capacity:")
            min_capacity = float(input(f"Minimum heat capacity ({self.properties['heat_capacity'].unit}): "))
            max_capacity = float(input(f"Maximum heat capacity ({self.properties['heat_capacity'].unit}): "))
            hc_importance = int(input("Importance (1-5, where 5 is most important): "))
            
            # Validate importance
            if not 1 <= hc_importance <= 5:
                raise ValueError("Importance must be between 1 and 5")
            
            # Density Input
            print("\nDensity:")
            min_density = float(input(f"Minimum density ({self.properties['density'].unit}): "))
            max_density = float(input(f"Maximum density ({self.properties['density'].unit}): "))
            density_importance = int(input("Importance (1-5, where 5 is most important): "))
            
            # Validate importance
            if not 1 <= density_importance <= 5:
                raise ValueError("Importance must be between 1 and 5")
            
            # Update the properties
            self.properties['heat_capacity'].range = (min_capacity, max_capacity)
            self.properties['heat_capacity'].importance = hc_importance
            self.properties['density'].range = (min_density, max_density)
            self.properties['density'].importance = density_importance
            
        except ValueError as e:
            print(f"Error: {e}")
            return None
        
        return self.properties

    def calculate_objective_score(self, heat_capacity: float, density: float) -> float:
        """
        Calculate weighted objective score for a given ionic liquid
        Returns score between 0 and 1, where 1 is best match
        """
        scores = []
        weights = []
        
        # Normalize importance weights
        total_importance = sum(prop.importance for prop in self.properties.values())
        
        # Heat Capacity Score
        hc_prop = self.properties['heat_capacity']
        hc_weight = hc_prop.importance / total_importance
        hc_score = self._calculate_property_score(
            heat_capacity, 
            hc_prop.range[0], 
            hc_prop.range[1]
        )
        scores.append(hc_score)
        weights.append(hc_weight)
        
        # Density Score
        density_prop = self.properties['density']
        density_weight = density_prop.importance / total_importance
        density_score = self._calculate_property_score(
            density, 
            density_prop.range[0], 
            density_prop.range[1]
        )
        scores.append(density_score)
        weights.append(density_weight)
        
        # Calculate weighted average
        final_score = sum(s * w for s, w in zip(scores, weights))
        
        return final_score

    @staticmethod
    def _calculate_property_score(value: float, min_val: float, max_val: float) -> float:
        """
        Calculate how well a property value fits within the desired range
        Returns score between 0 and 1
        """
        if min_val <= value <= max_val:
            # If within range, give perfect score
            return 1.0
        
        # Calculate distance from nearest boundary
        distance = min(abs(value - min_val), abs(value - max_val))
        # Use exponential decay for scores outside range
        range_width = max_val - min_val
        return np.exp(-distance / range_width)

def display_property_inputs():
    import streamlit as st
    
    st.title("Ionic Liquid Property Optimization")
    
    # Initialize optimizer
    optimizer = ParetoOptimizer()
    
    # Property inputs
    st.header("Property Constraints")
    
    # Heat Capacity
    st.subheader("Heat Capacity (J/mol·K)")
    hc_col1, hc_col2, hc_col3 = st.columns(3)
    with hc_col1:
        hc_min = st.number_input("Min", value=200.0, step=10.0, key="hc_min")
    with hc_col2:
        hc_max = st.number_input("Max", value=400.0, step=10.0, key="hc_max")
    with hc_col3:
        hc_importance = st.slider("Importance", 1, 5, 3, key="hc_importance")
    hc_strict = st.checkbox("Strict constraint", key="hc_strict")
    
    # Density
    st.subheader("Density (kg/m³)")
    d_col1, d_col2, d_col3 = st.columns(3)
    with d_col1:
        d_min = st.number_input("Min", value=800.0, step=50.0, key="d_min")
    with d_col2:
        d_max = st.number_input("Max", value=1500.0, step=50.0, key="d_max")
    with d_col3:
        d_importance = st.slider("Importance", 1, 5, 3, key="d_importance")
    d_strict = st.checkbox("Strict constraint", key="d_strict")
    
    # Update optimizer constraints
    optimizer.set_constraint("heat_capacity", hc_min, hc_max, hc_importance/5.0, hc_strict)
    optimizer.set_constraint("density", d_min, d_max, d_importance/5.0, d_strict)
    
    # Run optimization
    if st.button("Find Optimal Ionic Liquids"):
        with st.spinner("Optimizing..."):
            # Get combinations with properties
            combinations = combine_fragments_and_calculate_properties()
            
            # Get Pareto front and rankings
            pareto_front, ranked_solutions = optimizer.optimize_combinations(combinations)
            
            # Display results
            st.header("Optimization Results")
            
            # Plot Pareto front
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
            fig.add_shape(type="rect",
                x0=hc_min, y0=d_min, x1=hc_max, y1=d_max,
                line=dict(color="green", width=2),
                fillcolor="green", opacity=0.1
            )
            
            fig.update_layout(
                title="Pareto Front of Ionic Liquid Properties",
                xaxis_title="Heat Capacity (J/mol·K)",
                yaxis_title="Density (kg/m³)",
                showlegend=True
            )
            
            st.plotly_chart(fig)
            
            # Display top solutions
            st.subheader("Top Ionic Liquid Combinations")
            for i, solution in enumerate(ranked_solutions[:10]):
                with st.expander(f"{i+1}. {solution['name']} (Score: {solution['pareto_score']:.3f})"):
                    st.write(f"Heat Capacity: {solution['heat_capacity']:.1f} J/mol·K")
                    st.write(f"Density: {solution['density']:.1f} kg/m³")
                    if solution['in_ilthermo']:
                        st.write("✓ Found in ILThermo database")

def main():
    property_ranges = PropertyRanges()
    properties = property_ranges.get_user_input()
    
    if properties:
        print("\nSelected property ranges and importance:")
        for name, prop in properties.items():
            print(f"\n{name.replace('_', ' ').title()}:")
            print(f"Range: {prop.range[0]} to {prop.range[1]} {prop.unit}")
            print(f"Importance: {prop.importance}/5")
        
        # Example calculation
        print("\nExample objective score calculation:")
        test_hc = 250
        test_density = 1000
        score = property_ranges.calculate_objective_score(test_hc, test_density)
        print(f"Score for HC={test_hc}, density={test_density}: {score:.3f}")

if __name__ == "__main__":
    main()
