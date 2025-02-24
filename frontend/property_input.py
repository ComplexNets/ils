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
    optimize_higher: bool = True  # True if higher values are better, False if lower values are better

@dataclass
class ValidationCriteria:
    """Criteria for validating ionic liquid structures"""
    max_groups_per_chain: int = 6  # Default max groups allowed in an alkyl chain
    max_groups_per_type: int = 6   # Default max groups of a specific type allowed
    max_chain_length: int = 12     # Default max chain length
    min_chain_length: int = 1      # Default min chain length
    group_occurrence_upper: int = 3  # t1: Upper limit on group occurrences (Eq. 19)
    group_occurrence_lower: int = 0  # t2: Lower limit on group occurrences (Eq. 20)
    group_occurrence_exact: int = -1  # t3: Exact occurrence requirement, -1 means not enforced (Eq. 21)

class PropertyRanges:
    def __init__(self):
        self.properties = {
            'heat_capacity': PropertyCriteria(
                range=(100.0, 1200.0),  # Expanded range
                importance=3,
                unit="J/mol·K",
                optimize_higher=True  # Higher heat capacity is generally better
            ),
            'density': PropertyCriteria(
                range=(200.0, 2500.0),  # Expanded range for lower and higher density ILs
                importance=3,
                unit="kg/m³",
                optimize_higher=False  # Lower density might be preferred for some applications
            ),
            'toxicity': PropertyCriteria(
                range=(0.01, 200.0),  # Expanded for both highly toxic and less toxic ILs
                importance=3,
                unit="mM",
                optimize_higher=True  # Higher IC50 means less toxic
            ),
            'solubility': PropertyCriteria(
                range=(0.001, 5000.0),  # Expanded for near-zero and highly soluble ILs
                importance=3,
                unit="g/L",
                optimize_higher=True  # Higher solubility is generally better
            ),
            'hydrophobicity': PropertyCriteria(
                range=(-10.0, 15.0),  # Expanded for extremely hydrophilic and super hydrophobic ILs
                importance=3,
                unit="logP",
                optimize_higher=True
            ),
            'viscosity': PropertyCriteria(
                range=(0.1, 50000.0),  # Expanded for low viscosity and extremely viscous ILs
                importance=3,
                unit="cP",
                optimize_higher=False
            )
        }
        self.validation = ValidationCriteria()

    def update_property(self, property_name: str, range_values: tuple, weight: float = None, optimize_higher: bool = True):
        """Update a property's range, importance and optimization direction
        
        Args:
            property_name: Name of the property to update
            range_values: Tuple of (min, max) values
            weight: Optional weight/importance value (0-1)
            optimize_higher: True if higher values are better, False if lower values are better
        """
        if property_name not in self.properties:
            self.properties[property_name] = PropertyCriteria(
                range=range_values,
                importance=int(weight * 5) if weight is not None else 3,
                unit="",
                optimize_higher=optimize_higher
            )
        else:
            prop = self.properties[property_name]
            prop.range = range_values
            prop.optimize_higher = optimize_higher
            if weight is not None:
                prop.importance = int(weight * 5)

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
                
            # Toxicity Input
            print("\nToxicity (IC50):")
            min_toxicity = float(input(f"Minimum IC50 ({self.properties['toxicity'].unit}): "))
            max_toxicity = float(input(f"Maximum IC50 ({self.properties['toxicity'].unit}): "))
            toxicity_importance = int(input("Importance (1-5, where 5 is most important): "))
            
            # Validate importance
            if not 1 <= toxicity_importance <= 5:
                raise ValueError("Importance must be between 1 and 5")
            
            # Solubility Input
            print("\nSolubility:")
            min_solubility = float(input(f"Minimum solubility ({self.properties['solubility'].unit}): "))
            max_solubility = float(input(f"Maximum solubility ({self.properties['solubility'].unit}): "))
            solubility_importance = int(input("Importance (1-5, where 5 is most important): "))
            
            # Validate importance
            if not 1 <= solubility_importance <= 5:
                raise ValueError("Importance must be between 1 and 5")
            
            # Hydrophobicity Input
            print("\nHydrophobicity:")
            min_hydrophobicity = float(input(f"Minimum hydrophobicity ({self.properties['hydrophobicity'].unit}): "))
            max_hydrophobicity = float(input(f"Maximum hydrophobicity ({self.properties['hydrophobicity'].unit}): "))
            hydrophobicity_importance = int(input("Importance (1-5, where 5 is most important): "))
            
            # Validate importance
            if not 1 <= hydrophobicity_importance <= 5:
                raise ValueError("Importance must be between 1 and 5")
            
            # Viscosity Input
            print("\nViscosity:")
            min_viscosity = float(input(f"Minimum viscosity ({self.properties['viscosity'].unit}): "))
            max_viscosity = float(input(f"Maximum viscosity ({self.properties['viscosity'].unit}): "))
            viscosity_importance = int(input("Importance (1-5, where 5 is most important): "))
            
            # Validate importance
            if not 1 <= viscosity_importance <= 5:
                raise ValueError("Importance must be between 1 and 5")
            
            # Update the properties
            self.properties['heat_capacity'].range = (min_capacity, max_capacity)
            self.properties['heat_capacity'].importance = hc_importance
            self.properties['density'].range = (min_density, max_density)
            self.properties['density'].importance = density_importance
            self.properties['toxicity'].range = (min_toxicity, max_toxicity)
            self.properties['toxicity'].importance = toxicity_importance
            self.properties['solubility'].range = (min_solubility, max_solubility)
            self.properties['solubility'].importance = solubility_importance
            self.properties['hydrophobicity'].range = (min_hydrophobicity, max_hydrophobicity)
            self.properties['hydrophobicity'].importance = hydrophobicity_importance
            self.properties['viscosity'].range = (min_viscosity, max_viscosity)
            self.properties['viscosity'].importance = viscosity_importance
            
        except ValueError as e:
            print(f"Error: {e}")
            return None
        
        return self.properties

    def calculate_objective_score(self, heat_capacity: float, density: float, toxicity: float, solubility: float, hydrophobicity: float, viscosity: float) -> float:
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
            hc_prop.range[1],
            hc_prop.optimize_higher
        )
        scores.append(hc_score)
        weights.append(hc_weight)
        
        # Density Score
        density_prop = self.properties['density']
        density_weight = density_prop.importance / total_importance
        density_score = self._calculate_property_score(
            density, 
            density_prop.range[0], 
            density_prop.range[1],
            density_prop.optimize_higher
        )
        scores.append(density_score)
        weights.append(density_weight)
        
        # Toxicity Score
        toxicity_prop = self.properties['toxicity']
        toxicity_weight = toxicity_prop.importance / total_importance
        toxicity_score = self._calculate_property_score(
            toxicity, 
            toxicity_prop.range[0], 
            toxicity_prop.range[1],
            toxicity_prop.optimize_higher
        )
        scores.append(toxicity_score)
        weights.append(toxicity_weight)
        
        # Solubility Score
        solubility_prop = self.properties['solubility']
        solubility_weight = solubility_prop.importance / total_importance
        solubility_score = self._calculate_property_score(
            solubility, 
            solubility_prop.range[0], 
            solubility_prop.range[1],
            solubility_prop.optimize_higher
        )
        scores.append(solubility_score)
        weights.append(solubility_weight)
        
        # Hydrophobicity Score
        hydrophobicity_prop = self.properties['hydrophobicity']
        hydrophobicity_weight = hydrophobicity_prop.importance / total_importance
        hydrophobicity_score = self._calculate_property_score(
            hydrophobicity, 
            hydrophobicity_prop.range[0], 
            hydrophobicity_prop.range[1],
            hydrophobicity_prop.optimize_higher
        )
        scores.append(hydrophobicity_score)
        weights.append(hydrophobicity_weight)
        
        # Viscosity Score
        viscosity_prop = self.properties['viscosity']
        viscosity_weight = viscosity_prop.importance / total_importance
        viscosity_score = self._calculate_property_score(
            viscosity, 
            viscosity_prop.range[0], 
            viscosity_prop.range[1],
            viscosity_prop.optimize_higher
        )
        scores.append(viscosity_score)
        weights.append(viscosity_weight)
        
        # Calculate weighted average
        final_score = sum(s * w for s, w in zip(scores, weights))
        
        return final_score

    @staticmethod
    def _calculate_property_score(value: float, min_val: float, max_val: float, optimize_higher: bool) -> float:
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
        if optimize_higher:
            return np.exp(-distance / range_width)
        else:
            return np.exp(-(range_width - distance) / range_width)

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
        hc_min = st.number_input("Min", value=100.0, step=10.0, key="hc_min")
    with hc_col2:
        hc_max = st.number_input("Max", value=1200.0, step=10.0, key="hc_max")
    with hc_col3:
        hc_importance = st.slider("Importance", 1, 5, 3, key="hc_importance")
    
    # Density
    st.subheader("Density (kg/m³)")
    d_col1, d_col2, d_col3 = st.columns(3)
    with d_col1:
        d_min = st.number_input("Min", value=200.0, step=50.0, key="d_min")
    with d_col2:
        d_max = st.number_input("Max", value=2500.0, step=50.0, key="d_max")
    with d_col3:
        d_importance = st.slider("Importance", 1, 5, 3, key="d_importance")
    
    # Toxicity
    st.subheader("Toxicity (IC50, mM)")
    t_col1, t_col2, t_col3 = st.columns(3)
    with t_col1:
        t_min = st.number_input("Min", value=0.01, step=0.01, key="t_min")
    with t_col2:
        t_max = st.number_input("Max", value=200.0, step=1.0, key="t_max")
    with t_col3:
        t_importance = st.slider("Importance", 1, 5, 3, key="t_importance")
    
    # Solubility
    st.subheader("Solubility (g/L)")
    s_col1, s_col2, s_col3 = st.columns(3)
    with s_col1:
        s_min = st.number_input("Min", value=0.001, step=0.001, key="s_min")
    with s_col2:
        s_max = st.number_input("Max", value=5000.0, step=1.0, key="s_max")
    with s_col3:
        s_importance = st.slider("Importance", 1, 5, 3, key="s_importance")
    
    # Hydrophobicity
    st.subheader("Hydrophobicity (logP)")
    h_col1, h_col2, h_col3 = st.columns(3)
    with h_col1:
        h_min = st.number_input("Min", value=-10.0, step=0.1, key="h_min")
    with h_col2:
        h_max = st.number_input("Max", value=15.0, step=0.1, key="h_max")
    with h_col3:
        h_importance = st.slider("Importance", 1, 5, 3, key="h_importance")
    
    # Viscosity
    st.subheader("Viscosity (cP)")
    v_col1, v_col2, v_col3 = st.columns(3)
    with v_col1:
        v_min = st.number_input("Min", value=0.1, step=0.1, key="v_min")
    with v_col2:
        v_max = st.number_input("Max", value=50000.0, step=100.0, key="v_max")
    with v_col3:
        v_importance = st.slider("Importance", 1, 5, 3, key="v_importance")
    
    # Update optimizer constraints
    optimizer.set_constraint("heat_capacity", hc_min, hc_max, hc_importance/5.0)
    optimizer.set_constraint("density", d_min, d_max, d_importance/5.0)
    optimizer.set_constraint("toxicity", t_min, t_max, t_importance/5.0)
    optimizer.set_constraint("solubility", s_min, s_max, s_importance/5.0)
    optimizer.set_constraint("hydrophobicity", h_min, h_max, h_importance/5.0)
    optimizer.set_constraint("viscosity", v_min, v_max, v_importance/5.0)
    
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
                    st.write(f"Toxicity (IC50): {solution['toxicity']:.1f} mM")
                    st.write(f"Solubility: {solution['solubility']:.1f} g/L")
                    st.write(f"Hydrophobicity: {solution['hydrophobicity']:.1f} logP")
                    st.write(f"Viscosity: {solution['viscosity']:.1f} cP")
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
        test_toxicity = 50
        test_solubility = 500
        test_hydrophobicity = 2
        test_viscosity = 100
        score = property_ranges.calculate_objective_score(test_hc, test_density, test_toxicity, test_solubility, test_hydrophobicity, test_viscosity)
        print(f"Score for HC={test_hc}, density={test_density}, toxicity={test_toxicity}, solubility={test_solubility}, hydrophobicity={test_hydrophobicity}, viscosity={test_viscosity}: {score:.3f}")

if __name__ == "__main__":
    main()
