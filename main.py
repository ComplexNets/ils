import sys
import os
from contextlib import contextmanager
from combine_fragDensity import main as density_main
from combine_fragHeatCapacity import main as heat_capacity_main
from dataclasses import dataclass
from typing import Tuple, Dict, List
import numpy as np
from collections import defaultdict
from utils import generate_il_name, is_in_il_thermo
from core.pareto_optimizer import ParetoOptimizer

@dataclass
class PropertyCriteria:
    range: Tuple[float, float]
    importance: int
    unit: str

class MultiCriteriaOptimizer:
    def __init__(self):
        self.pareto_optimizer = ParetoOptimizer()
        self.properties = {
            'heat_capacity': PropertyCriteria(
                range=(200, 400),
                importance=3,
                unit="J/mol·K"
            ),
            'density': PropertyCriteria(
                range=(800, 1500),
                importance=3,
                unit="kg/m³"
            )
        }
        
    def set_criteria(self, property_name: str, range_vals: Tuple[float, float], importance: int):
        """Update criteria for a specific property"""
        # Convert importance (1-5) to weight (0-1)
        weight = importance / 5.0
        self.pareto_optimizer.set_constraint(
            property_name=property_name,
            min_val=range_vals[0],
            max_val=range_vals[1],
            weight=weight,
            is_strict=(importance >= 4)  # Make constraint strict if importance >= 4
        )
        self.properties[property_name].range = range_vals
        self.properties[property_name].importance = importance
    
    def calculate_score(self, heat_capacity: float, density: float) -> float:
        """Calculate score for a combination of properties"""
        solution = {
            'heat_capacity': heat_capacity,
            'density': density
        }
        # Use Pareto ranking for single solution
        ranked = self.pareto_optimizer.rank_solutions([solution])
        return ranked[0]['pareto_score'] if ranked else 0.0

def combine_fragments_and_calculate_properties():
    """Combine fragments and calculate properties with Pareto optimization"""
    optimizer = ParetoOptimizer()
    
    # Get initial fragments with properties
    fragments_data = get_filtered_fragments()
    
    # Screen fragments using relaxed Pareto optimization
    property_ranges = {
        'heat_capacity': (200, 400),
        'density': (800, 1500)
    }
    
    screened_fragments = {}
    for frag_type, frags in fragments_data.items():
        screened = optimizer.screen_fragments(frags, property_ranges)
        screened_fragments[frag_type] = screened
        
    # Generate valid combinations
    valid_combinations = []
    for cation in screened_fragments['Cation']:
        for anion in screened_fragments['Anion']:
            for alkyl in screened_fragments['Alkyl Chain']:
                # Generate combination name
                il_name = generate_il_name(cation, anion, alkyl)
                
                # Calculate properties
                heat_capacity = calculate_ionic_liquid_heat_capacity({
                    'cation': cation,
                    'anion': anion,
                    'alkyl': alkyl
                })
                
                density = calculate_ionic_liquid_density({
                    'cation': cation,
                    'anion': anion,
                    'alkyl': alkyl
                })
                
                if heat_capacity is not None and density is not None:
                    combination = {
                        'name': il_name,
                        'heat_capacity': heat_capacity,
                        'density': density,
                        'in_ilthermo': is_in_il_thermo(il_name)
                    }
                    valid_combinations.append(combination)
    
    # Get Pareto-optimal solutions and rankings
    pareto_front, ranked_solutions = optimizer.optimize_combinations(valid_combinations)
    
    return ranked_solutions  # Return ranked solutions for compatibility

def main(hc_importance, density_importance):
    # Example: Call the function that combines fragments and calculates properties
    valid_combinations = combine_fragments_and_calculate_properties()
    
    # Filter and score the combinations based on importance
    scored_combinations = []
    for combo in valid_combinations:
        score = (hc_importance * combo['heat_capacity_score'] +
                 density_importance * combo['density_score'])
        scored_combinations.append({
            'name': combo['name'],
            'score': score,
            'heat_capacity': combo['heat_capacity'],
            'density': combo['density'],
            'in_ilthermo': combo['in_ilthermo']
        })
    
    # Sort by score
    scored_combinations.sort(key=lambda x: x['score'], reverse=True)
    
    return scored_combinations

def process_ils(hc_importance, density_importance):
    # Example: Call the function that combines fragments and calculates properties
    valid_combinations = combine_fragments_and_calculate_properties()
    
    # Filter and score the combinations based on importance
    for combination in valid_combinations:
        # Calculate the score based on importance
        score = (combination['heat_capacity_score'] * hc_importance +
                 combination['density_score'] * density_importance)
        combination['score'] = score
        
        # Yield each result as it is calculated
        yield combination
