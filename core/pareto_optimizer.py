"""
Pareto optimization module for multi-objective ionic liquid property optimization.
Handles both fragment screening and final ionic liquid optimization.
"""
from typing import List, Dict, Tuple, Optional
import numpy as np
from dataclasses import dataclass

@dataclass
class PropertyConstraint:
    min_value: float
    max_value: float
    weight: float  # Importance weight (0-1)
    optimize_higher: bool = True  # If True, higher values are better

class ParetoOptimizer:
    def __init__(self):
        # Initialize with empty properties dict - constraints will be set via set_constraint
        self.properties = {}
        
    def set_constraint(self, property_name: str, min_val: float, max_val: float, 
                      weight: float, optimize_higher: bool = True):
        """Update constraints for a specific property"""
        self.properties[property_name] = PropertyConstraint(
            min_value=min_val,
            max_value=max_val,
            weight=weight,
            optimize_higher=optimize_higher
        )

    def _normalize_property(self, value: float, prop_constraint: PropertyConstraint) -> float:
        """Normalize property value to [0,1] range, considering optimization direction"""
        range_width = prop_constraint.max_value - prop_constraint.min_value
        if range_width == 0:
            return 0.0
            
        # First normalize to [0,1] range
        normalized = (value - prop_constraint.min_value) / range_width
        normalized = max(0.0, min(1.0, normalized))
        
        # If we want to minimize this property, invert the score
        if not prop_constraint.optimize_higher:
            normalized = 1.0 - normalized
            
        return normalized

    def _calculate_dominance(self, solution1: Dict, solution2: Dict) -> bool:
        """
        Check if solution1 dominates solution2 using weighted comparison
        Returns True if solution1 dominates solution2, False otherwise
        """
        total_weight = 0
        weighted_score1 = 0
        weighted_score2 = 0
        
        for prop_name, constraint in self.properties.items():
            val1 = solution1.get(prop_name, 0)
            val2 = solution2.get(prop_name, 0)
            
            # Normalize values considering optimization direction
            norm1 = self._normalize_property(val1, constraint)
            norm2 = self._normalize_property(val2, constraint)
            
            # Add weighted scores
            weighted_score1 += norm1 * constraint.weight
            weighted_score2 += norm2 * constraint.weight
            total_weight += constraint.weight
            
        # Compare weighted averages
        if total_weight > 0:
            return (weighted_score1 / total_weight) > (weighted_score2 / total_weight)
        return False

    def get_pareto_front(self, solutions: List[Dict]) -> List[Dict]:
        """
        Identify Pareto-optimal solutions using weighted scoring
        Args:
            solutions: List of dictionaries containing property values
        Returns:
            List of solutions sorted by weighted score
        """
        if not solutions:
            return []
            
        # Calculate weighted scores for all solutions
        scored_solutions = []
        for solution in solutions:
            total_score = 0
            total_weight = 0
            
            for prop_name, constraint in self.properties.items():
                if prop_name in solution:
                    value = solution[prop_name]
                    norm_value = self._normalize_property(value, constraint)
                    total_score += norm_value * constraint.weight
                    total_weight += constraint.weight
                    
            if total_weight > 0:
                weighted_score = total_score / total_weight
                scored_solutions.append((weighted_score, solution))
                
        # Sort by weighted score
        scored_solutions.sort(reverse=True, key=lambda x: x[0])
        
        # Return sorted solutions
        return [solution for score, solution in scored_solutions]

    def rank_solutions(self, solutions: List[Dict]) -> List[Dict]:
        """
        Rank solutions using weighted sum of normalized properties
        Args:
            solutions: List of dictionaries containing property values
        Returns:
            List of solutions sorted by rank (best to worst)
        """
        ranked_solutions = []
        
        for solution in solutions:
            total_score = 0
            total_weight = 0
            
            for prop_name, constraint in self.properties.items():
                value = solution.get(prop_name, 0)
                normalized = self._normalize_property(value, constraint)
                total_score += normalized * constraint.weight
                total_weight += constraint.weight
                
            if total_weight > 0:
                solution['pareto_score'] = total_score / total_weight
                ranked_solutions.append(solution)
                
        return sorted(ranked_solutions, key=lambda x: x['pareto_score'], reverse=True)

    def screen_fragments(self, fragments: List[Dict], property_ranges: Dict[str, Tuple[float, float]]) -> List[Dict]:
        """
        Screen fragments using relaxed Pareto optimization
        Args:
            fragments: List of fragment dictionaries with properties
            property_ranges: Dictionary of target property ranges
        Returns:
            List of fragments that pass screening
        """
        # Use wider ranges for fragment screening (typically 60-70% of final range)
        screened_fragments = []
        
        for fragment in fragments:
            meets_criteria = True
            for prop_name, (min_val, max_val) in property_ranges.items():
                if prop_name not in fragment:
                    continue
                    
                # Use wider ranges for screening (±30%)
                range_width = max_val - min_val
                screening_min = min_val - (0.3 * range_width)
                screening_max = max_val + (0.3 * range_width)
                
                if not (screening_min <= fragment[prop_name] <= screening_max):
                    meets_criteria = False
                    break
                    
            if meets_criteria:
                screened_fragments.append(fragment)
                
        return screened_fragments

    def optimize_combinations(self, combinations: List[Dict]) -> Tuple[List[Dict], List[Dict]]:
        """
        Optimize ionic liquid combinations using Pareto optimization
        Args:
            combinations: List of ionic liquid combinations with properties
        Returns:
            Tuple of (pareto_optimal_solutions, ranked_solutions)
        """
        # Get Pareto-optimal solutions
        pareto_front = self.get_pareto_front(combinations)
        
        # Rank all solutions
        ranked_solutions = self.rank_solutions(combinations)
        
        return pareto_front, ranked_solutions
