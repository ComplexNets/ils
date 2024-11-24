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
    is_strict: bool  # If True, solutions must be within range

class ParetoOptimizer:
    def __init__(self):
        self.properties = {
            'heat_capacity': PropertyConstraint(
                min_value=200,
                max_value=400,
                weight=0.5,
                is_strict=False
            ),
            'density': PropertyConstraint(
                min_value=800,
                max_value=1500,
                weight=0.5,
                is_strict=False
            )
        }
        
    def set_constraint(self, property_name: str, min_val: float, max_val: float, 
                      weight: float, is_strict: bool = False):
        """Update constraints for a specific property"""
        if property_name in self.properties:
            self.properties[property_name] = PropertyConstraint(
                min_value=min_val,
                max_value=max_val,
                weight=weight,
                is_strict=is_strict
            )

    def _normalize_property(self, value: float, prop_constraint: PropertyConstraint) -> float:
        """Normalize property value to [0,1] range"""
        range_width = prop_constraint.max_value - prop_constraint.min_value
        normalized = (value - prop_constraint.min_value) / range_width
        return max(0.0, min(1.0, normalized))

    def _calculate_dominance(self, solution1: Dict, solution2: Dict) -> bool:
        """
        Check if solution1 dominates solution2
        Returns True if solution1 dominates solution2, False otherwise
        """
        at_least_one_better = False
        for prop_name, constraint in self.properties.items():
            val1 = solution1.get(prop_name, 0)
            val2 = solution2.get(prop_name, 0)
            
            # Normalize values
            norm1 = self._normalize_property(val1, constraint)
            norm2 = self._normalize_property(val2, constraint)
            
            if norm1 < norm2:
                return False
            elif norm1 > norm2:
                at_least_one_better = True
                
        return at_least_one_better

    def get_pareto_front(self, solutions: List[Dict]) -> List[Dict]:
        """
        Identify Pareto-optimal solutions
        Args:
            solutions: List of dictionaries containing property values
        Returns:
            List of non-dominated solutions
        """
        pareto_front = []
        
        for solution in solutions:
            # Check strict constraints first
            if not self._meets_strict_constraints(solution):
                continue
                
            is_dominated = False
            
            # Compare with existing Pareto front
            for idx, pareto_solution in enumerate(pareto_front):
                if self._calculate_dominance(pareto_solution, solution):
                    is_dominated = True
                    break
                elif self._calculate_dominance(solution, pareto_solution):
                    pareto_front.pop(idx)
                    
            if not is_dominated:
                pareto_front.append(solution)
                
        return pareto_front

    def _meets_strict_constraints(self, solution: Dict) -> bool:
        """Check if solution meets all strict constraints"""
        for prop_name, constraint in self.properties.items():
            if not constraint.is_strict:
                continue
                
            value = solution.get(prop_name, 0)
            if not (constraint.min_value <= value <= constraint.max_value):
                return False
        return True

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
            if not self._meets_strict_constraints(solution):
                continue
                
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
