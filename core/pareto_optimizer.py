"""
Pareto optimization module for ionic liquids.
Handles multi-objective optimization using Pareto fronts.
"""

from typing import Dict, List, Tuple, Optional
import numpy as np

class ParetoOptimizer:
    def __init__(self):
        """Initialize the Pareto optimizer with default constraints"""
        self.constraints = {}
        self.objectives = {}
    
    def set_constraint(self, property_name: str, min_val: float, max_val: float, 
                      weight: float = 0.5, optimize_higher: bool = True):
        """
        Set a constraint for a property
        
        Args:
            property_name: Name of the property
            min_val: Minimum acceptable value
            max_val: Maximum acceptable value
            weight: Importance weight (0-1 scale, where 0 means ignore this property)
            optimize_higher: If True, higher values are better. If False, lower values are better.
        """
        self.constraints[property_name] = {
            'min': min_val,
            'max': max_val,
            'weight': weight,
            'optimize_higher': optimize_higher
        }
    
    def normalize_value(self, value: float, min_val: float, max_val: float) -> float:
        """Normalize a value to 0-1 scale"""
        if max_val == min_val:
            return 1.0 if value == max_val else 0.0
        return (value - min_val) / (max_val - min_val)
    
    def calculate_objective_score(self, properties: Dict[str, float]) -> float:
        """
        Calculate the objective score for a set of properties
        
        Args:
            properties: Dictionary of property values
            
        Returns:
            float: Weighted sum of normalized property values
        """
        total_score = 0.0
        total_weight = 0.0
        
        for prop_name, constraint in self.constraints.items():
            if prop_name not in properties:
                continue
                
            value = properties[prop_name]
            weight = constraint['weight']
            
            # Skip if weight is 0 (property is ignored)
            if weight == 0:
                continue
                
            # Normalize the value
            norm_value = self.normalize_value(value, constraint['min'], constraint['max'])
            
            # Invert if we want to minimize
            if not constraint['optimize_higher']:
                norm_value = 1.0 - norm_value
            
            # Add to weighted sum
            total_score += norm_value * weight
            total_weight += weight
        
        # Return normalized score
        return total_score / total_weight if total_weight > 0 else 0.0
    
    def is_dominated(self, a: Dict[str, float], b: Dict[str, float]) -> bool:
        """
        Check if solution a is dominated by solution b
        
        Args:
            a: First solution's properties
            b: Second solution's properties
            
        Returns:
            bool: True if a is dominated by b
        """
        at_least_one_better = False
        
        for prop_name, constraint in self.constraints.items():
            if constraint['weight'] == 0:  # Skip properties with zero weight
                continue
                
            a_val = a.get(prop_name, 0)
            b_val = b.get(prop_name, 0)
            
            if constraint['optimize_higher']:
                if b_val < a_val:
                    return False
                if b_val > a_val:
                    at_least_one_better = True
            else:
                if b_val > a_val:
                    return False
                if b_val < a_val:
                    at_least_one_better = True
        
        return at_least_one_better
    
    def find_pareto_front(self, solutions: List[Dict[str, float]]) -> List[Dict[str, float]]:
        """
        Find the Pareto front from a list of solutions
        
        Args:
            solutions: List of dictionaries containing property values
            
        Returns:
            List[Dict]: Non-dominated solutions
        """
        pareto_front = []
        
        for solution in solutions:
            dominated = False
            
            # Compare with solutions already in Pareto front
            for pareto_solution in pareto_front[:]:
                if self.is_dominated(solution, pareto_solution):
                    dominated = True
                    break
                elif self.is_dominated(pareto_solution, solution):
                    pareto_front.remove(pareto_solution)
            
            if not dominated:
                pareto_front.append(solution)
        
        return pareto_front

    def get_pareto_front(self, solutions: List[Dict[str, float]]) -> List[Dict[str, float]]:
        """
        Alias for find_pareto_front for backward compatibility
        
        Args:
            solutions: List of dictionaries containing property values
            
        Returns:
            List[Dict]: Non-dominated solutions
        """
        return self.find_pareto_front(solutions)

    def rank_solutions(self, solutions: List[Dict[str, float]]) -> List[Dict[str, float]]:
        """
        Rank solutions using weighted sum of normalized properties
        
        Args:
            solutions: List of dictionaries containing property values
            
        Returns:
            List[Dict]: Solutions sorted by rank (best to worst)
        """
        ranked_solutions = []
        
        for solution in solutions:
            total_score = 0.0
            total_weight = 0.0
            
            for prop_name, constraint in self.constraints.items():
                value = solution.get(prop_name, 0)
                weight = constraint['weight']
                
                # Skip if weight is 0 (property is ignored)
                if weight == 0:
                    continue
                    
                # Normalize the value
                norm_value = self.normalize_value(value, constraint['min'], constraint['max'])
                
                # Invert if we want to minimize
                if not constraint['optimize_higher']:
                    norm_value = 1.0 - norm_value
                
                # Add to weighted sum
                total_score += norm_value * weight
                total_weight += weight
            
            # Return normalized score
            solution['pareto_score'] = total_score / total_weight if total_weight > 0 else 0.0
            ranked_solutions.append(solution)
        
        return sorted(ranked_solutions, key=lambda x: x['pareto_score'], reverse=True)

    def screen_fragments(self, fragments: List[Dict[str, float]], property_ranges: Dict[str, Tuple[float, float]]) -> List[Dict[str, float]]:
        """
        Screen fragments using relaxed Pareto optimization
        
        Args:
            fragments: List of fragment dictionaries with properties
            property_ranges: Dictionary of target property ranges
            
        Returns:
            List[Dict]: Fragments that pass screening
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

    def optimize_combinations(self, combinations: List[Dict[str, float]]) -> Tuple[List[Dict[str, float]], List[Dict[str, float]]]:
        """
        Optimize ionic liquid combinations using Pareto optimization
        
        Args:
            combinations: List of ionic liquid combinations with properties
            
        Returns:
            Tuple of (pareto_optimal_solutions, ranked_solutions)
        """
        # Get Pareto-optimal solutions
        pareto_front = self.find_pareto_front(combinations)
        
        # Rank all solutions
        ranked_solutions = self.rank_solutions(combinations)
        
        return pareto_front, ranked_solutions
