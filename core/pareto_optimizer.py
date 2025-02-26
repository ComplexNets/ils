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
        """
        Normalize a value to 0-1 scale with special handling for zero values and edge cases
        """
        # Handle case where min and max are the same
        if max_val == min_val:
            # If the target value is exactly what we want, return 1.0
            if value == max_val:
                return 1.0
            # For zero case, treat values close to zero as valid
            if max_val == 0 and abs(value) < 1e-10:
                return 1.0
            return 0.0
            
        # Normal case - normalize between min and max
        normalized = (value - min_val) / (max_val - min_val)
        # Clamp between 0 and 1 to handle floating point errors
        return max(0.0, min(1.0, normalized))
    
    def calculate_objective_score(self, properties: Dict[str, float]) -> float:
        """
        Calculate the objective score for a set of properties with improved handling of edge cases
        
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
                
            # Special handling for zero target values
            if constraint['min'] == 0 and constraint['max'] == 0:
                # If we want exact zero and get very close to zero
                if abs(value) < 1e-10:
                    norm_value = 1.0
                else:
                    norm_value = 0.0
            else:
                # Normal normalization
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
        Check if solution a is dominated by solution b with improved handling of edge cases
        
        Args:
            a: First solution's properties
            b: Second solution's properties
            
        Returns:
            bool: True if a is dominated by b
        """
        at_least_one_better = False
        epsilon = 1e-10  # Small tolerance for floating point comparisons
        
        for prop_name, constraint in self.constraints.items():
            if constraint['weight'] == 0:  # Skip properties with zero weight
                continue
                
            # Get property values, defaulting to the constraint bounds if not present
            a_val = a.get(prop_name, constraint.get('min', 0) if not constraint['optimize_higher'] else constraint.get('max', 0))
            b_val = b.get(prop_name, constraint.get('min', 0) if not constraint['optimize_higher'] else constraint.get('max', 0))
            
            # Special handling for zero values
            if abs(a_val) < epsilon and abs(b_val) < epsilon:
                continue  # Consider them equal if both are effectively zero
                
            if constraint['optimize_higher']:
                if b_val < a_val - epsilon:  # Use epsilon for comparison
                    return False
                if b_val > a_val + epsilon:  # Use epsilon for comparison
                    at_least_one_better = True
            else:
                if b_val > a_val + epsilon:  # Use epsilon for comparison
                    return False
                if b_val < a_val - epsilon:  # Use epsilon for comparison
                    at_least_one_better = True
        
        return at_least_one_better
    
    def find_pareto_front(self, solutions: List[Dict[str, float]]) -> List[Dict[str, float]]:
        """
        Find the Pareto front from a list of solutions with improved handling for edge cases
        
        Args:
            solutions: List of dictionaries containing property values
            
        Returns:
            List[Dict]: Non-dominated solutions
        """
        if not solutions:
            return []
            
        # First filter solutions that are within constraints
        valid_solutions = []
        for solution in solutions:
            is_valid = True
            for prop_name, constraint in self.constraints.items():
                if prop_name not in solution:
                    continue
                    
                value = solution[prop_name]
                min_val = constraint['min']
                max_val = constraint['max']
                
                # Special handling for zero constraints
                if min_val == 0 and max_val == 0:
                    if abs(value) > 1e-10:  # Use small tolerance
                        is_valid = False
                        break
                else:
                    # Add small tolerance to bounds for floating point comparison
                    if not (min_val - 1e-10 <= value <= max_val + 1e-10):
                        is_valid = False
                        break
            
            if is_valid:
                valid_solutions.append(solution)
        
        if not valid_solutions:
            return []
            
        pareto_front = []
        
        for solution in valid_solutions:
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
