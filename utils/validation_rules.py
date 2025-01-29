from typing import Dict, Tuple, List, Optional
import numpy as np
from itertools import product
import os
from utils.utils import get_fragment_properties
from rdkit import Chem

# Constants for chemical rules
MAX_BOND_CAPACITY = 4  # Maximum total bond capacity

class MolecularValidator:
    """Validates ionic liquid combinations based on chemical rules"""
    
    def __init__(self, max_total_groups=6, min_total_groups=0, 
                 max_groups_per_type=6, min_groups_per_type=0,
                 max_alkyl_chain_length=12, max_alkyl_chains=4,
                 max_valence=4, min_valence=-4,
                 max_groups_per_chain=2, min_groups_per_chain=0,
                 max_group_occurrences=1, min_group_occurrences=0,
                 max_total_group_occurrences=2, min_total_group_occurrences=0):
        """Initialize validator with configurable parameters."""
        self.max_total_groups = max_total_groups
        self.min_total_groups = min_total_groups
        self.max_groups_per_type = max_groups_per_type
        self.min_groups_per_type = min_groups_per_type
        self.max_valence = max_valence
        self.min_valence = min_valence
        self.max_alkyl_chain_length = max_alkyl_chain_length
        self.max_alkyl_chains = max_alkyl_chains
        self.max_groups_per_chain = max_groups_per_chain
        self.min_groups_per_chain = min_groups_per_chain
        self.max_group_occurrences = max_group_occurrences
        self.min_group_occurrences = min_group_occurrences
        self.max_total_group_occurrences = max_total_group_occurrences
        self.min_total_group_occurrences = min_total_group_occurrences

    def _validate_required_types(self, fragments: List[Dict]) -> Tuple[bool, str]:
        """
        Validate that all required fragment types are present
        Args:
            fragments: List of fragment dictionaries
        Returns:
            Tuple[bool, str]: (is_valid, message)
        """
        try:
            # Check for required fragment types
            fragment_types = [f['fragment_type'].lower() for f in fragments if f]
            
            # Must have exactly one cation and one anion
            if fragment_types.count('cation') != 1:
                return False, "Must have exactly one cation"
            if fragment_types.count('anion') != 1:
                return False, "Must have exactly one anion"
                
            # Must have at least one alkyl chain
            if 'alkyl_chain' not in fragment_types:
                return False, "Must have at least one alkyl chain"
                
            # Functional groups are optional
            return True, "All required fragment types present"
            
        except Exception as e:
            return False, f"Error validating required types: {str(e)}"

    def validate_combination(self, fragments: List[Dict]) -> Tuple[bool, str]:
        """
        Validate a combination of fragments against all rules
        Args:
            fragments: List of fragment dictionaries
        Returns:
            Tuple of (is_valid, message)
        """
        try:
            if not fragments:
                print("No fragments provided")
                return False, "No fragments provided"
                
            print(f"\nValidating combination: {[f['name'] for f in fragments]}")
            
            # First check if we have required types
            is_valid, message = self._validate_required_types(fragments)
            if not is_valid:
                print(f"Failed required types check: {message}")
                return False, message
                
            # Then check alkyl chain count
            is_valid, message = self._alkyl_chain_count(fragments)
            if not is_valid:
                print(f"Failed alkyl chain count check: {message}")
                return False, message
                
            # Then check total groups
            is_valid, message = self._validate_total_groups(fragments)
            if not is_valid:
                print(f"Failed total groups check: {message}")
                return False, message
                
            # Then check groups per chain
            is_valid, message = self._validate_groups_per_chain(fragments)
            if not is_valid:
                print(f"Failed groups per chain check: {message}")
                return False, message
                
            # Check group occurrences if we have alkyl chains
            alkyl_chains = [f for f in fragments if f['fragment_type'].lower() == 'alkyl_chain']
            for chain in alkyl_chains:
                is_valid, message = self._validate_group_occurrences(chain)
                if not is_valid:
                    print(f"Failed group occurrences check: {message}")
                    return False, message
                    
            # Check cation-wide group occurrences
            is_valid, message = self._validate_cation_group_occurrences(fragments)
            if not is_valid:
                print(f"Failed cation group occurrences check: {message}")
                return False, message
                
            # Finally check modified octet rule
            is_valid, message = self._validate_modified_octet_rule(fragments)
            if not is_valid:
                print(f"Failed modified octet rule check: {message}")
                return False, message
                
            print("Combination passed all validation checks!")
            return True, "Valid combination"
            
        except Exception as e:
            print(f"Error validating combination: {str(e)}")
            return False, f"Error validating combination: {str(e)}"

    def validate(self, cation: Dict, anion: Dict, alkyl: Dict) -> Tuple[bool, str]:
        """
        Validate a combination of fragments
        Returns:
            Tuple[bool, str]: (is_valid, message)
        """
        try:
            # Create fragments list
            fragments = [cation, anion, alkyl]
            
            # First check if we have required types
            is_valid, message = self._validate_required_types(fragments)
            if not is_valid:
                print(f"Failed required types check: {message}")
                return False, message
                
            # Then check alkyl chain count
            is_valid, message = self._alkyl_chain_count(fragments)
            if not is_valid:
                print(f"Failed alkyl chain count check: {message}")
                return False, message
                
            # Then check total groups
            is_valid, message = self._validate_total_groups(fragments)
            if not is_valid:
                print(f"Failed total groups check: {message}")
                return False, message
                
            # Then check groups per chain
            is_valid, message = self._validate_groups_per_chain(fragments)
            if not is_valid:
                print(f"Failed groups per chain check: {message}")
                return False, message
                
            # Check group occurrences if we have alkyl chains
            alkyl_chains = [f for f in fragments if f['fragment_type'].lower() == 'alkyl_chain']
            for chain in alkyl_chains:
                is_valid, message = self._validate_group_occurrences(chain)
                if not is_valid:
                    print(f"Failed group occurrences check: {message}")
                    return False, message
                    
            # Check cation-wide group occurrences
            is_valid, message = self._validate_cation_group_occurrences(fragments)
            if not is_valid:
                print(f"Failed cation group occurrences check: {message}")
                return False, message
                
            # Finally check modified octet rule
            is_valid, message = self._validate_modified_octet_rule(fragments)
            if not is_valid:
                print(f"Failed modified octet rule check: {message}")
                return False, message
                
            print("Combination passed all validation checks!")
            return True, "Valid combination"
            
        except Exception as e:
            print(f"Error validating combination: {str(e)}")
            return False, f"Error validating combination: {str(e)}"

    def _get_paper_valence(self, fragment_type: str, fragment_name: str) -> Optional[int]:
        """Get the paper-defined valence for a fragment based on its type and name"""
        try:
            # Normalize case for consistent lookup
            fragment_type = fragment_type.lower().strip()
            fragment_name = fragment_name.lower().strip()
            
            # Cations have specific valences based on type
            if fragment_type == 'cation':
                cation_valences = {
                    'imidazolium': 2,    # Can have 2 side chains
                    'pyridinium': 1,     # Can have 1 side chain
                    'pyrrolidinium': 2,  # Can have 2 side chains
                    'piperidinium': 2,   # Can have 2 side chains
                    'morpholinium': 2,   # Can have 2 side chains
                    'ammonium': 4,       # Can have 4 side chains
                    'phosphonium': 4,    # Can have 4 side chains
                    # Common variations of the names
                    'phosphonium cation': 4,
                    'ammonium cation': 4,
                    'imidazolium cation': 2,
                    'pyridinium cation': 1,
                    'tetramethylammonium': 4,
                    'methylimidazolium': 2,
                    'methylpyridinium': 1,
                    'propylpyridinium': 1,
                    'hexylpyridinium': 1,
                    'butyl-methylimidazolium': 2,
                    'octyl-methylimidazolium': 2
                }
                
                # First try exact match
                valence = cation_valences.get(fragment_name)
                if valence is not None:
                    print(f"Found exact valence match for {fragment_name}: {valence}")
                    return valence
                    
                # Try matching base types
                base_types = ['imidazolium', 'pyridinium', 'ammonium', 'phosphonium']
                for base in base_types:
                    if base in fragment_name:
                        valence = cation_valences[base]
                        print(f"Matched {fragment_name} to base type {base} with valence {valence}")
                        return valence
                
                print(f"WARNING: Unknown cation type '{fragment_name}', defaulting to valence=1")
                return 1
            
            # Alkyl chains always have valence=1 (one connection point)
            elif fragment_type == 'alkyl_chain':
                return 1
            
            # Functional groups typically have valence=1 unless specified
            elif fragment_type == 'functional_group':
                special_group_valences = {
                    # Add any special functional groups that have valence != 1
                }
                return special_group_valences.get(fragment_name, 1)
            
            # Anions don't need valence
            elif fragment_type == 'anion':
                return 0
                
            print(f"WARNING: Unknown fragment type {fragment_type}")
            return None
            
        except Exception as e:
            print(f"Error getting paper valence: {str(e)}")
            return None

    def _alkyl_chain_count(self, fragments: List[Dict]) -> Tuple[bool, str]:
        """Validate that the number of alkyl chains matches the cation valence"""
        try:
            # Get cation from fragments
            cation = next((f for f in fragments if f['fragment_type'].lower() == 'cation'), None)
            if not cation:
                return False, "No cation found in fragments"
            
            print(f"DEBUG: Processing cation: {cation['name']}")
            print(f"DEBUG: Cation SMILES: {cation['smiles']}")
            
            # Get cation valence
            cation_valence = self._get_paper_valence(cation['fragment_type'], cation['name'])
            if cation_valence is None:
                return False, f"Could not calculate valence for cation: {cation['name']}"
                
            # Count alkyl chains
            alkyl_chains = [f for f in fragments if f['fragment_type'].lower() == 'alkyl_chain']
            alkyl_chain_count = len(alkyl_chains)
            
            # Validate chain count matches cation valence
            if alkyl_chain_count != cation_valence:
                return False, f"Number of alkyl chains ({alkyl_chain_count}) does not match cation valence ({cation_valence})"
            
            return True, f"Valid alkyl chain count ({alkyl_chain_count})"
            
        except Exception as e:
            return False, f"Error validating alkyl chain count: {str(e)}"

    def _validate_total_groups(self, fragments: List[Dict]) -> Tuple[bool, str]:
        """Validate total number of functional groups"""
        try:
            # Get functional groups
            functional_groups = [f for f in fragments if f['fragment_type'].lower() == 'functional_group']
            total_groups = len(functional_groups)
            
            # Validate total groups
            if total_groups > self.max_total_groups:
                return False, f"Total functional groups ({total_groups}) exceeds maximum ({self.max_total_groups})"
            if total_groups < self.min_total_groups:
                return False, f"Total functional groups ({total_groups}) is less than minimum ({self.min_total_groups})"
            
            return True, f"Valid total groups ({total_groups})"
            
        except Exception as e:
            return False, f"Error validating total groups: {str(e)}"

    def _validate_groups_per_chain(self, fragments: List[Dict]) -> Tuple[bool, str]:
        """Validate number of functional groups per alkyl chain"""
        try:
            # Get max groups per chain from validation criteria or use default
            max_groups_per_chain = self.max_groups_per_chain
            min_groups_per_chain = self.min_groups_per_chain
            
            # Get alkyl chains and functional groups
            alkyl_chains = [f for f in fragments if f['fragment_type'].lower() == 'alkyl_chain']
            functional_groups = [f for f in fragments if f['fragment_type'].lower() == 'functional_group']
            
            # For now, assume functional groups are evenly distributed
            # This will need to be updated when we have chain-group assignment information
            if len(alkyl_chains) > 0:
                groups_per_chain = len(functional_groups) / len(alkyl_chains)
                if groups_per_chain > max_groups_per_chain:
                    return False, f"Too many groups per chain ({groups_per_chain:.1f} > {max_groups_per_chain})"
                if groups_per_chain < min_groups_per_chain:
                    return False, f"Too few groups per chain ({groups_per_chain:.1f} < {min_groups_per_chain})"
            
            return True, "Valid groups per chain"
            
        except Exception as e:
            return False, f"Error validating groups per chain: {str(e)}"

    def _validate_group_occurrences(self, alkyl_chain: Dict) -> Tuple[bool, str]:
        """Validate occurrence constraints for functional groups on an alkyl chain"""
        try:
            # Get constraints from validation criteria or use defaults
            max_occurrences = self.max_group_occurrences
            min_occurrences = self.min_group_occurrences
            
            # For now, just check that we don't exceed max occurrences
            # This will need to be updated when we have actual group occurrence information
            return True, "Valid group occurrences"
            
        except Exception as e:
            return False, f"Error validating group occurrences: {str(e)}"

    def _validate_cation_group_occurrences(self, fragments: List[Dict]) -> Tuple[bool, str]:
        """Validate cation-wide occurrence constraints for functional groups"""
        try:
            # Get constraints from validation criteria or use defaults
            max_total_occurrences = self.max_total_group_occurrences
            min_total_occurrences = self.min_total_group_occurrences
            
            # Get functional groups
            functional_groups = [f for f in fragments if f['fragment_type'].lower() == 'functional_group']
            
            # Count occurrences of each group type
            group_counts = {}
            for group in functional_groups:
                group_name = group['name']
                group_counts[group_name] = group_counts.get(group_name, 0) + 1
                
                if group_counts[group_name] > max_total_occurrences:
                    return False, f"Group {group_name} occurs {group_counts[group_name]} times (max: {max_total_occurrences})"
                if group_counts[group_name] < min_total_occurrences:
                    return False, f"Group {group_name} occurs {group_counts[group_name]} times (min: {min_total_occurrences})"
            
            return True, "Valid group occurrences"
            
        except Exception as e:
            return False, f"Error validating group occurrences: {str(e)}"

    def _validate_modified_octet_rule(self, fragments: List[Dict]) -> Tuple[bool, str]:
        """Validate that the combination follows the modified octet rule"""
        try:
            # Calculate total valence
            total_valence = 0
            
            # Add cation valence
            cation = next((f for f in fragments if f['fragment_type'].lower() == 'cation'), None)
            if cation:
                cation_valence = self._get_paper_valence(cation['fragment_type'], cation['name'])
                if cation_valence is not None:
                    total_valence += cation_valence
            
            # Add alkyl chain valences
            alkyl_chains = [f for f in fragments if f['fragment_type'].lower() == 'alkyl_chain']
            for chain in alkyl_chains:
                chain_valence = self._get_paper_valence(chain['fragment_type'], chain['name'])
                if chain_valence is not None:
                    total_valence -= chain_valence
            
            # Add functional group valences
            functional_groups = [f for f in fragments if f['fragment_type'].lower() == 'functional_group']
            for group in functional_groups:
                group_valence = self._get_paper_valence(group['fragment_type'], group['name'])
                if group_valence is not None:
                    total_valence -= group_valence
            
            # Validate total valence
            if total_valence > self.max_valence:
                return False, f"Invalid total valence: {total_valence} (max: {self.max_valence})"
            if total_valence < self.min_valence:
                return False, f"Invalid total valence: {total_valence} (min: {self.min_valence})"
            
            return True, "Valid valence"
            
        except Exception as e:
            return False, f"Error validating modified octet rule: {str(e)}"
