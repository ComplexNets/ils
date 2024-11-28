from typing import Dict, Tuple, List
import numpy as np
from itertools import product
import os
from utils.utils import get_fragment_properties

# Constants for chemical rules
MAX_BOND_CAPACITY = 4  # Maximum total bond capacity
OCTET_TOTAL = 8  # Octet rule requirement

class MolecularValidator:
    """Validates ionic liquid combinations based on chemical rules"""
    
    def __init__(self):
        """Initialize validator"""
        pass
    
    def _get_fragment_valence(self, fragment: Dict) -> int:
        """Get fragment valence based on fragment type"""
        fragment_type = fragment.get('fragment_type', '').lower()
        if fragment_type == 'cation':
            return 1
        elif fragment_type == 'anion':
            return -1
        elif fragment_type == 'alkyl_chain':
            return 0
        return 0
    
    def _get_fragment_bond_capacity(self, fragment: Dict) -> int:
        """Get fragment bond capacity based on properties"""
        props = get_fragment_properties(fragment['name'], fragment['fragment_type'])
        if not props:
            return 0
        return min(props.get('rotatable_bond_count', 0) + 1, MAX_BOND_CAPACITY)
    
    def validate(self, cation: Dict, anion: Dict, alkyl: Dict) -> Tuple[bool, str]:
        """
        Validate a combination of fragments
        Returns:
            Tuple[bool, str]: (is_valid, message)
        """
        try:
            # Check if we can get properties for all fragments
            fragments = [
                (cation, "Cation"),
                (anion, "Anion"),
                (alkyl, "Alkyl chain")
            ]
            
            for frag, frag_type in fragments:
                props = get_fragment_properties(frag['name'], frag['fragment_type'])
                if not props:
                    return False, f"{frag_type} {frag.get('name')} properties not found"
            
            # Get valences
            cat_valence = self._get_fragment_valence(cation)
            an_valence = self._get_fragment_valence(anion)
            alkyl_valence = self._get_fragment_valence(alkyl)
            
            # Check valence balance
            total_valence = cat_valence + an_valence + alkyl_valence
            if total_valence != 0:
                return False, f"Invalid valence balance: {total_valence}"
            
            # Get and check bond capacities
            cat_bonds = self._get_fragment_bond_capacity(cation)
            an_bonds = self._get_fragment_bond_capacity(anion)
            alkyl_bonds = self._get_fragment_bond_capacity(alkyl)
            
            total_bonds = cat_bonds + an_bonds + alkyl_bonds
            if total_bonds > MAX_BOND_CAPACITY:
                return False, f"Too many bonds: {total_bonds}"
            
            # All checks passed
            return True, "Valid combination"
            
        except Exception as e:
            return False, f"Validation error: {str(e)}"

def generate_valid_combinations(fragments_list):
    """Generate valid combinations from the fragments list"""
    # Separate fragments by type
    cations = [f for f in fragments_list if f.get('fragment_type', '').lower() == 'cation']
    anions = [f for f in fragments_list if f.get('fragment_type', '').lower() == 'anion']
    alkyl_chains = [f for f in fragments_list if f.get('fragment_type', '').lower() == 'alkyl_chain']
    
    validator = MolecularValidator()
    valid_combinations = []
    
    # Calculate total combinations for progress tracking
    total_combinations = len(cations) * len(anions) * len(alkyl_chains)
    checked_combinations = 0
    
    # Generate all possible combinations
    for cation in cations:
        for anion in anions:
            for alkyl in alkyl_chains:
                checked_combinations += 1
                
                # Update progress in streamlit if available
                if 'st' in globals():
                    progress = checked_combinations / total_combinations
                    st.progress(progress)
                
                is_valid, message = validator.validate(cation, anion, alkyl)
                if is_valid:
                    valid_combinations.append({
                        'cation': cation,
                        'anion': anion,
                        'alkyl_chain': alkyl
                    })
    
    return valid_combinations

if __name__ == "__main__":
    # Test the validation
    from models.shortList_frag import fragments
    valid_combinations = generate_valid_combinations(fragments)
