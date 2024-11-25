from typing import Dict, Tuple, List
import numpy as np
from itertools import product
from models.shortList_frag import fragments

# Constants for chemical rules
MAX_BOND_CAPACITY = 4  # Maximum total bond capacity
OCTET_TOTAL = 8  # Octet rule requirement

class MolecularValidator:
    """Validates ionic liquid combinations based on chemical rules"""
    
    def __init__(self):
        """Initialize validator with fragment properties"""
        self.fragment_properties = {
            'cation': {
                'Ethylimidazolium': {"valence": 1, "bond_capacity": 2},
                'Methylimidazolium': {"valence": 1, "bond_capacity": 2},
                'Propylimidazolium': {"valence": 1, "bond_capacity": 2},
                '1-Butylpyridinium': {"valence": 1, "bond_capacity": 2},
                '1-Ethylpyridinium': {"valence": 1, "bond_capacity": 2},
                '1-Propylpyridinium': {"valence": 1, "bond_capacity": 2}
            },
            'anion': {
                'Tetrafluoroborate': {"valence": -1, "bond_capacity": 1},
                'Bis(trifluoromethanesulfonyl)imide': {"valence": -1, "bond_capacity": 1},
                'Chloride': {"valence": -1, "bond_capacity": 1},
                'Bromide': {"valence": -1, "bond_capacity": 1},
                'Trifluoromethanesulfonate': {"valence": -1, "bond_capacity": 1},
                'Dicyanamide': {"valence": -1, "bond_capacity": 1}
            },
            'alkyl_chain': {
                'Methyl': {"valence": 0, "bond_capacity": 1},
                'Ethyl': {"valence": 0, "bond_capacity": 1},
                'Propyl': {"valence": 0, "bond_capacity": 1},
                'Butyl': {"valence": 0, "bond_capacity": 1},
                'Pentyl': {"valence": 0, "bond_capacity": 1},
                'Hexyl': {"valence": 0, "bond_capacity": 1},
                'Octyl': {"valence": 0, "bond_capacity": 1}
            }
        }
    
    def _check_fragment_exists(self, fragment: Dict) -> bool:
        """Check if fragment exists in validation rules"""
        try:
            fragment_name = fragment.get('name')  
            fragment_type = fragment.get('fragment_type')  
            if not fragment_name or not fragment_type:
                return False
            return fragment_name in self.fragment_properties.get(fragment_type.lower(), {})
        except Exception as e:
            print(f"Error checking fragment existence: {e}")
            return False
    
    def validate(self, cation: Dict, anion: Dict, alkyl: Dict) -> Tuple[bool, str]:
        """
        Validate a combination of fragments
        Returns:
            Tuple[bool, str]: (is_valid, message)
        """
        try:
            # Check if fragments exist in our properties
            if not self._check_fragment_exists(cation):
                return False, f"Cation {cation.get('name')} not found in validation rules"
            if not self._check_fragment_exists(anion):
                return False, f"Anion {anion.get('name')} not found in validation rules"
            if not self._check_fragment_exists(alkyl):
                return False, f"Alkyl chain {alkyl.get('name')} not found in validation rules"
            
            # Get fragment properties
            cat_props = self.fragment_properties['cation'][cation['name']]
            an_props = self.fragment_properties['anion'][anion['name']]
            alkyl_props = self.fragment_properties['alkyl_chain'][alkyl['name']]
            
            # Check valence balance
            total_valence = (cat_props['valence'] + 
                           an_props['valence'] + 
                           alkyl_props['valence'])
            
            if total_valence != 0:
                return False, f"Invalid valence balance: {total_valence}"
            
            # Check bond capacity
            total_bonds = (cat_props['bond_capacity'] + 
                         an_props['bond_capacity'] + 
                         alkyl_props['bond_capacity'])
            
            if total_bonds > MAX_BOND_CAPACITY:
                return False, f"Too many bonds: {total_bonds}"
            
            # All checks passed
            return True, "Valid combination"
            
        except Exception as e:
            return False, f"Validation error: {str(e)}"

def generate_valid_combinations(fragments_list):
    """Generate valid combinations from the fragments list"""
    # Separate fragments by type
    cations = [f for f in fragments_list if f.get('fragment_type').lower() == 'cation']
    anions = [f for f in fragments_list if f.get('fragment_type').lower() == 'anion']
    alkyls = [f for f in fragments_list if f.get('fragment_type').lower() == 'alkyl_chain']
    
    valid_combinations = []
    validator = MolecularValidator()
    for cation, anion, alkyl in product(cations, anions, alkyls):
        is_valid, message = validator.validate(cation, anion, alkyl)
        if is_valid:
            valid_combinations.append({
                'name': f"{cation['name']} {anion['name']} with {alkyl['name']}",
                'cation': cation,
                'anion': anion,
                'alkyl': alkyl,
            })
    
    return valid_combinations

if __name__ == "__main__":
    # Test the validation
    valid_combinations = generate_valid_combinations(fragments)
    print(f"Found {len(valid_combinations)} valid combinations:")
    for combo in valid_combinations:
        print(f"\nName: {combo['name']}")
