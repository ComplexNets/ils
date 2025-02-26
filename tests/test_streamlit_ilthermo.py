"""
Test script to check if the in_ilthermo property is being correctly set in the Streamlit app
"""
import sys
import os

# Get absolute path to project root
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
print(f"Project root: {project_root}")

# Add project root to Python path if not already there
if project_root not in sys.path:
    sys.path.insert(0, project_root)
    print(f"Added {project_root} to Python path")

from utils.utils import is_in_il_thermo, generate_il_name
from core.combine_fragments import validate_combination
from utils.validation_rules import MolecularValidator

def test_validate_combination():
    """Test if the in_ilthermo property is being correctly set in validate_combination"""
    
    # Create a simple combination with a known IL
    cation = {
        "name": "1-Ethyl-3-methylimidazolium", 
        "fragment_type": "cation", 
        "smiles": "CC[n+]1cc[nH]c1C",
        "charge": 1,
        "molecular_weight": 125.17
    }
    anion = {
        "name": "Tetrafluoroborate", 
        "fragment_type": "anion", 
        "smiles": "F[B-](F)(F)F",
        "charge": -1,
        "molecular_weight": 86.80
    }
    alkyl = {
        "name": "Ethyl", 
        "fragment_type": "alkyl_chain", 
        "smiles": "CC",
        "charge": 0,
        "molecular_weight": 29.06
    }
    functional_group = None
    
    # Create the combination tuple
    combination_tuple = (cation, anion, alkyl, functional_group)
    
    # Create a validator
    validator = MolecularValidator()
    
    # Validate the combination
    result = validate_combination(combination_tuple, validator)
    
    # Check the result
    if result:
        print(f"Combination validated successfully: {result['name']}")
        print(f"in_ilthermo value: {result['in_ilthermo']}")
        
        # Double-check with direct call to is_in_il_thermo
        direct_check = is_in_il_thermo(result['name'])
        print(f"Direct check with is_in_il_thermo: {direct_check}")
        
        if result['in_ilthermo'] != direct_check:
            print("WARNING: Mismatch between stored value and direct check!")
    else:
        print("Combination validation failed")
        
    # Try a direct test without validation
    il_name = generate_il_name(cation, anion, alkyl, functional_group)
    print(f"\nDirect test with generate_il_name: {il_name}")
    in_ilthermo = is_in_il_thermo(il_name)
    print(f"in_ilthermo value: {in_ilthermo}")

if __name__ == "__main__":
    test_validate_combination()
