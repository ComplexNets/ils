"""
Test script to verify that the ionic liquid naming function works correctly with the ILThermo database
"""
import sys
import os

# Get absolute path to project root
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
print(f"Project root: {project_root}")

# Add project root to Python path
sys.path.append(project_root)
print(f"Added {project_root} to Python path")

from utils.utils import generate_il_name, is_in_il_thermo

def test_il_naming_ilthermo():
    """Test the IL naming function with the ILThermo database"""
    print("\n=== Testing Ionic Liquid Naming with ILThermo ===\n")
    
    # Test cases for different cation types
    test_cases = [
        # Imidazolium cations with different anions
        {
            'description': 'Imidazolium cation with tetrafluoroborate anion',
            'cation': {'name': 'Imidazolium', 'fragment_type': 'cation'},
            'anion': {'name': 'Tetrafluoroborate', 'fragment_type': 'anion'},
            'alkyl': {'name': 'Butyl', 'fragment_type': 'alkyl_chain'},
            'functional_group': None
        },
        {
            'description': 'Imidazolium cation with hexafluorophosphate anion',
            'cation': {'name': 'Imidazolium', 'fragment_type': 'cation'},
            'anion': {'name': 'Hexafluorophosphate', 'fragment_type': 'anion'},
            'alkyl': {'name': 'Ethyl', 'fragment_type': 'alkyl_chain'},
            'functional_group': None
        },
        
        # Pyridinium cations
        {
            'description': 'Pyridinium cation with chloride anion',
            'cation': {'name': 'Pyridinium', 'fragment_type': 'cation'},
            'anion': {'name': 'Chloride', 'fragment_type': 'anion'},
            'alkyl': {'name': 'Ethyl', 'fragment_type': 'alkyl_chain'},
            'functional_group': None
        },
        {
            'description': 'Pyridinium cation with bromide anion',
            'cation': {'name': 'Pyridinium', 'fragment_type': 'cation'},
            'anion': {'name': 'Bromide', 'fragment_type': 'anion'},
            'alkyl': {'name': 'Butyl', 'fragment_type': 'alkyl_chain'},
            'functional_group': None
        },
        
        # Ammonium cations
        {
            'description': 'Ammonium cation with bis(trifluoromethanesulfonyl)imide anion',
            'cation': {'name': 'Ammonium', 'fragment_type': 'cation'},
            'anion': {'name': 'Bis(trifluoromethanesulfonyl)imide', 'fragment_type': 'anion'},
            'alkyl': {'name': 'Hexyl', 'fragment_type': 'alkyl_chain'},
            'functional_group': None
        },
        {
            'description': 'Ammonium cation with dicyanamide anion',
            'cation': {'name': 'Ammonium', 'fragment_type': 'cation'},
            'anion': {'name': 'Dicyanamide', 'fragment_type': 'anion'},
            'alkyl': {'name': 'Methyl', 'fragment_type': 'alkyl_chain'},
            'functional_group': None
        },
        
        # Phosphonium cations
        {
            'description': 'Phosphonium cation with bromide anion',
            'cation': {'name': 'Phosphonium', 'fragment_type': 'cation'},
            'anion': {'name': 'Bromide', 'fragment_type': 'anion'},
            'alkyl': {'name': 'Octyl', 'fragment_type': 'alkyl_chain'},
            'functional_group': None
        },
        {
            'description': 'Phosphonium cation with acetate anion',
            'cation': {'name': 'Phosphonium', 'fragment_type': 'cation'},
            'anion': {'name': 'Acetate', 'fragment_type': 'anion'},
            'alkyl': {'name': 'Butyl', 'fragment_type': 'alkyl_chain'},
            'functional_group': None
        },
        
        # With functional groups
        {
            'description': 'Imidazolium cation with hydroxyl functional group',
            'cation': {'name': 'Imidazolium', 'fragment_type': 'cation'},
            'anion': {'name': 'Tetrafluoroborate', 'fragment_type': 'anion'},
            'alkyl': {'name': 'Butyl', 'fragment_type': 'alkyl_chain'},
            'functional_group': {'name': 'Hydroxyl', 'fragment_type': 'functional_group'}
        },
        {
            'description': 'Pyridinium cation with amino functional group',
            'cation': {'name': 'Pyridinium', 'fragment_type': 'cation'},
            'anion': {'name': 'Chloride', 'fragment_type': 'anion'},
            'alkyl': {'name': 'Ethyl', 'fragment_type': 'alkyl_chain'},
            'functional_group': {'name': 'Amino', 'fragment_type': 'functional_group'}
        }
    ]
    
    # Run tests
    found_in_ilthermo = 0
    not_found_in_ilthermo = 0
    
    for i, test in enumerate(test_cases):
        print(f"Test {i+1}: {test.get('description', '')}")
        
        # Generate name
        il_name = generate_il_name(test['cation'], test['anion'], test['alkyl'], test['functional_group'])
        
        # Truncate long names for display if necessary
        display_name = il_name
        if len(display_name) > 50:
            display_name = display_name[:47] + "..."
        
        # Check if name is in ILThermo
        in_ilthermo = is_in_il_thermo(il_name)
        
        print(f"IL name: '{display_name}'")
        print(f"In ILThermo: {in_ilthermo}")
        
        if in_ilthermo:
            found_in_ilthermo += 1
        else:
            not_found_in_ilthermo += 1
        
        # Add a newline after each test for better readability
        print()
    
    # Print summary
    print("=== Test Summary ===")
    print(f"Total tests: {len(test_cases)}")
    print(f"Found in ILThermo: {found_in_ilthermo}")
    print(f"Not found in ILThermo: {not_found_in_ilthermo}")
    print(f"Found rate: {found_in_ilthermo/len(test_cases)*100:.1f}%")
    print("=== Testing Complete ===")

if __name__ == "__main__":
    test_il_naming_ilthermo()
