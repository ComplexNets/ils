import sys
import os

# Add project root to Python path
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from utils.utils import generate_il_name, standardize_il_name

def test_naming():
    """Test the IL naming function with various fragment combinations"""
    print("\n=== Testing Ionic Liquid Naming ===\n")
    
    # Test cases for different cation types
    test_cases = [
        # Imidazolium cations with different anions
        {
            'description': 'Imidazolium cation with tetrafluoroborate anion',
            'cation': {'name': 'Imidazolium', 'fragment_type': 'cation'},
            'anion': {'name': 'Tetrafluoroborate', 'fragment_type': 'anion'},
            'alkyl': {'name': 'Butyl', 'fragment_type': 'alkyl_chain'},
            'expected': '1-butyl-3-methylimidazolium tetrafluoroborate'
        },
        {
            'description': 'Imidazolium cation with hexafluorophosphate anion',
            'cation': {'name': 'Imidazolium', 'fragment_type': 'cation'},
            'anion': {'name': 'Hexafluorophosphate', 'fragment_type': 'anion'},
            'alkyl': {'name': 'Ethyl', 'fragment_type': 'alkyl_chain'},
            'expected': '1-ethyl-3-methylimidazolium hexafluorophosphate'
        },
        
        # Pyridinium cations
        {
            'description': 'Pyridinium cation with chloride anion',
            'cation': {'name': 'Pyridinium', 'fragment_type': 'cation'},
            'anion': {'name': 'Chloride', 'fragment_type': 'anion'},
            'alkyl': {'name': 'Ethyl', 'fragment_type': 'alkyl_chain'},
            'expected': '1-ethylpyridinium chloride'
        },
        {
            'description': 'Pyridinium cation with bromide anion',
            'cation': {'name': 'Pyridinium', 'fragment_type': 'cation'},
            'anion': {'name': 'Bromide', 'fragment_type': 'anion'},
            'alkyl': {'name': 'Butyl', 'fragment_type': 'alkyl_chain'},
            'expected': '1-butylpyridinium bromide'
        },
        
        # Ammonium cations
        {
            'description': 'Ammonium cation with bis(trifluoromethanesulfonyl)imide anion',
            'cation': {'name': 'Ammonium', 'fragment_type': 'cation'},
            'anion': {'name': 'Bis(trifluoromethanesulfonyl)imide', 'fragment_type': 'anion'},
            'alkyl': {'name': 'Hexyl', 'fragment_type': 'alkyl_chain'},
            'expected': 'N,N,N-trimethyl-N-hexylammonium bis(trifluoromethanesulfonyl)imide'
        },
        {
            'description': 'Ammonium cation with dicyanamide anion',
            'cation': {'name': 'Ammonium', 'fragment_type': 'cation'},
            'anion': {'name': 'Dicyanamide', 'fragment_type': 'anion'},
            'alkyl': {'name': 'Methyl', 'fragment_type': 'alkyl_chain'},
            'expected': 'N,N,N-trimethyl-N-methylammonium dicyanamide'
        },
        
        # Phosphonium cations
        {
            'description': 'Phosphonium cation with bromide anion',
            'cation': {'name': 'Phosphonium', 'fragment_type': 'cation'},
            'anion': {'name': 'Bromide', 'fragment_type': 'anion'},
            'alkyl': {'name': 'Octyl', 'fragment_type': 'alkyl_chain'},
            'expected': 'tetraoctylphosphonium bromide'
        },
        {
            'description': 'Phosphonium cation with acetate anion',
            'cation': {'name': 'Phosphonium', 'fragment_type': 'cation'},
            'anion': {'name': 'Acetate', 'fragment_type': 'anion'},
            'alkyl': {'name': 'Butyl', 'fragment_type': 'alkyl_chain'},
            'expected': 'tetrabutylphosphonium acetate'
        },
        
        # Pyrrolidinium cations
        {
            'description': 'Pyrrolidinium cation with dicyanamide anion',
            'cation': {'name': 'Pyrrolidinium', 'fragment_type': 'cation'},
            'anion': {'name': 'Dicyanamide', 'fragment_type': 'anion'},
            'alkyl': {'name': 'Butyl', 'fragment_type': 'alkyl_chain'},
            'expected': 'N-butyl-N-methylpyrrolidinium dicyanamide'
        },
        {
            'description': 'Pyrrolidinium cation with bis(trifluoromethanesulfonyl)imide anion',
            'cation': {'name': 'Pyrrolidinium', 'fragment_type': 'cation'},
            'anion': {'name': 'Bis(trifluoromethanesulfonyl)imide', 'fragment_type': 'anion'},
            'alkyl': {'name': 'Ethyl', 'fragment_type': 'alkyl_chain'},
            'expected': 'N-ethyl-N-methylpyrrolidinium bis(trifluoromethanesulfonyl)imide'
        },
        
        # Piperidinium cations
        {
            'description': 'Piperidinium cation with acetate anion',
            'cation': {'name': 'Piperidinium', 'fragment_type': 'cation'},
            'anion': {'name': 'Acetate', 'fragment_type': 'anion'},
            'alkyl': {'name': 'Propyl', 'fragment_type': 'alkyl_chain'},
            'expected': 'N-propyl-N-methylpiperidinium acetate'
        },
        {
            'description': 'Piperidinium cation with tetrafluoroborate anion',
            'cation': {'name': 'Piperidinium', 'fragment_type': 'cation'},
            'anion': {'name': 'Tetrafluoroborate', 'fragment_type': 'anion'},
            'alkyl': {'name': 'Butyl', 'fragment_type': 'alkyl_chain'},
            'expected': 'N-butyl-N-methylpiperidinium tetrafluoroborate'
        },
        
        # With functional groups
        {
            'description': 'Imidazolium cation with hydroxyl functional group',
            'cation': {'name': 'Imidazolium', 'fragment_type': 'cation'},
            'anion': {'name': 'Tetrafluoroborate', 'fragment_type': 'anion'},
            'alkyl': {'name': 'Butyl', 'fragment_type': 'alkyl_chain'},
            'functional_group': {'name': 'Hydroxyl', 'fragment_type': 'functional_group'},
            'expected': '1-(2-hydroxybutyl)-3-methylimidazolium tetrafluoroborate'
        },
        {
            'description': 'Pyridinium cation with amino functional group',
            'cation': {'name': 'Pyridinium', 'fragment_type': 'cation'},
            'anion': {'name': 'Chloride', 'fragment_type': 'anion'},
            'alkyl': {'name': 'Ethyl', 'fragment_type': 'alkyl_chain'},
            'functional_group': {'name': 'Amino', 'fragment_type': 'functional_group'},
            'expected': '1-(3-aminoethyl)pyridinium chloride'
        }
    ]
    
    # Run tests
    passed = 0
    failed = 0
    
    for i, test in enumerate(test_cases):
        print(f"Test {i+1}: {test.get('description', '')}")
        
        # Get functional group if present
        functional_group = test.get('functional_group', None)
        
        # Generate name
        il_name = generate_il_name(test['cation'], test['anion'], test['alkyl'], functional_group)
        
        # Truncate long names for display if necessary
        display_name = il_name
        if len(display_name) > 50:
            display_name = display_name[:47] + "..."
        
        # Check if name matches expected
        if il_name == test['expected']:
            print(f"✅ PASS: Generated name: '{display_name}'")
            passed += 1
        else:
            print(f"❌ FAIL: Generated name: '{display_name}'")
            print(f"         Expected name: '{test['expected']}'")
            failed += 1
        
        # Add a newline after each test for better readability
        print()
    
    # Print summary
    print("=== Test Summary ===")
    print(f"Total tests: {len(test_cases)}")
    print(f"Passed: {passed}")
    print(f"Failed: {failed}")
    print(f"Success rate: {passed/len(test_cases)*100:.1f}%")
    print("=== Testing Complete ===")

if __name__ == "__main__":
    test_naming()
