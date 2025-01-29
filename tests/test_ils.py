from core.combine_fragments import combine_fragments, get_filtered_fragments
from utils.validation_rules import MolecularValidator
from itertools import product

def test_fragment_loading():
    """Test that fragments are loaded correctly"""
    print("\n=== Testing Fragment Loading ===")
    fragments = get_filtered_fragments('short')
    
    # Print counts for each fragment type
    for frag_type, frags in fragments.items():
        print(f"\n{frag_type.upper()}:")
        print(f"Count: {len(frags)}")
        print("Names:", [f['name'] for f in frags])

def test_combinations():
    """Test the combination generation and validation"""
    print("\n=== Testing Combination Generation ===")
    
    # Get fragments and create combinations manually
    fragments = get_filtered_fragments('short')
    combinations = list(product(
        fragments['cation'],
        fragments['anion'],
        fragments['alkyl_chain'],
        fragments.get('functional_group', [None])
    ))
    
    print(f"\nTotal possible combinations: {len(combinations)}")
    
    # Test a few combinations with the validator
    print("\n=== Testing Validation ===")
    validator = MolecularValidator()
    
    for i, combo in enumerate(combinations[:3]):  # Test first 3 combinations
        cation, anion, alkyl, func = combo
        fragments = [cation, anion, alkyl]
        if func:
            fragments.append(func)
            
        print(f"\nTesting combination {i+1}:")
        print(f"Cation: {cation['name']}")
        print(f"Anion: {anion['name']}")
        print(f"Alkyl: {alkyl['name']}")
        print(f"Functional Group: {func['name'] if func else 'None'}")
        
        is_valid, message = validator.validate_combination(fragments)
        print(f"Valid: {is_valid}")
        print(f"Message: {message}")

def test_full_pipeline():
    """Test the full combination and validation pipeline"""
    print("\n=== Testing Full Pipeline ===")
    valid_combinations = combine_fragments(list_name='short')
    
    print(f"\nTotal valid combinations: {len(valid_combinations)}")
    
    # Print some sample valid combinations
    print("\n=== Sample Valid Combinations ===")
    for i, combo in enumerate(valid_combinations[:5]):  # Show first 5 valid combinations
        print(f"\nValid Combination {i+1}:")
        print(f"Name: {combo['name']}")
        print(f"Cation: {combo['cation']['name']}")
        print(f"Anion: {combo['anion']['name']}")
        print(f"Alkyl Chain: {combo['alkyl_chain']['name']}")
        if 'functional_group' in combo and combo['functional_group']:
            print(f"Functional Group: {combo['functional_group']['name']}")
        print(f"In ILThermo: {combo['in_ilthermo']}")

if __name__ == "__main__":
    # Run all tests
    test_fragment_loading()
    test_combinations()
    test_full_pipeline()
