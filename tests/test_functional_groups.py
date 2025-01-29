from core.combine_fragments import get_filtered_fragments, combine_fragments
from utils.validation_rules import MolecularValidator
from pprint import pprint

def test_functional_groups():
    print("\n=== Testing Functional Group Integration ===\n")
    
    # 1. Test fragment loading
    print("1. Testing fragment loading...")
    fragments = get_filtered_fragments('short')
    
    # Print functional groups
    print("\nFunctional Groups:")
    for fg in fragments['functional_group']:
        print(f"  - {fg['name']}")
    
    # 2. Test combinations with functional groups
    print("\n2. Testing combinations with functional groups...")
    valid_combinations = combine_fragments(list_name='short')
    
    # Count combinations by functional group
    fg_counts = {}
    for combo in valid_combinations:
        fg_name = combo.get('functional_group', {}).get('name', 'None')
        fg_counts[fg_name] = fg_counts.get(fg_name, 0) + 1
    
    print("\nCombinations by functional group:")
    for fg, count in fg_counts.items():
        print(f"  - {fg}: {count} combinations")
    
    # 3. Print sample combinations for each functional group
    print("\n3. Sample combinations for each functional group:")
    for fg in fragments['functional_group']:
        print(f"\nWith {fg['name']}:")
        # Find first combination with this functional group
        for combo in valid_combinations:
            if combo.get('functional_group', {}).get('name') == fg['name']:
                print(f"  Cation: {combo['cation']['name']}")
                print(f"  Anion: {combo['anion']['name']}")
                print(f"  Alkyl Chain: {combo['alkyl_chain']['name']}")
                print(f"  Functional Group: {combo['functional_group']['name']}")
                print(f"  Name: {combo['name']}")
                break
    
    # 4. Validate a specific combination
    print("\n4. Testing specific combination validation...")
    validator = MolecularValidator()
    test_combo = [
        fragments['cation'][0],  # First cation
        fragments['anion'][0],   # First anion
        fragments['alkyl_chain'][0],  # First alkyl chain
        fragments['functional_group'][0]  # First functional group
    ]
    
    is_valid, message = validator.validate_combination(test_combo)
    print(f"\nValidation result for test combination:")
    print(f"  Valid: {is_valid}")
    print(f"  Message: {message}")
    print(f"  Components: {[f['name'] for f in test_combo]}")

if __name__ == "__main__":
    test_functional_groups()
