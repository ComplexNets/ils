from core.fragment_combiner import get_filtered_fragments
from itertools import product

def debug_combinations():
    # Get all fragments
    print("\n=== Getting Fragments ===")
    fragments = get_filtered_fragments('short')
    
    # Print fragment counts
    print("\n=== Fragment Counts ===")
    for frag_type, frags in fragments.items():
        print(f"{frag_type}: {len(frags)} fragments")
        print("Names:", [f['name'] for f in frags])
    
    # Create all possible combinations
    print("\n=== Creating Combinations ===")
    combinations = list(product(
        fragments['cation'],
        fragments['anion'],
        fragments['alkyl_chain'],
        fragments.get('functional_group', [None])  # Use [None] if no functional groups
    ))
    
    print(f"\nTotal possible combinations: {len(combinations)}")
    
    # Print first few combinations
    print("\n=== Sample Combinations ===")
    for i, combo in enumerate(combinations[:5]):  # Show first 5 combinations
        print(f"\nCombination {i+1}:")
        cation, anion, alkyl, func = combo
        print(f"Cation: {cation['name']}")
        print(f"Anion: {anion['name']}")
        print(f"Alkyl Chain: {alkyl['name']}")
        print(f"Functional Group: {func['name'] if func else 'None'}")

if __name__ == "__main__":
    debug_combinations()
