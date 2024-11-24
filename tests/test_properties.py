import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core.density import test_density_calculations, screen_fragments_by_density
from core.heat_capacity import test_heat_capacity_calculations, screen_fragments_by_heat_capacity
from core.toxicity import test_toxicity_calculations, screen_fragments_by_toxicity
from core.combine_fragments import get_filtered_fragments, combine_fragments

def run_all_tests():
    """Run tests for all property calculation modules"""
    print("\n=== Running All Property Tests ===\n")
    
    # Get test data
    print("Getting test data...")
    fragments_data = get_filtered_fragments()
    combinations = combine_fragments()
    
    if not fragments_data or not combinations:
        print("Error: Could not get test data")
        return
        
    print(f"Got {sum(len(frags) for frags in fragments_data.values())} fragments")
    print(f"Got {len(combinations)} valid combinations")
    
    # Test density calculations
    print("\n" + "="*50)
    print("Testing Density Calculations")
    print("="*50)
    test_density_calculations(fragments_data, combinations)
    
    # Test heat capacity calculations
    print("\n" + "="*50)
    print("Testing Heat Capacity Calculations")
    print("="*50)
    test_heat_capacity_calculations(fragments_data, combinations)
    
    # Test toxicity calculations
    print("\n" + "="*50)
    print("Testing Toxicity Calculations")
    print("="*50)
    test_toxicity_calculations(fragments_data, combinations)
    
    # Test property ranges for all combinations
    print("\n" + "="*50)
    print("Property Ranges for All Combinations")
    print("="*50)
    
    density_range = (800, 2000)  # kg/m³
    cp_range = (100, 400)       # J/mol·K
    tox_range = (5.0, 100.0)    # mM (IC50)
    
    density_results = screen_fragments_by_density(fragments_data, density_range)
    cp_results = screen_fragments_by_heat_capacity(fragments_data, cp_range)
    tox_results = screen_fragments_by_toxicity(fragments_data, tox_range)
    
    print("\nFragments meeting all criteria:")
    for frag_type in fragments_data.keys():
        print(f"\n{frag_type.capitalize()}:")
        density_names = {f['name'] for f in density_results[frag_type]}
        cp_names = {f['name'] for f in cp_results[frag_type]}
        tox_names = {f['name'] for f in tox_results[frag_type]}
        
        # Find fragments that meet all criteria
        common_fragments = density_names & cp_names & tox_names
        
        if common_fragments:
            for name in common_fragments:
                print(f"  - {name}")
                # Get detailed properties
                density_frag = next((f for f in density_results[frag_type] if f['name'] == name), None)
                cp_frag = next((f for f in cp_results[frag_type] if f['name'] == name), None)
                tox_frag = next((f for f in tox_results[frag_type] if f['name'] == name), None)
                
                if all([density_frag, cp_frag, tox_frag]):
                    print(f"    Density: {density_frag.get('estimated_density', 'N/A'):.2f} kg/m³")
                    print(f"    Heat Capacity: {cp_frag.get('estimated_heat_capacity', 'N/A'):.2f} J/mol·K")
                    print(f"    IC50: {tox_frag.get('estimated_ic50', 'N/A'):.2f} mM")
        else:
            print("  No fragments meet all criteria")

if __name__ == "__main__":
    run_all_tests()
