"""
Test script to check if a specific ionic liquid is in ILThermo
"""
from utils.utils import is_in_il_thermo

def test_specific_il():
    """Test if a specific ionic liquid is in ILThermo"""
    
    # The specific ionic liquid we're interested in
    target_il = "1-ethyl-3-methylimidazolium tetrafluoroborate"
    
    print(f"Testing if '{target_il}' is in ILThermo:")
    print("-" * 50)
    
    result = is_in_il_thermo(target_il)
    print(f"Result: {'✅ Found in ILThermo' if result else '❌ Not found in ILThermo'}")
    
    # Try with different capitalization and spacing
    variant = "1-Ethyl-3-methylimidazolium    tetrafluoroborate"
    print(f"\nTesting variant: '{variant}'")
    result = is_in_il_thermo(variant)
    print(f"Result: {'✅ Found in ILThermo' if result else '❌ Not found in ILThermo'}")
    
    # Try with abbreviated name
    abbreviated = "EMIM BF4"
    print(f"\nTesting abbreviated: '{abbreviated}'")
    result = is_in_il_thermo(abbreviated)
    print(f"Result: {'✅ Found in ILThermo' if result else '❌ Not found in ILThermo'}")

if __name__ == "__main__":
    test_specific_il()
