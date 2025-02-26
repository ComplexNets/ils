"""
Test script for the is_in_il_thermo function
"""
from utils.utils import is_in_il_thermo

def test_is_in_ilthermo():
    """Test the is_in_il_thermo function with various ILs"""
    test_cases = [
        "1-ethyl-3-methylimidazolium tetrafluoroborate",
        "1-butyl-3-methylimidazolium tetrafluoroborate",
        "1-propyl-3-methylimidazolium tetrafluoroborate",
        "ethyltrimethylammonium tetrafluoroborate"
    ]
    
    for il_name in test_cases:
        result = is_in_il_thermo(il_name)
        print(f"IL: {il_name} - In ILThermo: {'Yes' if result else 'No'}")
        print("-" * 50)

if __name__ == "__main__":
    test_is_in_ilthermo()
