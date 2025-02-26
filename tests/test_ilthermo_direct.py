"""
Direct test of the is_in_il_thermo function with specific ionic liquids
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

from utils.utils import is_in_il_thermo

def test_ilthermo_direct():
    """Test the is_in_il_thermo function directly with specific ionic liquids"""
    
    # List of ionic liquids to test
    test_ils = [
        "1-ethyl-3-methylimidazolium tetrafluoroborate",
        "1-butyl-3-methylimidazolium tetrafluoroborate",
        "1-butyl-3-methylimidazolium hexafluorophosphate",
        "1-ethyl-3-methylimidazolium bis(trifluoromethanesulfonyl)imide",
        "1-butyl-3-methylimidazolium bis(trifluoromethanesulfonyl)imide",
        "1-hexyl-3-methylimidazolium bis(trifluoromethanesulfonyl)imide",
        "1-octyl-3-methylimidazolium bis(trifluoromethanesulfonyl)imide",
        "1-butyl-3-methylimidazolium chloride",
        "1-ethyl-3-methylimidazolium acetate",
        "1-butyl-3-methylimidazolium acetate",
        "methyltrimethylammonium bis(trifluoromethanesulfonyl)imide",
        "ethyltrimethylammonium bis(trifluoromethanesulfonyl)imide",
        "butyltrimethylammonium bis(trifluoromethanesulfonyl)imide",
    ]
    
    print("Testing is_in_il_thermo function directly:")
    print("-" * 50)
    
    for il in test_ils:
        result = is_in_il_thermo(il)
        print(f"IL: {il}")
        print(f"Result: {'✅ Found in ILThermo' if result else '❌ Not found in ILThermo'}")
        print("-" * 50)
    
    # Test a non-existent IL
    non_existent = "non-existent-cation non-existent-anion"
    result = is_in_il_thermo(non_existent)
    print(f"IL: {non_existent}")
    print(f"Result: {'✅ Found in ILThermo' if result else '❌ Not found in ILThermo'}")
    print("-" * 50)

if __name__ == "__main__":
    test_ilthermo_direct()
