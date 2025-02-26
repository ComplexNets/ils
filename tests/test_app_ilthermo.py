"""
Test script to check if the Streamlit app correctly identifies ILs in ILThermo
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

from core.combine_fragments import combine_fragments
from utils.utils import is_in_il_thermo

# Create a mock class for Streamlit objects
class MockStreamlitObject:
    def write(self, text):
        print(text)
    
    def progress(self, value):
        pass

def test_app_ilthermo():
    """Test if the Streamlit app correctly identifies ILs in ILThermo"""
    
    # Create mock Streamlit objects
    status_text = MockStreamlitObject()
    progress_bar = MockStreamlitObject()
    
    # Set validator params
    validator_params = {
        'min_total_groups': 0,
        'max_total_groups': 6
    }
    
    # Call combine_fragments with short list (faster)
    combinations = combine_fragments(
        status_text=status_text,
        progress_bar=progress_bar,
        validator_params=validator_params,
        list_name='short'
    )
    
    # Check if combinations were generated
    if combinations:
        print(f"\nGenerated {len(combinations)} combinations")
        
        # Count how many are in ILThermo
        in_ilthermo_count = sum(1 for combo in combinations if combo.get('in_ilthermo', False))
        print(f"Found {in_ilthermo_count} combinations in ILThermo")
        
        # Print the first 5 combinations that are in ILThermo
        print("\nFirst 5 combinations in ILThermo:")
        count = 0
        for combo in combinations:
            if combo.get('in_ilthermo', False):
                print(f"- {combo['name']}")
                count += 1
                if count >= 5:
                    break
        
        # Check for 1-ethyl-3-methylimidazolium tetrafluoroborate specifically
        target_found = False
        for combo in combinations:
            name = combo['name']
            in_ilthermo = combo.get('in_ilthermo', False)
            
            if "ethyl" in name.lower() and "methylimidazolium" in name.lower() and "tetrafluoroborate" in name.lower():
                print(f"\nFound target IL: {name}")
                print(f"in_ilthermo value: {in_ilthermo}")
                
                # Double-check with direct call to is_in_il_thermo
                direct_check = is_in_il_thermo(name)
                print(f"Direct check with is_in_il_thermo: {direct_check}")
                
                if in_ilthermo != direct_check:
                    print("WARNING: Mismatch between stored value and direct check!")
                
                target_found = True
                break
        
        if not target_found:
            print("\nTarget IL not found in generated combinations")
    else:
        print("No combinations generated")

if __name__ == "__main__":
    test_app_ilthermo()
