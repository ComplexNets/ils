"""
Test script to check if in_ilthermo property is being correctly set in the combination object
"""
from core.combine_fragments import combine_fragments
from utils.utils import is_in_il_thermo

# Create a mock class for Streamlit objects
class MockStreamlitObject:
    def write(self, text):
        print(text)
    
    def progress(self, value):
        pass

def test_combine_fragments():
    """Test the combine_fragments function with a simple fragment list"""
    # Use a simple fragment list with known ILThermo compounds
    fragment_list = "test_fragments.csv"
    
    # Create mock Streamlit objects
    status_text = MockStreamlitObject()
    progress_bar = MockStreamlitObject()
    
    # Set validator params - use the correct parameters for MolecularValidator
    validator_params = {
        'min_total_groups': 0,
        'max_total_groups': 6
    }
    
    # Call combine_fragments
    combinations = combine_fragments(
        status_text=status_text,
        progress_bar=progress_bar,
        validator_params=validator_params,
        list_name=fragment_list
    )
    
    # Check if combinations were generated
    if combinations:
        print(f"Generated {len(combinations)} combinations")
        
        # Check for 1-ethyl-3-methylimidazolium tetrafluoroborate
        for combo in combinations:
            name = combo['name']
            in_ilthermo = combo.get('in_ilthermo', False)
            
            # Check if this is the IL we're looking for
            if "ethyl" in name.lower() and "methylimidazolium" in name.lower() and "tetrafluoroborate" in name.lower():
                print(f"Found target IL: {name}")
                print(f"in_ilthermo value: {in_ilthermo}")
                
                # Double-check with direct call to is_in_il_thermo
                direct_check = is_in_il_thermo(name)
                print(f"Direct check with is_in_il_thermo: {direct_check}")
                
                if in_ilthermo != direct_check:
                    print("WARNING: Mismatch between stored value and direct check!")
    else:
        print("No combinations generated")

if __name__ == "__main__":
    test_combine_fragments()
