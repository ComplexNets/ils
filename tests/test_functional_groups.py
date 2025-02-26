import sys
import os
import pandas as pd
import re

# Get absolute path to project root
project_root = os.path.dirname(os.path.abspath(__file__))
print(f"Project root: {project_root}")

# Add project root to Python path if not already there
if project_root not in sys.path:
    sys.path.insert(0, project_root)
    print(f"Added {project_root} to Python path")

from core.combine_fragments import combine_fragments
from models.shortList_frag import fragments as short_fragments
from utils.validation_rules import MolecularValidator
from analysis.il_decon import IonicLiquidAnalyzer, process_excel_file

def analyze_functional_groups():
    """Analyze why functional groups aren't being included in the output"""
    print("\n=== Analyzing Functional Groups ===")
    
    # Print functional groups from fragment list
    print("\nFunctional groups in short_fragments:")
    func_groups = [f for f in short_fragments if f['fragment_type'] == 'functional_group']
    for fg in func_groups:
        print(f"  {fg['name']} (SMILES: {fg['smiles']})")
    
    # Get all combinations
    print("\nGetting combinations...")
    combinations = combine_fragments(list_name='short')
    print(f"Found {len(combinations)} combinations")
    
    # Check how many have functional groups
    with_func_groups = [c for c in combinations if 'functional_group' in c and c['functional_group'] is not None]
    print(f"Combinations with functional groups: {len(with_func_groups)}/{len(combinations)}")
    
    # If we have combinations with functional groups, print some examples
    if with_func_groups:
        print("\nExamples with functional groups:")
        for i, combo in enumerate(with_func_groups[:3]):
            print(f"\nCombination {i+1}:")
            print(f"  Name: {combo['name']}")
            print(f"  Functional group: {combo['functional_group']['name']}")
            
        # Create a test input file with these combinations
        create_test_input_file(with_func_groups[:10])
    else:
        print("\nNo combinations with functional groups found!")

def create_test_input_file(combinations):
    """Create a test input file with combinations that have functional groups"""
    print("\nCreating test input file...")
    
    # Create dataframe with the combinations
    data = []
    for combo in combinations:
        data.append({
            'Name': combo['name'],
            'Heat Capacity (J/mol·K)': 400.0,
            'Density (kg/m³)': 1000.0,
            'Toxicity (IC50 in mM)': 0.5,
            'Solubility (g/L)': 100.0,
            'Hydrophobicity (logP)': 0.0,
            'Viscosity (Pa·s)': 0.05,  # Valid value within the range of 0.01 to 10.0 Pa·s
            'Pareto Score': 0.0,
            'In ILThermo': 'No'
        })
    
    # Create dataframe
    df = pd.DataFrame(data)
    
    # Save to Excel
    test_input_file = 'analysis/il_prop_test_input.xlsx'
    df.to_excel(test_input_file, index=False)
    print(f"Created test input file: {test_input_file}")
    
    # Process the file with the analyzer
    test_output_file = 'analysis/il_prop_test_output.xlsx'
    process_excel_file(test_input_file, test_output_file, verbose=True)
    print(f"Processed test file and saved to: {test_output_file}")
    
    # Read the output file and check for functional groups
    df_output = pd.read_excel(test_output_file)
    print("\nFunctional groups in output file:")
    if 'functional_group' in df_output.columns:
        fg_counts = df_output['functional_group'].value_counts()
        if len(fg_counts) > 0:
            print(fg_counts.to_string())
        else:
            print("No functional groups found in output file!")
            
            # Test the analyzer directly on some example names
            print("\nTesting analyzer directly:")
            analyzer = IonicLiquidAnalyzer()
            for name in df['Name']:
                analysis = analyzer.analyze_ionic_liquid(name)
                print(f"Name: {name}")
                print(f"  Functional group: {analysis['functional_group']}")
                
                # Check if the name contains any of the functional group patterns
                for pattern, result in analyzer.functional_group_patterns:
                    if any(re.search(pattern, part, re.IGNORECASE) for part in name.split()):
                        print(f"  Pattern '{pattern}' should match '{result}'")
    else:
        print("No functional_group column in output file!")

if __name__ == "__main__":
    import re
    analyze_functional_groups()
