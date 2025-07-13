import sys
import os

# Get absolute path to project root
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
print(f"Project root: {project_root}")
print(f"Python path before: {sys.path}")

# Add project root to Python path if not already there
if project_root not in sys.path:
    sys.path.insert(0, project_root)
    print(f"Added {project_root} to Python path")

print(f"Python path after: {sys.path}")
print(f"Looking for models in: {os.path.join(project_root, 'models')}")

from dotenv import load_dotenv
from typing import Dict, List, Tuple
from rdkit import Chem # Import RDKit
from utils.validation_rules import MolecularValidator
from utils.utils import (generate_il_name, get_molecular_weight,
                       get_fragment_properties)
from models import short_fragments, medium_fragments, long_fragments
import concurrent.futures
import multiprocessing
from functools import partial
from itertools import product

def get_filtered_fragments(list_name: str = 'short') -> Dict:
    """
    Get fragments from specified list (short, medium, or long)
    Args:
        list_name: Name of fragment list to use ('short', 'medium', or 'long')
    Returns:
        Dictionary of filtered fragments by type
    """
    try:
        print("\n=== Loading Fragment List ===")
        print(f"List name: {list_name}")
        
        # Get fragments based on list name
        if list_name == 'short':
            fragments = short_fragments
        elif list_name == 'medium':
            fragments = medium_fragments
        else:
            fragments = long_fragments
        
        print(f"Found {len(fragments)} total fragments")
        print("Fragment types:")
        for f in fragments:
            print(f"- {f['fragment_type']}: {f['name']} (SMILES: {f['smiles']})")
        
        # Initialize filtered fragments dictionary
        filtered_fragments = {
            'cation': [],
            'anion': [],
            'alkyl_chain': [],
            'functional_group': []
        }
        
        print("\n=== Getting Fragment Properties ===")
        # Get properties and filter fragments
        for fragment in fragments:
            print(f"\nChecking {fragment['fragment_type']} {fragment['name']}:")
            print(f"SMILES: {fragment['smiles']}")
            
            # Set default properties based on fragment type
            default_props = {
                'molecular_weight': 0.0,
                'density': None,
                'specific_heat_capacity': None,
                'hydrogen_bond_donor_count': 0,
                'hydrogen_bond_acceptor_count': 0,
                'rotatable_bond_count': 0,
                'charge': 1 if fragment['fragment_type'] == 'cation' else (-1 if fragment['fragment_type'] == 'anion' else 0),
                'heavy_atom_count': 0,
                'valence': 1 if fragment['fragment_type'] == 'cation' else (-1 if fragment['fragment_type'] == 'anion' else 0)
            }
            
            # Get additional properties from CSV if available
            csv_props = get_fragment_properties(fragment['name'], fragment['fragment_type'])
            if csv_props:
                print(f"Properties found in CSV:")
                for k, v in csv_props.items():
                    print(f"  {k}: {v}")
                # Update default props with CSV props
                default_props.update(csv_props)
            else:
                print(f"Warning: No properties found for {fragment['name']}, using defaults")
            
            # Add properties to fragment
            fragment.update(default_props)
            
            # Add fragment to appropriate list
            frag_type = fragment['fragment_type'].lower()
            if frag_type in filtered_fragments:
                filtered_fragments[frag_type].append(fragment)
                print(f"Added to {frag_type} list")
            else:
                print(f"WARNING: Unknown fragment type {frag_type}")
        
        print("\n=== Filtered Fragment Counts ===")
        for frag_type, frags in filtered_fragments.items():
            print(f"{frag_type}: {len(frags)} fragments")
            if frag_type != 'functional_group' and len(frags) == 0:
                print(f"ERROR: No valid {frag_type} fragments found!")
                return {}
        
        return filtered_fragments
        
    except Exception as e:
        print(f"ERROR in get_filtered_fragments: {str(e)}")
        import traceback
        traceback.print_exc()
        return {}

def screen_fragments(fragments_data: Dict, density_range: tuple = (800, 2000), cp_range: tuple = (100, 400)) -> Dict:
    """Screen fragments based on simple addition of properties"""
    # Initialize screened_fragments with all fragment types from input
    screened_fragments = {frag_type: [] for frag_type in fragments_data.keys()}
    
    # Screen each fragment type
    for frag_type in fragments_data:
        # Pass through functional groups without screening
        if frag_type == 'functional_group':
            screened_fragments[frag_type] = fragments_data[frag_type]
            continue
            
        for fragment in fragments_data[frag_type]:
            # Get properties
            props = get_fragment_properties(fragment['name'], fragment['fragment_type'])
            if not props:
                continue
            
            # Add to screened list if it meets criteria
            screened_fragments[frag_type].append(fragment)
    
    return screened_fragments

def validate_combination(combination_tuple, validator=None):
    """
    Validate a single combination of fragments
    Args:
        combination_tuple: Tuple of (cation, anion, alkyl, functional_group)
        validator: MolecularValidator instance (optional)
    Returns:
        Valid combination dict or None if invalid
    """
    try:
        # Use the validator instance passed as an argument
        if validator is None:
             # Fallback if somehow no validator was passed (shouldn't happen with current code)
             print("WARNING: No validator instance provided to validate_combination, creating default.")
             validator = MolecularValidator()
            
        cation, anion, alkyl, functional_group = combination_tuple
        
        # Create list of fragments for validation
        fragments = [cation, anion, alkyl]
        if functional_group:
            fragments.append(functional_group)

        # Prepare arguments for the validator.validate method
        alkyl_chains_list = [alkyl]
        # Functional groups need to be in a list of lists, one per chain
        functional_groups_list = [[functional_group]] if functional_group else [[]]

        # Validate using the validator's main validate method
        is_valid, message = validator.validate(cation, anion, alkyl_chains_list, functional_groups_list)
        if not is_valid:
            # print(f"Validation failed for {cation['name']}/{anion['name']}/{alkyl['name']}: {message}") # More detailed logging
            return None
            
        # Calculate total molecular weight
        total_mw = sum(get_molecular_weight(f['name'], f['fragment_type']) for f in fragments)
        
        # Generate IL name
        il_name = generate_il_name(cation, anion, alkyl, functional_group)

        # --- Generate Canonical SMILES for components ---
        canonical_smiles_dict = {}
        for frag_type, frag in [('cation', cation), ('anion', anion), ('alkyl_chain', alkyl)]:
            try:
                mol = Chem.MolFromSmiles(frag.get('smiles', ''))
                if mol:
                    canonical_smiles_dict[f"{frag_type}_canonical_smiles"] = Chem.MolToSmiles(mol, canonical=True)
                else:
                    canonical_smiles_dict[f"{frag_type}_canonical_smiles"] = None
            except Exception as e:
                print(f"Warning: RDKit error processing SMILES for {frag_type} {frag.get('name')}: {e}")
                canonical_smiles_dict[f"{frag_type}_canonical_smiles"] = None

        if functional_group:
             try:
                mol = Chem.MolFromSmiles(functional_group.get('smiles', ''))
                if mol:
                    canonical_smiles_dict["functional_group_canonical_smiles"] = Chem.MolToSmiles(mol, canonical=True)
                else:
                     canonical_smiles_dict["functional_group_canonical_smiles"] = None
             except Exception as e:
                print(f"Warning: RDKit error processing SMILES for functional_group {functional_group.get('name')}: {e}")
                canonical_smiles_dict["functional_group_canonical_smiles"] = None
        # --- End Canonical SMILES generation ---

        # Create the combination object
        result = {
            'name': il_name,
            'cation': cation,
            'anion': anion,
            'alkyl_chain': alkyl,
            'molecular_weight': total_mw,
            **canonical_smiles_dict # Add canonical smiles to the result
        }
        if functional_group:
            result['functional_group'] = functional_group
        
        print(f"Valid combination found: {il_name}")
        return result
        
    except Exception as e:
        print(f"Error in validate_combination: {str(e)}")
        return None

def combine_fragments(status_text=None, progress_bar=None, validator_params: Dict = None, list_name: str = 'short') -> List[Dict]:
    """
    Combines fragments into valid ionic liquid combinations using parallel processing.
    Args:
        status_text: Streamlit text element for status updates
        progress_bar: Streamlit progress bar element
        validator_params: Dictionary of validation parameters from UI
        list_name: Name of fragment list to use ('short', 'medium', or 'long')
    Returns:
        List of valid ionic liquid combinations
    """
    print("\n=== Starting Fragment Combination ===")
    
    # Get filtered fragments
    filtered_fragments = get_filtered_fragments(list_name)
    if not filtered_fragments:
        print("No fragments found after filtering!")
        return []

    # Initialize validator. Force default initialization to avoid TypeError due to environment issues.
    # NOTE: This will ignore validator_params from the UI.
    print("WARNING: Forcing default MolecularValidator initialization due to persistent TypeError. UI parameters ignored.")
    validator = MolecularValidator()

    # Generate all possible combinations
    cations = filtered_fragments.get('cation', [])
    anions = filtered_fragments.get('anion', [])
    alkyl_chains = filtered_fragments.get('alkyl_chain', [])
    functional_groups = filtered_fragments.get('functional_group', [])
    
    print(f"\nGenerating combinations from:")
    print(f"- {len(cations)} cations")
    print(f"- {len(anions)} anions")
    print(f"- {len(alkyl_chains)} alkyl chains")
    print(f"- {len(functional_groups)} functional groups")
    
    # Add None to functional groups to make them optional
    functional_groups.append(None)
    
    total_combinations = len(cations) * len(anions) * len(alkyl_chains) * len(functional_groups)
    print(f"\nTotal possible combinations: {total_combinations}")
    
    try:
        # Create all possible combinations including functional groups
        combinations = list(product(
            cations,
            anions,
            alkyl_chains,
            functional_groups  # Include functional groups if they exist
        ))
        
        if status_text:
            status_text.write(f"Validating combinations: 0/{total_combinations}")
        
        # Process combinations in parallel
        valid_combinations = []
        processed = 0
        
        print("\nStarting combination validation...")
        # Use ThreadPoolExecutor instead of ProcessPoolExecutor for Streamlit compatibility
        with concurrent.futures.ThreadPoolExecutor(max_workers=multiprocessing.cpu_count()) as executor:
            # Submit all combinations for processing, ensuring the correctly initialized validator is passed
            # The 'validator' variable here holds the instance created with validation_criteria=validator_params
            future_to_combo = {executor.submit(validate_combination, combo, validator): combo
                             for combo in combinations}

            # Process results as they complete
            for future in concurrent.futures.as_completed(future_to_combo):
                processed += 1
                if processed % 100 == 0:
                    print(f"Processed {processed}/{total_combinations} combinations")
                    
                if status_text and progress_bar:
                    progress = processed / total_combinations
                    progress_bar.progress(progress)
                    status_text.write(f"Validating combinations: {processed}/{total_combinations}")
                
                # Get the result
                result = future.result()
                if result is not None:
                    valid_combinations.append(result)
                    print(f"Found valid combination: {result['name']}")
        
        print(f"\nValidation complete. Found {len(valid_combinations)} valid combinations")
        
        if status_text:
            status_text.write(f"✅ Found {len(valid_combinations)} valid combinations")
        
        return valid_combinations
    
    except Exception as e:
        print(f"Error in combine_fragments: {str(e)}")
        return []

if __name__ == "__main__":
    print("\n=== Testing Fragment Combination ===\n")
    
    # Let user choose fragment list
    available_lists = ['short', 'medium', 'long']
    print(f"Available fragment lists: {', '.join(available_lists)}")
    list_choice = input("Enter fragment list to use (short/medium/long): ").lower()
    
    try:
        combinations = combine_fragments(list_name=list_choice)
        print(f"\nFound {len(combinations)} valid combinations")
    except ValueError as e:
        print(f"Error: {str(e)}")
