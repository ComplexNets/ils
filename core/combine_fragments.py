import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pymysql
from dotenv import load_dotenv
from typing import Dict, List, Tuple
from utils.validation_rules import MolecularValidator
from utils.utils import (generate_il_name, is_in_il_thermo, get_molecular_weight, 
                       get_fragment_properties)
from models.shortList_frag import fragments
import concurrent.futures
import multiprocessing
from functools import partial
from itertools import product

def get_filtered_fragments() -> Dict:
    """Get fragments from shortList_frag.py"""
    fragments_data = {'cation': [], 'anion': [], 'alkyl_chain': []}
    
    print("\nGetting fragment properties:")
    for fragment in fragments:
        frag_type = fragment['fragment_type']
        if frag_type in fragments_data:
            print(f"\nProcessing {fragment['name']} ({frag_type}):")
            
            # Create processed fragment with default properties
            processed_fragment = fragment.copy()
            processed_fragment.update({
                'name': fragment['name'],
                'molecular_weight': 0.0,
                'log_p': 0.0,
                'hydrogen_bond_donor_count': 0,
                'hydrogen_bond_acceptor_count': 0,
                'rotatable_bond_count': 0,
                'complexity': 0.0,
                'charge': 0,
                'heavy_atom_count': 0,
                'tpsa': 0.0
            })
            
            # Try to get properties from database first
            try:
                props = get_fragment_properties(fragment['name'], frag_type)
                if props:
                    print(f"  Found database properties for {fragment['name']}")
                    processed_fragment.update(props)
            except Exception as e:
                print(f"  Database error: {str(e)}")
            
            # If database properties not found or incomplete, use RDKit
            if not props or not all(k in props for k in ['molecular_weight', 'hydrogen_bond_donor_count', 'hydrogen_bond_acceptor_count']):
                print(f"  No complete database properties for {fragment['name']}, using RDKit...")
                if 'smiles' in fragment:
                    print(f"  Found SMILES for {fragment['name']}: {fragment['smiles']}")
                    # You can add RDKit calculations here if needed
                else:
                    print(f"  Warning: No SMILES found for {fragment['name']}")
            
            fragments_data[frag_type].append(processed_fragment)
            print(f"  Added to {frag_type} list")
    
    print("\nFragment counts:")
    for frag_type, frags in fragments_data.items():
        print(f"  {frag_type}: {len(frags)}")
    
    return fragments_data

def screen_fragments(fragments_data: Dict, density_range: tuple = (800, 2000), cp_range: tuple = (100, 400)) -> Dict:
    """Screen fragments based on simple addition of properties"""
    screened_fragments = {
        'cation': [],
        'anion': [],
        'alkyl_chain': []
    }
    
    # Screen each fragment type
    for frag_type in fragments_data:
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
        combination_tuple: Tuple of (cation, anion, alkyl)
        validator: MolecularValidator instance (optional)
    Returns:
        Valid combination dict or None if invalid
    """
    if validator is None:
        validator = MolecularValidator()
        
    cation, anion, alkyl = combination_tuple
    is_valid, message = validator.validate(cation, anion, alkyl)
    
    if is_valid:
        # Generate the IL name
        il_name = generate_il_name(cation, anion, alkyl)
        
        # Check if this IL exists in ILThermo
        in_ilthermo = is_in_il_thermo(il_name)
        
        # Create the combination object
        return {
            'name': il_name,
            'cation': cation,
            'anion': anion,
            'alkyl_chain': alkyl,
            'in_ilthermo': in_ilthermo
        }
    return None

def combine_fragments(status_text=None, progress_bar=None) -> List[Dict]:
    """
    Combines fragments into valid ionic liquid combinations using parallel processing.
    Args:
        status_text: Streamlit text element for status updates
        progress_bar: Streamlit progress bar element
    Returns:
        List of valid ionic liquid combinations
    """
    try:
        # Get all fragments with properties
        fragments_data = get_filtered_fragments()
        
        # Calculate total possible combinations
        total_possible = (len(fragments_data['cation']) * 
                        len(fragments_data['anion']) * 
                        len(fragments_data['alkyl_chain']))
        
        if status_text:
            status_text.write(f"Starting validation of {total_possible} possible combinations...")
        
        # Create all possible combinations
        combinations = list(product(
            fragments_data['cation'],
            fragments_data['anion'],
            fragments_data['alkyl_chain']
        ))
        
        # Process combinations in parallel with progress tracking
        valid_combinations = []
        processed_count = 0
        
        # Use ThreadPoolExecutor instead of ProcessPoolExecutor for Streamlit compatibility
        with concurrent.futures.ThreadPoolExecutor(max_workers=multiprocessing.cpu_count()) as executor:
            # Submit all combinations for processing
            future_to_combo = {executor.submit(validate_combination, combo): combo 
                             for combo in combinations}
            
            # Process results as they complete
            for future in concurrent.futures.as_completed(future_to_combo):
                processed_count += 1
                
                # Update progress
                if status_text and progress_bar:
                    progress = processed_count / total_possible
                    progress_bar.progress(progress)
                    status_text.write(f"Validating combinations: {processed_count}/{total_possible}")
                
                # Get the result
                result = future.result()
                if result is not None:
                    valid_combinations.append(result)
                    print(f"Valid combination: {result['name']}")
        
        # Show final statistics if UI elements available
        if status_text:
            valid_percent = (len(valid_combinations) / total_possible) * 100
            status_text.write(f"✅ Validation complete: Found {len(valid_combinations)} valid combinations ({valid_percent:.1f}% of total possible)")
        
        return valid_combinations
        
    except Exception as e:
        if status_text:
            status_text.error(f"Error during combination process: {str(e)}")
        print(f"Error during combination process: {str(e)}")
        return []

if __name__ == "__main__":
    print("\n=== Testing Fragment Combination ===\n")
    combinations = combine_fragments()
    print(f"\nFound {len(combinations)} valid combinations")
