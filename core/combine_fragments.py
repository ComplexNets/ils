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
                'h_bond_donor_count': 0,
                'h_bond_acceptor_count': 0,
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
            if not props or not all(k in props for k in ['molecular_weight', 'h_bond_donor_count', 'h_bond_acceptor_count']):
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
    screened_fragments = {'cation': [], 'anion': [], 'alkyl_chain': []}
    
    print("\nScreening fragments:")
    print(f"Density range: {density_range}")
    print(f"Heat capacity range: {cp_range}")
    
    # First calculate total properties for each possible combination
    for cation in fragments_data['cation']:
        for anion in fragments_data['anion']:
            for alkyl in fragments_data['alkyl_chain']:
                print(f"\nChecking combination:")
                print(f"  Cation: {cation['name']}")
                print(f"  Anion: {anion['name']}")
                print(f"  Alkyl: {alkyl['name']}")
                
                # Simple addition of molecular weights
                total_mw = (
                    cation['molecular_weight'] +
                    anion['molecular_weight'] +
                    alkyl['molecular_weight']
                )
                print(f"  Total MW: {total_mw}")
                
                # Estimate density based on molecular weight and heavy atoms
                total_heavy_atoms = (
                    cation['heavy_atom_count'] +
                    anion['heavy_atom_count'] +
                    alkyl['heavy_atom_count']
                )
                # Rough density estimate: 20 g/mol per heavy atom
                estimated_density = (total_mw / (total_heavy_atoms * 20.0)) * 1000  # Convert to kg/m³
                print(f"  Estimated density: {estimated_density}")
                
                # Estimate heat capacity based on molecular weight and rotatable bonds
                total_rotatable_bonds = (
                    cation['rotatable_bond_count'] +
                    anion['rotatable_bond_count'] +
                    alkyl['rotatable_bond_count']
                )
                # Rough heat capacity estimate: 2 J/mol·K per g/mol + 5 J/mol·K per rotatable bond
                estimated_cp = 2 * total_mw + 5 * total_rotatable_bonds
                print(f"  Estimated heat capacity: {estimated_cp}")
                
                # If properties are within range, keep all fragments
                if (density_range[0] <= estimated_density <= density_range[1] and
                    cp_range[0] <= estimated_cp <= cp_range[1]):
                    print("  Combination passes screening!")
                    if cation not in screened_fragments['cation']:
                        screened_fragments['cation'].append(cation)
                    if anion not in screened_fragments['anion']:
                        screened_fragments['anion'].append(anion)
                    if alkyl not in screened_fragments['alkyl_chain']:
                        screened_fragments['alkyl_chain'].append(alkyl)
                else:
                    print("  Combination fails screening")
                    if estimated_density < density_range[0]:
                        print("    Density too low")
                    elif estimated_density > density_range[1]:
                        print("    Density too high")
                    if estimated_cp < cp_range[0]:
                        print("    Heat capacity too low")
                    elif estimated_cp > cp_range[1]:
                        print("    Heat capacity too high")
    
    print("\nScreening results:")
    print(f"Screened cations: {[c['name'] for c in screened_fragments['cation']]}")
    print(f"Screened anions: {[a['name'] for a in screened_fragments['anion']]}")
    print(f"Screened alkyls: {[a['name'] for a in screened_fragments['alkyl_chain']]}")
    
    return screened_fragments

def combine_fragments() -> List[Dict]:
    """
    Combines fragments into valid ionic liquid combinations.
    Returns:
        List of valid ionic liquid combinations
    """
    try:
        # Get all fragments with properties
        fragments_data = get_filtered_fragments()
        
        valid_combinations = []
        validator = MolecularValidator()
        
        print("\nGenerating combinations:")
        # Generate combinations from fragments
        for cation in fragments_data['cation']:
            for anion in fragments_data['anion']:
                for alkyl in fragments_data['alkyl_chain']:
                    print(f"\nChecking combination:")
                    print(f"  Cation: {cation['name']}")
                    print(f"  Anion: {anion['name']}")
                    print(f"  Alkyl: {alkyl['name']}")
                    
                    # Validate the combination using chemical rules
                    is_valid, message = validator.validate(cation, anion, alkyl)
                    
                    if is_valid:
                        print("  Valid combination!")
                        # Generate the IL name
                        il_name = generate_il_name(cation, anion, alkyl)
                        
                        # Check if this IL exists in ILThermo
                        in_ilthermo = is_in_il_thermo(il_name)
                        
                        # Create the combination object
                        combination = {
                            'name': il_name,
                            'cation': cation,
                            'anion': anion,
                            'alkyl_chain': alkyl,
                            'in_ilthermo': in_ilthermo
                        }
                        
                        valid_combinations.append(combination)
                    else:
                        print(f"  Invalid combination: {message}")
        
        print(f"\nFound {len(valid_combinations)} valid combinations")
        return valid_combinations
        
    except Exception as e:
        print(f"Error in combine_fragments: {str(e)}")
        return []

if __name__ == "__main__":
    print("\n=== Testing Fragment Combination ===\n")
    combinations = combine_fragments()
    print(f"\nFound {len(combinations)} valid combinations")
