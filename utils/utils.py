# il_naming.py - functions to standardize and generate IL names 
import re
import sys
import os
import requests
import pandas as pd
from typing import Dict, Optional
import json # Import json for JSONDecodeError

# Add project root to Python path
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path:
    sys.path.append(project_root)

from utils.rdkit_utils import get_rdkit_properties
from models.shortList_frag import fragments as short_fragments
from models.mediumList_frag import fragments as medium_fragments
from models.longList_frag import fragments as long_fragments

# Cache for fragment data DataFrame
_fragment_data_cache = None

__all__ = ['standardize_il_name', 'generate_il_name', 'get_molecular_weight', 'get_fragment_properties', 'is_in_il_thermo', 'clear_fragment_cache']

def standardize_il_name(name: str) -> str:
    """Standardize ionic liquid name to match IL Thermo format with lowercase conventions"""
    replacements = {
        # Alkyl Chains
        'butyl': 'butyl',
        'methyl': 'methyl',
        'ethyl': 'ethyl',
        'propyl': 'propyl',
        'hexyl': 'hexyl',
        'octyl': 'octyl',
        
        # Anions
        'tetrafluoroborate': 'tetrafluoroborate',
        'hexafluorophosphate': 'hexafluorophosphate',
        'chloride': 'chloride',
        'bromide': 'bromide',
        'iodide': 'iodide',
        'nitrate': 'nitrate',
        'trifluoromethanesulfonate': 'trifluoromethanesulfonate',
        'bis(trifluoromethanesulfonyl)imide': 'bis(trifluoromethanesulfonyl)imide',
        'acetate': 'acetate',
        'dicyanamide': 'dicyanamide',
        'hydrogen sulfate': 'hydrogensulfate'
    }

    # Standardize to lowercase and apply replacements
    standardized = name.lower()
    for old, new in replacements.items():
        standardized = re.sub(rf'\b{old}\b', new, standardized, flags=re.IGNORECASE)
    
    return standardized

def generate_il_name(cation: Dict, anion: Dict, alkyl: Dict, functional_group: Dict = None) -> str:
    """
    Generate standardized IUPAC-compliant name for ionic liquid combination
    
    Following proper naming conventions:
    - Imidazolium: 1-alkyl-3-methylimidazolium (or 1,3-dialkylimidazolium for identical substituents)
    - Pyridinium: N-alkylpyridinium or 1-alkylpyridinium
    - Ammonium: alkyltrimethylammonium (or tetramethylammonium for methyl)
    - Pyrrolidinium: N-alkyl-N-methylpyrrolidinium
    - Piperidinium: N-alkyl-N-methylpiperidinium
    - Functional groups should specify position: 1-(2-hydroxyethyl)-3-methylimidazolium
    """
    try:
        cation_name = cation['name'].lower()
        anion_name = anion['name'].lower()
        base_alkyl_name = alkyl['name'].lower() # Store original alkyl name

        # === START NEW Functional Group Handling ===
        # If functional group exists, modify the alkyl name *before* cation formatting
        if functional_group:
            functional_group_name = functional_group['name'].lower()
            # Determine the functionalized alkyl part string (assuming positions)
            if functional_group_name == "hydroxyl":
                # Assume position 2 for hydroxyl based on common examples
                alkyl_name = f"(2-hydroxy{base_alkyl_name})"
            elif functional_group_name == "amino":
                # Assume position 3 for amino based on common examples
                alkyl_name = f"(3-amino{base_alkyl_name})"
            else:
                # Generic fallback - just prepend functional group name
                # This might need refinement based on specific functional groups
                alkyl_name = f"({functional_group_name}{base_alkyl_name})"
        else:
             # No functional group, use the base alkyl name
             alkyl_name = base_alkyl_name
        # === END NEW Functional Group Handling ===

        # Format cation naming by checking specific cation type (using potentially modified alkyl_name)
        if "imidazolium" in cation_name:
            # Special case for identical substituents (check base_alkyl_name here)
            if base_alkyl_name == "methyl" and not functional_group: # Only if methyl is the *only* substituent
                 formatted_cation = f"1,3-dimethylimidazolium"
            # The 'else' now correctly handles non-methyl base alkyls AND functionalized alkyls
            else:
                 formatted_cation = f"1-{alkyl_name}-3-methylimidazolium"
        elif "pyridinium" in cation_name:
            # Both N- and 1- prefixes are acceptable for pyridinium
            formatted_cation = f"1-{alkyl_name}pyridinium"
        elif "ammonium" in cation_name:
            # Updated naming convention for quaternary ammonium used in IL literature
            # Special case for tetramethylammonium (only if base is methyl and no FG)
            if base_alkyl_name == "methyl" and not functional_group:
                formatted_cation = "tetramethylammonium"
            else:
                # Handles non-methyl base alkyls AND functionalized alkyls
                formatted_cation = f"{alkyl_name}trimethylammonium"
        elif "phosphonium" in cation_name:
            # Proper naming for phosphonium
            formatted_cation = f"tetra{alkyl_name}phosphonium"
            # Alternative: P,P,P-trihexyl-P-{alkyl_name}phosphonium
        elif "pyrrolidinium" in cation_name:
            # Pyrrolidinium cations
            formatted_cation = f"N-{alkyl_name}-N-methylpyrrolidinium"
        elif "piperidinium" in cation_name:
            # Piperidinium cations
            formatted_cation = f"N-{alkyl_name}-N-methylpiperidinium"
        else:
            # Generic fallback for other cation types
            formatted_cation = f"{alkyl_name}-{cation_name}"
        
        # Anion mapping with standardized names (functional group logic moved earlier)
        anion_mapping = {
            'chloride': 'chloride',
            'bromide': 'bromide',
            'iodide': 'iodide',
            'tetrafluoroborate': 'tetrafluoroborate',
            'hexafluorophosphate': 'hexafluorophosphate',
            'bis(trifluoromethanesulfonyl)imide': 'bis(trifluoromethanesulfonyl)imide',
            'trifluoromethanesulfonate': 'trifluoromethanesulfonate',
            'acetate': 'acetate',
            'nitrate': 'nitrate',
            'dicyanamide': 'dicyanamide',
            'hydrogen sulfate': 'hydrogensulfate',
            'thiocyanate': 'thiocyanate',
            'trifluoroacetate': 'trifluoroacetate',
            'tosylate': 'tosylate',
            'methanesulfonate': 'methanesulfonate',
            'trifluoromethanesulfonimide': 'trifluoromethanesulfonimide'
        }
        
        formatted_anion = anion_mapping.get(anion_name, standardize_il_name(anion_name))
        
        # Combine the formatted names with a space between cation and anion
        il_name = f"{formatted_cation} {formatted_anion}"
        return il_name
        
    except KeyError as e:
        print(f"Error in IL name generation: Missing key {e}")
        return "unknown-il"
    except Exception as e:
        print(f"Unexpected error: {e}")
        return "unknown-il"

def _get_alkyl_carbon_count(alkyl_name: str) -> int:
    """Helper function to get the carbon count from an alkyl name"""
    alkyl_carbon_counts = {
        'methyl': 1,
        'ethyl': 2,
        'propyl': 3,
        'butyl': 4,
        'pentyl': 5,
        'hexyl': 6,
        'heptyl': 7,
        'octyl': 8,
        'nonyl': 9,
        'decyl': 10,
        'undecyl': 11,
        'dodecyl': 12
    }
    return alkyl_carbon_counts.get(alkyl_name.lower(), 0)

def get_molecular_weight(fragment_name: str, fragment_type: str) -> float:
    """Get molecular weight for a fragment from the CSV data"""
    try:
        # Get properties from CSV
        props = get_fragment_properties(fragment_name, fragment_type)
        if props and 'molecular_weight' in props:
            return float(props['molecular_weight'])
        return 0.0
        
    except Exception as e:
        print(f"Error getting molecular weight: {e}")
        return 0.0

def load_fragment_data_from_csv() -> pd.DataFrame:
    """Load fragment data from local CSV file with caching"""
    global _fragment_data_cache
    
    # Return cached data if available
    if _fragment_data_cache is not None:
        return _fragment_data_cache
        
    try:
        csv_path = os.path.join(project_root, 'fragment_data', 'autono17_ilselect_db.csv')
        print(f"\nLooking for fragment data at: {csv_path}")
        
        if not os.path.exists(csv_path):
            print(f"ERROR: Fragment data CSV file not found at {csv_path}")
            print("Please ensure the fragment data file exists at this location")
            return pd.DataFrame()
            
        df = pd.read_csv(csv_path)
        print(f"Successfully loaded fragment data: {len(df)} rows")
        
        # Convert numeric columns that might have NaN to float
        numeric_columns = ['molecular_weight', 'density', 'specific_heat_capacity', 
                         'hydrogen_bond_donor_count', 'hydrogen_bond_acceptor_count',
                         'rotatable_bond_count', 'charge', 'heavy_atom_count', 'tpsa']
        for col in numeric_columns:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors='coerce')
            else:
                print(f"WARNING: Expected column {col} not found in fragment data")
        
        # Store in cache
        _fragment_data_cache = df
                
        return df
    except Exception as e:
        print(f"ERROR loading fragment data from CSV: {str(e)}")
        return pd.DataFrame()

def get_fragment_properties(fragment_name: str, fragment_type: str) -> Optional[Dict]:
    """
    Get fragment properties from cached CSV file or calculate them using RDKit
    
    Uses a cached version of the CSV data for improved performance. The cache is
    loaded only once during program execution, significantly reducing disk I/O.
    """
    try:
        print(f"\nGetting properties for {fragment_type} {fragment_name}")
        # First try to get from CSV
        df = load_fragment_data_from_csv()
        if not df.empty:
            # Convert fragment name to lowercase for comparison
            df['name'] = df['name'].str.lower()
            result = df[df['name'] == fragment_name.lower()]
            if not result.empty:
                props = {
                    'molecular_weight': float(result['molecular_weight'].iloc[0]) if not pd.isna(result['molecular_weight'].iloc[0]) else 0.0,
                    'density': float(result['density'].iloc[0]) if not pd.isna(result['density'].iloc[0]) else None,
                    'specific_heat_capacity': float(result['specific_heat_capacity'].iloc[0]) if not pd.isna(result['specific_heat_capacity'].iloc[0]) else None,
                    'hydrogen_bond_donor_count': int(result['hydrogen_bond_donor_count'].iloc[0]) if not pd.isna(result['hydrogen_bond_donor_count'].iloc[0]) else 0,
                    'hydrogen_bond_acceptor_count': int(result['hydrogen_bond_acceptor_count'].iloc[0]) if not pd.isna(result['hydrogen_bond_acceptor_count'].iloc[0]) else 0,
                    'rotatable_bond_count': int(result['rotatable_bond_count'].iloc[0]) if not pd.isna(result['rotatable_bond_count'].iloc[0]) else 0,
                    'charge': int(result['charge'].iloc[0]) if not pd.isna(result['charge'].iloc[0]) else 0,
                    'heavy_atom_count': int(result['heavy_atom_count'].iloc[0]) if not pd.isna(result['heavy_atom_count'].iloc[0]) else 0,
                    'valence': int(result['charge'].iloc[0]) if not pd.isna(result['charge'].iloc[0]) else 0
                }
                print(f"Found properties in database for {fragment_name}")
                return props
                
        # If not in CSV, try to calculate using RDKit
        print(f"Calculating properties using RDKit for {fragment_name}")
        # Find the fragment in our fragment lists
        fragment = None
        for fragments in [short_fragments, medium_fragments, long_fragments]:
            for frag in fragments:
                if frag['name'].lower() == fragment_name.lower() and frag['fragment_type'].lower() == fragment_type.lower():
                    fragment = frag
                    break
            if fragment:
                break
                
        if fragment is None:
            print(f"Error: Could not find fragment {fragment_name} of type {fragment_type} in any fragment lists")
            return None
            
        if fragment_type.lower() == 'cation':
            charge = 1
        elif fragment_type.lower() == 'anion':
            charge = -1
        else:
            charge = 0
            
        rdkit_props = get_rdkit_properties(fragment['smiles'])
        if rdkit_props:
            props = {
                'molecular_weight': rdkit_props.get('molecular_weight', 0.0),
                'density': None,  # Can't calculate
                'specific_heat_capacity': None,  # Can't calculate
                'hydrogen_bond_donor_count': rdkit_props.get('hydrogen_bond_donor_count', 0),
                'hydrogen_bond_acceptor_count': rdkit_props.get('hydrogen_bond_acceptor_count', 0),
                'rotatable_bond_count': rdkit_props.get('rotatable_bond_count', 0),
                'charge': charge,  # Use expected charge based on type
                'heavy_atom_count': rdkit_props.get('heavy_atom_count', 0),
                'valence': charge,  # Use charge as valence for now
                'smiles': fragment['smiles']  # Add SMILES to properties
            }
            print(f"Calculated properties for {fragment_name}")
            return props
            
        print(f"Error: Failed to calculate properties for {fragment_name}")
        return None
            
    except Exception as e:
        print(f"Error getting fragment properties: {str(e)}")
        return None

def clear_fragment_cache():
    """Clear the fragment data cache to force reloading from CSV file"""
    global _fragment_data_cache
    _fragment_data_cache = None
    print("Fragment data cache cleared")

def is_in_il_thermo(il_name):
    """Check if ionic liquid exists in IL Thermo database."""
    # Known ILThermo compounds (common ones)
    known_ilthermo_compounds = [
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
    # Standardize the name using the dedicated function
    search_name = standardize_il_name(il_name)

    # Check for exact match in our known list
    for known_compound in known_ilthermo_compounds:
        # Standardize the known compound name as well for consistent comparison
        if search_name == standardize_il_name(known_compound):
            print(f"DEBUG: ILThermo match found in local list: {search_name}")
            return True
    
    # Removed flawed partial match logic

    # If not in known list, try the API
    try:
        base_url = "https://ilthermo.boulder.nist.gov/ILT2/ilsearch"
        headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36',
            'Accept': 'application/json'
        }
        params = {'cmp': search_name, 'orderby': 'T', 'output': 'json'}
        
        print(f"DEBUG: Querying ILThermo API for: {search_name}")
        # Use a longer timeout and specific exception handling
        response = requests.get(base_url, params=params, headers=headers, timeout=10) # Increased timeout to 10s
        
        response.raise_for_status() # Raise HTTPError for bad responses (4xx or 5xx)

        data = response.json()
        if 'res' in data and data['res']:
            print(f"DEBUG: ILThermo API found match for: {search_name}")
            return True
        else:
            print(f"DEBUG: ILThermo API found no match for: {search_name}")
            return False

    except requests.exceptions.Timeout:
        print(f"Warning: ILThermo API query timed out for '{search_name}'. Assuming not present.")
        return False
    except requests.exceptions.RequestException as e:
        print(f"Warning: ILThermo API query failed for '{search_name}': {e}. Assuming not present.")
        return False
    except json.JSONDecodeError:
        print(f"Warning: Failed to decode ILThermo API response for '{search_name}'. Assuming not present.")
        return False
    except Exception as e:
        # Catch any other unexpected errors during the API call
        print(f"Warning: Unexpected error during ILThermo check for '{search_name}': {e}. Assuming not present.")
        return False

if __name__ == "__main__":
    # Test database connection
    pass
