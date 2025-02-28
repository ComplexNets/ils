# il_naming.py - functions to standardize and generate IL names 
import re
import sys
import os
import requests
import pandas as pd
from typing import Dict, Optional

# Add project root to Python path
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path:
    sys.path.append(project_root)

from utils.rdkit_utils import get_rdkit_properties
from models.shortList_frag import fragments as short_fragments
from models.mediumList_frag import fragments as medium_fragments
from models.longList_frag import fragments as long_fragments

__all__ = ['standardize_il_name', 'generate_il_name', 'get_molecular_weight', 'get_fragment_properties', 'is_in_il_thermo']

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
        alkyl_name = alkyl['name'].lower()
        
        # Format cation naming by checking specific cation type
        if "imidazolium" in cation_name:
            # Special case for identical substituents
            if alkyl_name == "methyl":
                formatted_cation = f"1,3-dimethylimidazolium"
            else:
                formatted_cation = f"1-{alkyl_name}-3-methylimidazolium"
        elif "pyridinium" in cation_name:
            # Both N- and 1- prefixes are acceptable for pyridinium
            formatted_cation = f"1-{alkyl_name}pyridinium"
        elif "ammonium" in cation_name:
            # Updated naming convention for quaternary ammonium used in IL literature
            if alkyl_name == "methyl":
                formatted_cation = "tetramethylammonium"
            else:
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
        
        # Add functional group if present
        if functional_group:
            functional_group_name = functional_group['name'].lower()
            
            # Proper functional group naming based on type
            if functional_group_name == "hydroxyl":
                # Format as 1-(2-hydroxyethyl)-3-methylimidazolium
                if "imidazolium" in cation_name:
                    # Extract the alkyl part and add hydroxy prefix
                    alkyl_carbon_count = _get_alkyl_carbon_count(alkyl_name)
                    if alkyl_carbon_count > 0:
                        formatted_cation = f"1-(2-hydroxy{alkyl_name})-3-methylimidazolium"
                    else:
                        formatted_cation = f"1-(hydroxy{alkyl_name})-3-methylimidazolium"
                elif "pyridinium" in cation_name:
                    formatted_cation = f"1-(2-hydroxy{alkyl_name})pyridinium"
                elif "ammonium" in cation_name:
                    formatted_cation = f"(2-hydroxy{alkyl_name})trimethylammonium"
                elif "pyrrolidinium" in cation_name or "piperidinium" in cation_name:
                    formatted_cation = f"N-(2-hydroxy{alkyl_name})-N-methylpyrrolidinium"
                else:
                    # Generic fallback
                    formatted_cation = f"(2-hydroxy{alkyl_name})-{cation_name}"
            elif functional_group_name == "amino":
                # Format as 1-(3-aminopropyl)-3-methylimidazolium
                if "imidazolium" in cation_name:
                    # Extract the alkyl part and add amino prefix
                    alkyl_carbon_count = _get_alkyl_carbon_count(alkyl_name)
                    if alkyl_carbon_count > 0:
                        formatted_cation = f"1-(3-amino{alkyl_name})-3-methylimidazolium"
                    else:
                        formatted_cation = f"1-(amino{alkyl_name})-3-methylimidazolium"
                elif "pyridinium" in cation_name:
                    formatted_cation = f"1-(3-amino{alkyl_name})pyridinium"
                elif "ammonium" in cation_name:
                    formatted_cation = f"(3-amino{alkyl_name})trimethylammonium"
                elif "pyrrolidinium" in cation_name or "piperidinium" in cation_name:
                    formatted_cation = f"N-(3-amino{alkyl_name})-N-methylpyrrolidinium"
                else:
                    # Generic fallback
                    formatted_cation = f"(3-amino{alkyl_name})-{cation_name}"
            else:
                # Generic functional group handling
                formatted_cation = f"{functional_group_name}-{formatted_cation}"

        # Anion mapping with standardized names
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
    """Load fragment data from local CSV file"""
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
                
        return df
    except Exception as e:
        print(f"ERROR loading fragment data from CSV: {str(e)}")
        return pd.DataFrame()

def get_fragment_properties(fragment_name: str, fragment_type: str) -> Optional[Dict]:
    """
    Get fragment properties from CSV file or calculate them using RDKit
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
    
    # Standardize the name: ensure proper spacing and lowercase
    search_name = ' '.join(il_name.lower().split())
    
    # Special case for known compounds
    if ("ethyl" in search_name and "methylimidazolium" in search_name and "tetrafluoroborate" in search_name) or \
       ("butyl" in search_name and "methylimidazolium" in search_name and "tetrafluoroborate" in search_name):
        return True
    
    # Check for exact match in our known list
    for known_compound in known_ilthermo_compounds:
        if search_name == known_compound.lower():
            return True
    
    # Check for partial matches in our known list
    for known_compound in known_ilthermo_compounds:
        # Split into parts and check if all parts are in the search name
        known_parts = known_compound.lower().split()
        search_parts = search_name.split()
        
        # Check if all known parts are in search parts
        if all(part in search_parts for part in known_parts):
            return True
    
    # If not in known list, try the API
    try:
        base_url = "https://ilthermo.boulder.nist.gov/ILT2/ilsearch"
        headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36',
            'Accept': 'application/json'
        }
        params = {'cmp': search_name, 'orderby': 'T', 'output': 'json'}
        
        # Use a short timeout to avoid hanging
        response = requests.get(base_url, params=params, headers=headers, timeout=3)
        
        if response.status_code == 200:
            try:
                data = response.json()
                if 'res' in data and data['res']:
                    return True
                return False
            except:
                return False
        else:
            return False
            
    except Exception:
        # If API fails, just rely on our known list (which we already checked)
        return False

if __name__ == "__main__":
    # Test database connection
    pass