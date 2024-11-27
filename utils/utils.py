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
from models.shortList_frag import fragments

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

def generate_il_name(cation: Dict, anion: Dict, alkyl: Dict) -> str:
    """Generate standardized lowercase name for ionic liquid combination"""
    try:
        cation_name = cation['name'].lower()
        anion_name = anion['name'].lower()
        alkyl_name = alkyl['name'].lower()
        
        # Format cation naming by checking specific cation type
        if "imidazolium" in cation_name:
            formatted_cation = f"1-{alkyl_name}-3-methylimidazolium"
        elif "pyridinium" in cation_name:
            formatted_cation = f"1-{alkyl_name}pyridinium"
        elif "ammonium" in cation_name:
            formatted_cation = f"{alkyl_name}trimethylammonium"
        elif "phosphonium" in cation_name:
            formatted_cation = f"{alkyl_name}phosphonium"
        else:
            formatted_cation = f"{alkyl_name}-{cation_name}"

        # Retrieve mapped name if available or fallback
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
            'hydrogen sulfate': 'hydrogensulfate'
        }
        
        formatted_anion = anion_mapping.get(anion_name, standardize_il_name(anion_name))
        
        # Combine the formatted names
        il_name = f"{formatted_cation} {formatted_anion}"
        return il_name
        
    except KeyError as e:
        print(f"Error in IL name generation: Missing key {e}")
        return "unknown-il"
    except Exception as e:
        print(f"Unexpected error: {e}")
        return "unknown-il"

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
        if not os.path.exists(csv_path):
            print(f"Error: CSV file not found at {csv_path}")
            return pd.DataFrame()
            
        df = pd.read_csv(csv_path)
        # Convert numeric columns that might have NaN to float
        numeric_columns = ['molecular_weight', 'density', 'specific_heat_capacity', 
                         'hydrogen_bond_donor_count', 'hydrogen_bond_acceptor_count',
                         'rotatable_bond_count', 'charge', 'heavy_atom_count', 'tpsa']
        for col in numeric_columns:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors='coerce')
        return df
    except Exception as e:
        print(f"Error loading fragment data from CSV: {str(e)}")
        return pd.DataFrame()

def get_fragment_properties(fragment_name: str, fragment_type: str) -> Optional[Dict]:
    """
    Get fragment properties from CSV file, return None if not found or error occurs
    """
    try:
        # Load data from CSV
        df = load_fragment_data_from_csv()
        if df.empty:
            return None
            
        # Query for properties with case-insensitive name matching
        result = df[df['name'].str.lower() == fragment_name.lower()]
        if not result.empty:
            # Convert numeric values to their proper types
            props = {
                'molecular_weight': float(result['molecular_weight'].iloc[0]) if not pd.isna(result['molecular_weight'].iloc[0]) else 0.0,
                'density': float(result['density'].iloc[0]) if not pd.isna(result['density'].iloc[0]) else None,
                'specific_heat_capacity': float(result['specific_heat_capacity'].iloc[0]) if not pd.isna(result['specific_heat_capacity'].iloc[0]) else None,
                'hydrogen_bond_donor_count': int(result['hydrogen_bond_donor_count'].iloc[0]) if not pd.isna(result['hydrogen_bond_donor_count'].iloc[0]) else 0,
                'hydrogen_bond_acceptor_count': int(result['hydrogen_bond_acceptor_count'].iloc[0]) if not pd.isna(result['hydrogen_bond_acceptor_count'].iloc[0]) else 0,
                'rotatable_bond_count': int(result['rotatable_bond_count'].iloc[0]) if not pd.isna(result['rotatable_bond_count'].iloc[0]) else 0,
                'charge': int(result['charge'].iloc[0]) if not pd.isna(result['charge'].iloc[0]) else 0,
                'heavy_atom_count': int(result['heavy_atom_count'].iloc[0]) if not pd.isna(result['heavy_atom_count'].iloc[0]) else 0,
                'topological_polar_surface_area': float(result['tpsa'].iloc[0]) if not pd.isna(result['tpsa'].iloc[0]) else 0.0,
                'smiles': result['smiles'].iloc[0] if not pd.isna(result['smiles'].iloc[0]) else None
            }
            return props
            
        print(f"  No properties found for {fragment_name}")
        return None
        
    except Exception as e:
        print(f"Error getting fragment properties: {str(e)}")
        return None

def is_in_il_thermo(il_name):
    """Check if ionic liquid exists in IL Thermo database."""
    base_url = "https://ilthermo.boulder.nist.gov/ILT2/ilsearch"
    
    # Standardize the name: ensure proper spacing and lowercase
    search_name = ' '.join(il_name.lower().split())
    
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36',
        'Accept': 'application/json'
    }
    
    try:
        params = {'cmp': search_name, 'orderby': 'T', 'output': 'json'}
        response = requests.get(base_url, params=params, headers=headers)
        response.raise_for_status()
        
        # Parse JSON response
        data = response.json()
        
        # Check if we got any results
        if 'res' in data and data['res']:
            return True
        return False
            
    except requests.RequestException as e:
        print(f"Error accessing IL Thermo for {search_name}: {e}")
        return False
    except ValueError as e:
        print(f"Error parsing response for {search_name}: {e}")
        return False

if __name__ == "__main__":
    # Test database connection
    pass