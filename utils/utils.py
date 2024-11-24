# il_naming.py - functions to standardize and generate IL names 
import re
import sys
import os
import requests
import pymysql
from bs4 import BeautifulSoup
from dotenv import load_dotenv
from typing import Dict, Optional

# Add project root to Python path
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path:
    sys.path.append(project_root)

from utils.rdkit_utils import get_rdkit_properties
from models.shortList_frag import fragments

__all__ = ['standardize_il_name', 'generate_il_name', 'get_molecular_weight', 'get_fragment_properties', 'is_in_il_thermo']

def connect_to_database():
    """Connect to MySQL database using environment variables"""
    load_dotenv()
    
    try:
        connection = pymysql.connect(
            host=os.getenv('DB_HOST'),
            user=os.getenv('DB_USER'),
            password=os.getenv('DB_PASSWORD'),
            database=os.getenv('DB_NAME')
        )
        cursor = connection.cursor()
        return connection, cursor
        
    except pymysql.Error as e:
        print(f"Error connecting to MySQL database: {e}")
        return None, None

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
    """Get molecular weight for a fragment from the database"""
    try:
        conn, cursor = connect_to_database()
        if not conn:
            return 0.0

        # Standardize the fragment name
        fragment_name = standardize_il_name(fragment_name)
        
        # Query the IonicLiquids table
        query = "SELECT MolecularWeight FROM IonicLiquids WHERE Name = %s AND Type = %s"
        cursor.execute(query, (fragment_name, fragment_type))
        result = cursor.fetchone()
        
        if result and result[0]:
            return float(result[0])
        return 0.0
        
    except Exception as e:
        print(f"Error getting molecular weight: {e}")
        return 0.0
    finally:
        if conn:
            conn.close()

def get_fragment_properties(fragment_name: str, fragment_type: str) -> Optional[Dict]:
    """
    Get fragment properties from database, falling back to RDKit if not found
    Args:
        fragment_name: Name of the fragment
        fragment_type: Type of fragment (cation, anion, or alkyl_chain)
    Returns:
        Dictionary of properties or None if not found
    """
    # First try database
    props = get_db_properties(fragment_name, fragment_type)
    if props and all(props.get(key) is not None for key in ['molecular_weight', 'heavy_atoms']):
        print(f"Found complete database properties for {fragment_name}")
        return props
        
    print(f"No complete database properties for {fragment_name}, trying RDKit...")
    
    # Try RDKit
    # Get SMILES for the fragment
    smiles = None
    fragment_name_lower = fragment_name.lower()
    fragment_type_lower = fragment_type.lower()
    
    for frag in fragments:
        if (frag['name'].lower() == fragment_name_lower and 
            frag['fragment_type'].lower() == fragment_type_lower):
            smiles = frag['smiles']
            print(f"Found SMILES for {fragment_name}: {smiles}")
            break
    
    if not smiles:
        print(f"No SMILES found for {fragment_name}")
        return None
        
    # Get RDKit properties
    rdkit_props = get_rdkit_properties(smiles)
    if not rdkit_props:
        print(f"Failed to calculate RDKit properties for {fragment_name}")
        return None
        
    # If we have partial database properties, merge them with RDKit properties
    if props:
        rdkit_props.update(props)
        print(f"Merged database and RDKit properties for {fragment_name}")
    
    return rdkit_props

def get_db_properties(fragment_name: str, fragment_type: str) -> Optional[Dict]:
    """Get fragment properties from database"""
    try:
        # Connect to database
        connection, cursor = connect_to_database()
        if not connection or not cursor:
            return None
            
        try:
            # Get properties from database
            sql = """
                SELECT 
                    MolecularWeight, HeavyAtomCount, RotatableBondCount,
                    HydrogenBondDonorCount, HydrogenBondAcceptorCount,
                    Charge, SMILES
                FROM Fragments
                WHERE LOWER(Name) = LOWER(%s) AND LOWER(Type) = LOWER(%s)
            """
            cursor.execute(sql, (fragment_name, fragment_type))
            result = cursor.fetchone()
            
            if result:
                return {
                    'molecular_weight': float(result[0]) if result[0] is not None else None,
                    'heavy_atoms': int(result[1]) if result[1] is not None else None,
                    'rotatable_bonds': int(result[2]) if result[2] is not None else None,
                    'h_bond_donors': int(result[3]) if result[3] is not None else None,
                    'h_bond_acceptors': int(result[4]) if result[4] is not None else None,
                    'charge': int(result[5]) if result[5] is not None else None,
                    'smiles': result[6]
                }
                
        finally:
            cursor.close()
            connection.close()
            
    except Exception as e:
        print(f"Database error: {e}")
        
    return None

def is_in_il_thermo(il_name):
    """Check if ionic liquid exists in IL Thermo database."""
    base_url = "https://ilthermo.boulder.nist.gov/ILT2/ilsearch"
    search_name = ' '.join(il_name.lower().split())
    params = {'cmp': search_name, 'orderby': 'T', 'output': 'json'}
    
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36',
        'Accept': 'application/json'
    }
    
    try:
        response = requests.get(base_url, params=params, headers=headers)
        response.raise_for_status()
        
        # Parse JSON response
        data = response.json()
        
        # Check if we got any results
        if 'res' in data and data['res']:
            return True
        else:
            return False
            
    except requests.RequestException as e:
        print(f"Error accessing IL Thermo: {e}")
        return False
    except ValueError as e:
        print(f"Error parsing response: {e}")
        return False

if __name__ == "__main__":
    # Test database connection
    conn, cursor = connect_to_database()
    if conn:
        conn.close()