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
    try:
        # Load environment variables
        load_dotenv()
        
        # Database connection parameters
        db_params = {
            'host': os.getenv('DB_HOST'),
            'user': os.getenv('DB_USER'),
            'password': os.getenv('DB_PASSWORD'),
            'db': os.getenv('DB_NAME'),
            'charset': 'utf8mb4',
            'cursorclass': pymysql.cursors.DictCursor  # Use DictCursor to return results as dictionaries
        }
        
        # Connect to database
        connection = pymysql.connect(**db_params)
        cursor = connection.cursor()
        
        return connection, cursor
        
    except Exception as e:
        print(f"Database connection error: {str(e)}")
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
        query = "SELECT molecular_weight FROM IonicLiquids WHERE name = %s AND type = %s"
        cursor.execute(query, (fragment_name, fragment_type))
        result = cursor.fetchone()
        
        if result and result['molecular_weight']:
            return float(result['molecular_weight'])
        return 0.0
        
    except Exception as e:
        print(f"Error getting molecular weight: {e}")
        return 0.0
    finally:
        if conn:
            conn.close()

def get_fragment_properties(fragment_name: str, fragment_type: str) -> Optional[Dict]:
    """
    Get fragment properties from database, return None if not found or error occurs
    """
    try:
        # Connect to database
        conn, cursor = connect_to_database()
        if not conn:
            return None
            
        # Query for properties with case-insensitive type matching
        sql = """
            SELECT molecular_weight, density, specific_heat_capacity,
                   hydrogen_bond_donor_count, hydrogen_bond_acceptor_count,
                   rotatable_bond_count, charge, heavy_atom_count,
                   tpsa, smiles
            FROM IonicLiquids
            WHERE name = %s AND LOWER(type) = LOWER(%s)
        """
        
        cursor.execute(sql, (fragment_name, fragment_type))
        result = cursor.fetchone()
        
        if result:
            # Convert numeric values to their proper types
            props = {
                'molecular_weight': float(result['molecular_weight']) if result['molecular_weight'] else 0.0,
                'density': float(result['density']) if result['density'] else None,
                'specific_heat_capacity': float(result['specific_heat_capacity']) if result['specific_heat_capacity'] else None,
                'hydrogen_bond_donor_count': int(result['hydrogen_bond_donor_count']) if result['hydrogen_bond_donor_count'] else 0,
                'hydrogen_bond_acceptor_count': int(result['hydrogen_bond_acceptor_count']) if result['hydrogen_bond_acceptor_count'] else 0,
                'rotatable_bond_count': int(result['rotatable_bond_count']) if result['rotatable_bond_count'] else 0,
                'charge': int(result['charge']) if result['charge'] else 0,
                'heavy_atom_count': int(result['heavy_atom_count']) if result['heavy_atom_count'] else 0,
                'topological_polar_surface_area': float(result['tpsa']) if result['tpsa'] else 0.0,
                'smiles': result['smiles']
            }
            return props
        
        # If not in database, try to find SMILES in fragments data
        matching_fragments = [f for f in fragments if f['fragment_type'] == fragment_type and f['name'] == fragment_name]
        if matching_fragments:
            frag = matching_fragments[0]
            print(f"  No database entry for {fragment_name}, calculating with RDKit...")
            return get_rdkit_properties(frag['smiles'])
            
        print(f"  No properties found for {fragment_name}")
        return None
        
    except Exception as e:
        print(f"Database error: {str(e)}")
        return None
    finally:
        if cursor:
            cursor.close()
        if conn:
            conn.close()

def get_db_properties(cursor, name):
    """Get properties from database for a given fragment"""
    try:
        # Using the new snake_case column names
        query = """
        SELECT name, smiles, type, molecular_weight, density, 
               specific_heat_capacity, toxicity_rating, ld50
        FROM IonicLiquids 
        WHERE name = %s
        """
        cursor.execute(query, (name,))
        result = cursor.fetchone()
        
        if result:
            # Convert numeric values to float where applicable
            return {
                'name': result['name'],
                'smiles': result['smiles'],
                'type': result['type'],
                'molecular_weight': float(result['molecular_weight']) if result['molecular_weight'] else None,
                'density': float(result['density']) if result['density'] else None,
                'specific_heat_capacity': float(result['specific_heat_capacity']) if result['specific_heat_capacity'] else None,
                'toxicity_rating': result['toxicity_rating'],
                'ld50': result['ld50']
            }
        return None
        
    except Exception as e:
        print(f"Database error: {str(e)}")
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
    conn, cursor = connect_to_database()
    if conn:
        conn.close()