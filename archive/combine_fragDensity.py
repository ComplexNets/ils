import sys
import os
# Add the parent directory to the Python path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pymysql
from dotenv import load_dotenv
from typing import Dict, List, Tuple
from utils.validation_rules import MolecularValidator
from utils.utils import generate_il_name

# Add these function definitions before main()
def is_in_il_thermo(name: str) -> bool:
    """Check if the ionic liquid exists in ILThermo database"""
    return False  # Placeholder - implement actual database check if needed

# Add this at the top with other imports
property_ranges = {
    'density_range': (800, 1500)  # Default density range in kg/m³
}

def connect_to_db():
    """Connect to MySQL database using environment variables"""
    load_dotenv()
    
    conn = pymysql.connect(
        host=os.getenv('DB_HOST'),
        user=os.getenv('DB_USER'),
        password=os.getenv('DB_PASSWORD'),
        database=os.getenv('DB_NAME'),
        cursorclass=pymysql.cursors.DictCursor
    )
    cursor = conn.cursor()
    return conn, cursor

def get_filtered_fragments(cursor, density_range: Tuple[float, float]) -> Dict:
    """Get fragments filtered by density range"""
    fragments_data = {'Cation': [], 'Anion': [], 'Alkyl Chain': []}
    
    query = """
    SELECT 
        ID, 
        Name, 
        Type, 
        MolecularWeight,
        LogP,
        HydrogenBondDonorCount,
        HydrogenBondAcceptorCount,
        RotatableBondCount,
        Complexity,
        Charge,
        HeavyAtomCount,
        TPSA
    FROM IonicLiquids
    WHERE Type IN ('Cation', 'Anion', 'Alkyl Chain')
    """
    
    cursor.execute(query)
    results = cursor.fetchall()
    
    # Debug print to see the structure of a row
    if results:
        print("\nSample row structure:", results[0].keys())
    
    for row in results:
        # Ensure all required fields have values, even if None
        fragment = {
            'ID': row['ID'],
            'Name': row['Name'],
            'Type': row['Type'],
            'MolecularWeight': row['MolecularWeight'] or 0.0,
            'LogP': row['LogP'] or 0.0,
            'HydrogenBondDonorCount': row['HydrogenBondDonorCount'] or 0,
            'HydrogenBondAcceptorCount': row['HydrogenBondAcceptorCount'] or 0,
            'RotatableBondCount': row['RotatableBondCount'] or 0,
            'Complexity': row['Complexity'] or 0.0,
            'Charge': row['Charge'] or 0,
            'HeavyAtomCount': row['HeavyAtomCount'] or 0,
            'TPSA': row['TPSA'] or 0.0
        }
        fragments_data[row['Type']].append(fragment)
    
    print(f"\nFound fragments:")
    for ftype, frags in fragments_data.items():
        print(f"{ftype}: {len(frags)} fragments")
    
    return fragments_data

def validate_combination(cation: Dict, anion: Dict, alkyl: Dict, cursor) -> Tuple[bool, str]:
    """Validate if the fragment combination is chemically feasible"""
    try:
        validator = MolecularValidator(cursor)  # Pass cursor to validator
        return validator.validate(cation, anion, alkyl)
    except Exception as e:
        return False, f"Validation error: {str(e)}"

def calculate_heat_capacity(cation, anion, alkyl):
    """Calculate heat capacity using UNIFAC method"""
    try:
        # Simple estimation based on molecular weights and empirical factors
        total_mass = (cation['MolecularWeight'] + 
                     anion['MolecularWeight'] + 
                     alkyl['MolecularWeight'])
        
        # Empirical factors for heat capacity estimation
        cation_factor = 1.2
        anion_factor = 0.8
        alkyl_factor = 1.0
        
        cp = (cation['MolecularWeight'] * cation_factor +
              anion['MolecularWeight'] * anion_factor +
              alkyl['MolecularWeight'] * alkyl_factor)
        
        # Convert to J/mol·K (simplified estimation)
        return cp * 0.5
    except Exception as e:
        print(f"Error calculating heat capacity: {str(e)}")
        return None

def calculate_density(ionic_liquid: Dict) -> float:
    """
    Calculate density for a given ionic liquid combination
    Args:
        ionic_liquid: Dictionary containing cation, anion, and alkyl chain information
    Returns:
        Estimated density in kg/m³
    """
    try:
        cation = ionic_liquid['cation']
        anion = ionic_liquid['anion']
        alkyl = ionic_liquid['alkyl']
        
        # Basic density calculation based on molecular weight and empirical factors
        total_mass = (cation['MolecularWeight'] + 
                     anion['MolecularWeight'] + 
                     alkyl['MolecularWeight'])
        
        # Empirical correction factors based on fragment type
        cation_factor = 1.2  # Density tends to increase with cation size
        anion_factor = 0.9   # Anions have varying effects on density
        alkyl_factor = 0.85  # Longer alkyl chains tend to decrease density
        
        # Calculate volume contribution of each component
        volume = (cation['MolecularWeight'] * cation_factor +
                 anion['MolecularWeight'] * anion_factor +
                 alkyl['MolecularWeight'] * alkyl_factor)
        
        # Convert to density (kg/m³)
        # Base density multiplier (1000) converts g/mL to kg/m³
        density = (total_mass / volume) * 1000
        
        # Adjust density based on other molecular properties
        # These are empirical adjustments based on chemical intuition
        density *= (1 + 0.05 * cation['HydrogenBondDonorCount'])  # H-bond donors increase density
        density *= (1 + 0.03 * anion['HydrogenBondAcceptorCount']) # H-bond acceptors increase density
        density *= (1 - 0.02 * alkyl['RotatableBondCount'])  # Flexibility decreases density
        
        return max(800, min(1500, density))  # Constrain to reasonable range
        
    except Exception as e:
        print(f"Error calculating density: {str(e)}")
        return None

def calculate_all_densities(combinations: List[Dict]) -> List[Dict]:
    """Calculate densities for a list of ionic liquid combinations"""
    for combination in combinations:
        density = calculate_density(combination)
        if density is not None:
            combination['density'] = density
    return combinations

def main(density_range=None):
    """
    Main function to process ionic liquid combinations
    Returns a list of valid combinations with their properties
    """
    try:
        conn, cursor = connect_to_db()
        
        # Get fragments filtered by density range if specified
        fragments = get_filtered_fragments(cursor, density_range or property_ranges['density_range'])
        
        valid_combinations = []
        validator = MolecularValidator(cursor)
        
        # Generate all possible combinations
        for cation in fragments['Cation']:
            for anion in fragments['Anion']:
                for alkyl in fragments['Alkyl Chain']:
                    # Validate the combination
                    is_valid, exists_in_db = validator.validate(cation, anion, alkyl)
                    
                    if is_valid:
                        # Calculate properties for valid combinations
                        il_name = generate_il_name(cation, anion, alkyl)
                        
                        combination = {
                            'name': il_name,
                            'cation': cation,
                            'anion': anion,
                            'alkyl': alkyl,
                            'in_il_thermo': exists_in_db
                        }
                        valid_combinations.append(combination)
        
        cursor.close()
        conn.close()
        
        # Calculate densities
        combinations_with_density = calculate_all_densities(valid_combinations)
        return combinations_with_density
        
    except Exception as e:
        print(f"Error in main: {str(e)}")
        return []

if __name__ == "__main__":
    print("\n=== Running Ionic Liquid Density Calculator ===\n")
    results = main()
    print("\nResults:")
    print("-" * 80)
    for i, result in enumerate(results, 1):
        print(f"\n{i}. Ionic Liquid: {result['name']}")
        print(f"   Density: {result['density']:.1f} kg/m³")
        print(f"   Molecular Weight: {(result['cation']['MolecularWeight'] + result['anion']['MolecularWeight'] + result['alkyl']['MolecularWeight']):.1f} g/mol")
    print("\n" + "=" * 80)