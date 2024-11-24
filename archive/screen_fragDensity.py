import pymysql
from dotenv import load_dotenv
import os
from property_input import property_ranges
from shortList_frag import fragments

def connect_to_database():
    """Connect to the remote database using pymysql and environment variables"""
    load_dotenv()
    try:
        connection = pymysql.connect(
            host=os.getenv('DB_HOST'),
            user=os.getenv('DB_USER'),
            password=os.getenv('DB_PASSWORD'),
            database=os.getenv('DB_NAME'),
            cursorclass=pymysql.cursors.DictCursor
        )
        return connection
    except Exception as e:
        print(f"Error connecting to database: {e}")
        return None

def calculate_molecular_volume(molecular_weight, heavy_atoms, fragment_type):
    """Calculate molecular volume if not available in database"""
    # Approximate volume per heavy atom (in Å³)
    if fragment_type == 'Cation':
        vol_per_atom = 25.0
    elif fragment_type == 'Anion':
        vol_per_atom = 23.0
    else:  # Alkyl Chain
        vol_per_atom = 21.0
        
    # Calculate volume in cm³/mol (converting from Å³)
    volume = (heavy_atoms * vol_per_atom) * 1e-24 * 6.022e23
    return volume

def calculate_fragment_density(molecular_weight, molecular_volume):
    """Calculate density for a single fragment"""
    # Convert molecular weight from g/mol to kg/mol
    mw_kg = molecular_weight / 1000
    
    # molecular_volume should be in cm³/mol
    # Convert to m³/mol
    volume_m3 = molecular_volume / 1000000
    
    # Calculate density in kg/m³
    density = mw_kg / volume_m3
    
    return density

def update_density_in_database(cursor, fragment_name, fragment_type, estimated_density):
    """Update the estimated density in the database"""
    query = """
        UPDATE IonicLiquids 
        SET est_density = %s 
        WHERE Name = %s AND Type = %s
    """
    cursor.execute(query, (estimated_density, fragment_name, fragment_type))

def screen_density(cursor):
    """Get density information for all fragments and update database"""
    fragments_with_density = {
        'Cation': [],
        'Anion': [],
        'Alkyl Chain': []
    }
    
    # First, add est_density column if it doesn't exist
    try:
        cursor.execute("""
            ALTER TABLE IonicLiquids 
            ADD COLUMN IF NOT EXISTS est_density FLOAT
        """)
    except Exception as e:
        print(f"Error adding est_density column: {e}")
    
    allowed_fragments = {
        'Cation': [frag['name'] for frag in fragments if frag['fragment_type'] == 'Cation'],
        'Anion': [frag['name'] for frag in fragments if frag['fragment_type'] == 'Anion'],
        'Alkyl Chain': [frag['name'] for frag in fragments if frag['fragment_type'] == 'Alkyl Chain']
    }
    
    for frag_type, names in allowed_fragments.items():
        if not names:
            continue
            
        placeholders = ', '.join(['%s'] * len(names))
        query = f"""
            SELECT Name, Type, MolecularWeight, HeavyAtomCount, SMILES, est_density
            FROM IonicLiquids
            WHERE Type = %s 
            AND Name IN ({placeholders})
        """
        
        cursor.execute(query, (frag_type, *names))
        results = cursor.fetchall()
        
        for result in results:
            # Calculate molecular volume
            mol_volume = calculate_molecular_volume(
                result['MolecularWeight'],
                result['HeavyAtomCount'],
                result['Type']
            )
            
            # Calculate estimated density if not already in database
            if result['est_density'] is None:
                est_density = calculate_fragment_density(
                    result['MolecularWeight'],
                    mol_volume
                )
                # Update database with estimated density
                update_density_in_database(
                    cursor, 
                    result['Name'], 
                    result['Type'], 
                    est_density
                )
            else:
                est_density = result['est_density']
            
            fragments_with_density[result['Type']].append({
                'name': result['Name'],
                'molecular_weight': result['MolecularWeight'],
                'heavy_atoms': result['HeavyAtomCount'],
                'smiles': result['SMILES'],
                'molecular_volume': mol_volume,
                'density': est_density,
                'fragment_type': result['Type']
            })
    
    return fragments_with_density

def main():
    connection = connect_to_database()
    if not connection:
        return
    
    cursor = connection.cursor()
    
    try:
        fragments_with_density = screen_density(cursor)
        
        # Print results
        for frag_type, frags in fragments_with_density.items():
            print(f"\n{frag_type}s ({len(frags)} fragments):")
            for frag in frags:
                print(f"  - {frag['name']}")
                print(f"    Estimated Density: {frag['density']:.1f} kg/m³")
                print(f"    Molecular Weight: {frag['molecular_weight']:.1f} g/mol")
                print(f"    Molecular Volume: {frag['molecular_volume']:.1f} cm³/mol")
    
    except Exception as e:
        print(f"Error during density processing: {e}")
    finally:
        cursor.close()
        connection.close()

if __name__ == "__main__":
    main()
