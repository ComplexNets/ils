import pymysql
from dotenv import load_dotenv
import os
from property_input import property_ranges
from shortList_frag import fragments

def connect_to_db():
    """Establish database connection using environment variables"""
    load_dotenv()
    
    conn = pymysql.connect(
        host=os.getenv('DB_HOST'),
        user=os.getenv('DB_USER'),
        password=os.getenv('DB_PASSWORD'),
        database=os.getenv('DB_NAME')
    )
    return conn, conn.cursor()

def estimate_heat_capacity(cursor, fragment_name):
    """
    Basic estimation of specific heat capacity using:
    - Molecular weight
    - Number of atoms (HeavyAtomCount)
    - Number of rotatable bonds (RotatableBondCount)
    """
    cursor.execute('''
        SELECT MolecularWeight, HeavyAtomCount, RotatableBondCount
        FROM IonicLiquids
        WHERE Name = %s
    ''', (fragment_name,))
    
    result = cursor.fetchone()
    if not result:
        return None
    
    mw, heavy_atoms, rotatable_bonds = result
    
    if all(x is not None for x in [mw, heavy_atoms, rotatable_bonds]):
        # Simple linear estimation
        estimated_capacity = (
            0.5 * mw +              # Weight contribution
            10 * heavy_atoms +      # Atom contribution
            5 * rotatable_bonds +   # Flexibility contribution
            30                      # Base value
        )
        return estimated_capacity
    return None

def update_estimated_heat_capacity(cursor, fragment_name, estimated_capacity):
    """Update the database with estimated heat capacity in the new column"""
    try:
        sql = '''
            UPDATE IonicLiquids 
            SET est_heat_capacity = %s
            WHERE Name = %s 
        '''
        cursor.execute(sql, (estimated_capacity, fragment_name))
        return cursor.rowcount > 0
    except Exception as e:
        print(f"Error updating estimated heat capacity for {fragment_name}: {e}")
        return False

def get_fragments_by_heat_capacity(cursor, heat_capacity_range):
    """
    Get fragments from shortList_frag.py and estimate their heat capacities
    """
    min_capacity, max_capacity = heat_capacity_range
    
    # Initialize dictionary for organized fragments
    fragments_data = {
        'Cation': [],
        'Anion': [],
        'Alkyl Chain': []
    }
    
    print("\nEstimating heat capacities for fragments...")
    
    # Process fragments from shortList_frag.py
    for fragment in fragments:
        frag_type = fragment['fragment_type']
        if frag_type in fragments_data:
            # Get fragment properties from database
            cursor.execute('''
                SELECT MolecularWeight, HeavyAtomCount, RotatableBondCount, 
                       SpecificHeatCapacity, est_heat_capacity
                FROM IonicLiquids
                WHERE Name = %s AND Type = %s
            ''', (fragment['name'], frag_type))
            
            result = cursor.fetchone()
            if result:
                mw, atoms, bonds, actual_cp, est_cp = result
                
                estimated_capacity = estimate_heat_capacity(cursor, fragment['name'])
                if estimated_capacity is not None:
                    # Update database with estimate
                    if update_estimated_heat_capacity(cursor, fragment['name'], estimated_capacity):
                        print(f"Updated {fragment['name']} with estimated heat capacity: {estimated_capacity:.1f} J/mol·K")
                        if actual_cp:
                            print(f"  (Actual heat capacity: {actual_cp:.1f} J/mol·K)")
                    
                    fragments_data[frag_type].append({
                        'name': fragment['name'],
                        'smiles': fragment['smiles'],
                        'heat_capacity': estimated_capacity,
                        'actual_heat_capacity': actual_cp,
                        'molecular_weight': mw,
                        'heavy_atoms': atoms,
                        'rotatable_bonds': bonds
                    })
    
    print(f"\nFragment Analysis for Heat Capacity Range: {min_capacity} to {max_capacity} J/mol·K")
    print("\nTop candidates by type:")
    
    for frag_type, frags in fragments_data.items():
        if frags:
            print(f"\n{frag_type}s:")
            target_capacity = (min_capacity + max_capacity) / 2
            sorted_frags = sorted(frags, key=lambda x: abs(x['heat_capacity'] - target_capacity))
            for frag in sorted_frags[:3]:
                print(f"- {frag['name']}:")
                print(f"  Estimated Heat Capacity: {frag['heat_capacity']:.1f} J/mol·K")
                if frag['actual_heat_capacity']:
                    print(f"  Actual Heat Capacity: {frag['actual_heat_capacity']:.1f} J/mol·K")
                print(f"  MW: {frag['molecular_weight']}, Heavy Atoms: {frag['heavy_atoms']}, Rotatable Bonds: {frag['rotatable_bonds']}")
    
    return fragments_data

def main():
    conn, cursor = connect_to_db()
    
    # Fix the key to match property_input.py
    heat_capacity_range = property_ranges['heat_capacity']  # Changed from 'heat_capacity_range'
    
    fragments = get_fragments_by_heat_capacity(cursor, heat_capacity_range)
    
    if not any(fragments.values()):
        print("\nNo fragments found with sufficient data for heat capacity estimation.")
        print("Please ensure fragments have:")
        print("- Molecular Weight")
        print("- Heavy Atom Count")
        print("- Rotatable Bond Count")
    
    conn.close()

if __name__ == "__main__":
    main() 