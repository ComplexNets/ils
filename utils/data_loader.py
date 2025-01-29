"""
Utility functions for loading and processing fragment data from CSV.
"""
import pandas as pd
from typing import Dict, Optional

def load_fragment_data() -> Optional[pd.DataFrame]:
    """Load fragment data from CSV"""
    try:
        df = pd.read_csv('fragment_data/autono17_ilselect_db.csv')
        return df
    except Exception as e:
        print(f"Error loading fragment data: {str(e)}")
        return None

def get_filtered_fragments() -> Dict:
    """Get fragments from CSV database"""
    fragments_data = {'cation': [], 'anion': []}
    
    # Load data from CSV
    df = load_fragment_data()
    if df is None:
        return fragments_data
        
    print("\nGetting fragment properties from CSV:")
    for _, row in df.iterrows():
        # Determine fragment type based on charge
        frag_type = 'cation' if row['charge'] > 0 else 'anion'
        
        processed_fragment = {
            'name': row['name'],
            'fragment_type': frag_type,
            'molecular_weight': row['molecular_weight'],
            'log_p': row['log_p'],
            'hydrogen_bond_donor_count': row['hydrogen_bond_donor_count'],
            'hydrogen_bond_acceptor_count': row['hydrogen_bond_acceptor_count'],
            'tpsa': row['tpsa'],
            'rotatable_bond_count': row['rotatable_bond_count'],
            'complexity': row['complexity'],
            'charge': row['charge'],
            'heavy_atom_count': row['heavy_atom_count']
        }
        
        fragments_data[frag_type].append(processed_fragment)
        print(f"Added {row['name']} as {frag_type}")
    
    return fragments_data
