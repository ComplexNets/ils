import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

# Define our functional groups
functional_groups = [
    {
        'name': 'Hydroxyl',
        'smiles': 'O',  # -OH group
        'fragment_type': 'functional_group'
    },
    {
        'name': 'Amino',
        'smiles': 'N',  # -NH2 group
        'fragment_type': 'functional_group'
    },
    {
        'name': 'Cyano',
        'smiles': 'C#N',  # -CN group
        'fragment_type': 'functional_group'
    }
]

def calculate_properties(smiles):
    """Calculate molecular properties using RDKit"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
        
    return {
        'molecular_weight': Descriptors.ExactMolWt(mol),
        'log_p': Descriptors.MolLogP(mol),
        'hydrogen_bond_donor_count': rdMolDescriptors.CalcNumHBD(mol),
        'hydrogen_bond_acceptor_count': rdMolDescriptors.CalcNumHBA(mol),
        'rotatable_bond_count': Descriptors.NumRotatableBonds(mol),
        'complexity': Descriptors.BertzCT(mol),
        'charge': sum(atom.GetFormalCharge() for atom in mol.GetAtoms()),
        'heavy_atom_count': mol.GetNumHeavyAtoms(),
        'tpsa': Descriptors.TPSA(mol)
    }

def main():
    # Read existing CSV
    csv_path = 'fragment_data/autono17_ilselect_db.csv'
    df = pd.read_csv(csv_path)
    
    # Calculate properties for functional groups
    new_rows = []
    next_id = df['id'].max() + 1
    
    for fg in functional_groups:
        props = calculate_properties(fg['smiles'])
        if props:
            row = {
                'id': next_id,
                'name': fg['name'],
                'smiles': fg['smiles'],
                'type': None,
                'fragment_type': fg['fragment_type'],
                **props
            }
            # Fill remaining columns with NULL
            for col in df.columns:
                if col not in row:
                    row[col] = None
                    
            new_rows.append(row)
            next_id += 1
            
            print(f"\nCalculated properties for {fg['name']}:")
            for key, value in props.items():
                print(f"  {key}: {value}")
    
    # Append new rows to DataFrame
    df_new = pd.concat([df, pd.DataFrame(new_rows)], ignore_index=True)
    
    # Save updated DataFrame
    df_new.to_csv(csv_path, index=False)
    print(f"\nUpdated {csv_path} with {len(new_rows)} new functional groups")

if __name__ == "__main__":
    main()
