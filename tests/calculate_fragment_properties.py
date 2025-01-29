import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from models.shortList_frag import fragments

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
    # Create fragment_data directory if it doesn't exist
    os.makedirs('fragment_data', exist_ok=True)
    
    # Calculate properties for all fragments
    data = []
    for fragment in fragments:
        props = calculate_properties(fragment['smiles'])
        if props:
            row = {
                'name': fragment['name'],
                'smiles': fragment['smiles'],
                'fragment_type': fragment['fragment_type'],
                **props
            }
            data.append(row)
            print(f"Calculated properties for {fragment['name']}:")
            for key, value in props.items():
                print(f"  {key}: {value}")
        else:
            print(f"Failed to calculate properties for {fragment['name']}")
    
    # Create DataFrame and save to CSV
    df = pd.DataFrame(data)
    csv_path = 'fragment_data/fragment_properties.csv'
    df.to_csv(csv_path, index=False)
    print(f"\nSaved properties to {csv_path}")
    print("\nFragment counts by type:")
    print(df['fragment_type'].value_counts())

if __name__ == "__main__":
    main()
