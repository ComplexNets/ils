from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
from typing import Dict, Optional

def get_rdkit_properties(smiles: str) -> Optional[Dict]:
    """
    Calculate molecular properties using RDKit
    Args:
        smiles: SMILES string of the molecule
    Returns:
        Dictionary of properties or None if calculation fails
    """
    try:
        # Create RDKit molecule
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"Failed to parse SMILES: {smiles}")
            return None
            
        # Add hydrogens
        mol = Chem.AddHs(mol)
        
        # Calculate 3D coordinates (needed for some descriptors)
        try:
            AllChem.EmbedMolecule(mol, randomSeed=42)
        except:
            print("Warning: 3D coordinate generation failed, some properties may be approximate")
        
        # Calculate properties with standardized names and types
        props = {
            'molecular_weight': float(Descriptors.ExactMolWt(mol)),
            'heavy_atom_count': int(mol.GetNumHeavyAtoms()),
            'rotatable_bond_count': int(Descriptors.NumRotatableBonds(mol)),
            'hydrogen_bond_donor_count': int(Descriptors.NumHDonors(mol)),
            'hydrogen_bond_acceptor_count': int(Descriptors.NumHAcceptors(mol)),
            'charge': int(Chem.GetFormalCharge(mol)),
            'aromatic_ring_count': int(Descriptors.NumAromaticRings(mol)),
            'topological_polar_surface_area': float(Descriptors.TPSA(mol)),
            'smiles': Chem.MolToSmiles(mol)  # Canonical SMILES
        }
        
        print(f"\nRDKit properties for {smiles}:")
        for prop, value in props.items():
            print(f"  {prop}: {value}")
            
        return props
        
    except Exception as e:
        print(f"Error calculating RDKit properties: {e}")
        return None

def test_rdkit():
    """Test RDKit property calculations"""
    # Test with EMIM cation
    emim_smiles = "CCn1ccn(C)c1"  # EMIM cation
    print("\nTesting EMIM cation:")
    emim_props = get_rdkit_properties(emim_smiles)
    
    # Test with BF4 anion
    bf4_smiles = "[B-](F)(F)(F)F"  # BF4 anion
    print("\nTesting BF4 anion:")
    bf4_props = get_rdkit_properties(bf4_smiles)

if __name__ == "__main__":
    test_rdkit()
