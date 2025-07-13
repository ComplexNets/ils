fragments = [
    # Cations (Common cation bases in research and industry)
    {'name': 'Imidazolium', 'smiles': '[n+]1cc[nH]c1', 'fragment_type': 'cation'},  # Imidazolium cation with explicit charge
    {'name': 'Ammonium', 'smiles': '[NH4+]', 'fragment_type': 'cation'},  # Ammonium cation
    
    # Anions (Common anions in research and industry)
    {'name': 'Tetrafluoroborate', 'smiles': 'F[B-](F)(F)F', 'fragment_type': 'anion'},  # [BF4]
    {'name': 'Bis(trifluoromethanesulfonyl)imide', 'smiles': 'C(F)(F)(F)S(=O)(=O)[N-]S(=O)(=O)C(F)(F)F', 'fragment_type': 'anion'},  # [NTf2]
    
    # Alkyl Chains (Common alkyl substituents)
    {'name': 'Methyl', 'smiles': 'C', 'fragment_type': 'alkyl_chain'},  # -CH3
    {'name': 'Ethyl', 'smiles': 'CC', 'fragment_type': 'alkyl_chain'},  # -CH2CH3

    # === Functional Groups ===
    {'name': 'Ethoxyethyl', 'smiles': 'CCOCC', 'fragment_type': 'functional_group'},

    # === Anions ===
    {'name': 'Bis(trifluoromethanesulfonyl)imide', 'smiles': 'C(F)(F)(F)S(=O)(=O)[N-]S(=O)(=O)C(F)(F)F', 'fragment_type': 'anion'}
]