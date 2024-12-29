fragments = [
    # Cations (Common cation bases in research and industry)
    {'name': 'Imidazolium', 'smiles': '[n+]1cc[nH]c1', 'fragment_type': 'cation'},  # Imidazolium cation with explicit charge
    {'name': 'Pyridinium', 'smiles': '[n+]1ccccc1', 'fragment_type': 'cation'},  # Pyridinium cation with explicit charge
    {'name': 'Ammonium', 'smiles': '[NH4+]', 'fragment_type': 'cation'},  # Ammonium cation
    {'name': 'Phosphonium', 'smiles': '[PH4+]', 'fragment_type': 'cation'},  # Phosphonium cation
    
    # Anions (Common anions in research and industry)
    {'name': 'Tetrafluoroborate', 'smiles': 'F[B-](F)(F)F', 'fragment_type': 'anion'},  # [BF4]
    {'name': 'Bis(trifluoromethanesulfonyl)imide', 'smiles': 'C(F)(F)(F)S(=O)(=O)[N-]S(=O)(=O)C(F)(F)F', 'fragment_type': 'anion'},  # [NTf2]
    {'name': 'Trifluoromethanesulfonate', 'smiles': 'C(F)(F)(F)S(=O)(=O)[O-]', 'fragment_type': 'anion'},  # [OTf]
    {'name': 'Dicyanamide', 'smiles': '[N-](C#N)C#N', 'fragment_type': 'anion'},  # [DCA]
    
    # Alkyl Chains (Common alkyl substituents)
    {'name': 'Methyl', 'smiles': 'C', 'fragment_type': 'alkyl_chain'},  # -CH3
    {'name': 'Ethyl', 'smiles': 'CC', 'fragment_type': 'alkyl_chain'},  # -CH2CH3
    {'name': 'Propyl', 'smiles': 'CCC', 'fragment_type': 'alkyl_chain'},  # -CH2CH2CH3
    {'name': 'Butyl', 'smiles': 'CCCC', 'fragment_type': 'alkyl_chain'},  # -CH2CH2CH2CH3
]