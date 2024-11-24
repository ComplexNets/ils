fragments = [
    # Cations (Most common imidazolium and pyridinium cations)
    {'name': 'Ethylimidazolium', 'smiles': 'CCn1cc[n+]c1', 'fragment_type': 'cation'},  # [EMIM]
    {'name': '1-Butylpyridinium', 'smiles': 'CCCC[n+]1ccccc1', 'fragment_type': 'cation'},  # [BPy]

    # Anions (Most common anions in research and industry)
    {'name': 'Tetrafluoroborate', 'smiles': 'F[B-](F)(F)F', 'fragment_type': 'anion'},  # [BF4]
    {'name': 'Bis(trifluoromethanesulfonyl)imide', 'smiles': 'C(F)(F)(F)S(=O)(=O)[N-]S(=O)(=O)C(F)(F)F', 'fragment_type': 'anion'},  # [NTf2]

    # Alkyl Chains (Most common chain lengths)
    {'name': 'Ethyl', 'smiles': 'CC', 'fragment_type': 'alkyl_chain'},
    {'name': 'Butyl', 'smiles': 'CCCC', 'fragment_type': 'alkyl_chain'},
]