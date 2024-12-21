
fragments = [
    # Cations
    {'name': 'Ammonium', 'smiles': '[NH4+]', 'fragment_type': 'cation'},
    {'name': 'Phosphonium', 'smiles': '[PH4+]', 'fragment_type': 'cation'},
    {'name': 'Imidazolium', 'smiles': 'C1=CN=C[CNH+]1', 'fragment_type': 'cation'},
    {'name': 'Pyridinium', 'smiles': 'C1=CC=[NH+]C=C1', 'fragment_type': 'cation'},
    {'name': 'Methylimidazolium', 'smiles': 'C[nH+]1cccnc1', 'fragment_type': 'cation'},
    {'name': 'Tetramethylammonium', 'smiles': '[N+](C)(C)(C)C', 'fragment_type': 'cation'},
    {'name': '1-Butyl-3-methylimidazolium', 'smiles': 'CCCC[nH+]1cccnc1C', 'fragment_type': 'cation'},
    {'name': '1-Octyl-3-methylimidazolium', 'smiles': 'CCCCCCCC[nH+]1cccnc1C', 'fragment_type': 'cation'},
    {'name': '1-Methylpyridinium', 'smiles': 'C[nH+]1ccccc1', 'fragment_type': 'cation'},
    {'name': '1-Propylpyridinium', 'smiles': 'CCC[nH+]1ccccc1', 'fragment_type': 'cation'},
    {'name': '1-Hexylpyridinium', 'smiles': 'CCCCCC[nH+]1ccccc1', 'fragment_type': 'cation'},

    # Anions
    {'name': 'Chloride', 'smiles': '[Cl-]', 'fragment_type': 'anion'},
    {'name': 'Bromide', 'smiles': '[Br-]', 'fragment_type': 'anion'},
    {'name': 'Iodide', 'smiles': '[I-]', 'fragment_type': 'anion'},
    {'name': 'Dicyanamide', 'smiles': '[N-](C#N)C#N', 'fragment_type': 'anion'},
    {'name': 'Trifluoromethanesulfonate', 'smiles': 'C(F)(F)(F)S(=O)(=O)[O-]', 'fragment_type': 'anion'},
    {'name': 'Methanesulfonate', 'smiles': 'CS(=O)(=O)[O-]', 'fragment_type': 'anion'},
    {'name': 'Tosylate', 'smiles': 'CCc1ccccc1S(=O)(=O)[O-]', 'fragment_type': 'anion'},
    {'name': 'Trifluoroacetate', 'smiles': 'C(F)(F)F C(=O)[O-]', 'fragment_type': 'anion'},
    {'name': 'Phosphate', 'smiles': 'O[P](=O)([O-])[O-]', 'fragment_type': 'anion'},
    {'name': 'Bis(trifluoromethanesulfonyl)imide', 'smiles': 'C(F)(F)(F)S(=O)(=O)N(S(=O)(=O)C(F)(F)F)[-]', 'fragment_type': 'anion'},

    # Alkyl Chains
    {'name': 'Methyl', 'smiles': 'C', 'fragment_type': 'alkyl_chain'},
    {'name': 'Ethyl', 'smiles': 'CC', 'fragment_type': 'alkyl_chain'},
    {'name': 'Propyl', 'smiles': 'CCC', 'fragment_type': 'alkyl_chain'},
    {'name': 'Butyl', 'smiles': 'CCCC', 'fragment_type': 'alkyl_chain'},
    {'name': 'Pentyl', 'smiles': 'CCCCC', 'fragment_type': 'alkyl_chain'},
    {'name': 'Hexyl', 'smiles': 'CCCCCC', 'fragment_type': 'alkyl_chain'}
]