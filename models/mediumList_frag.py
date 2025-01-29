fragments = [
    # Cations
    {'name': 'Ammonium', 'smiles': '[NH4+]', 'fragment_type': 'cation'},
    {'name': 'Phosphonium', 'smiles': '[PH4+]', 'fragment_type': 'cation'},
    {'name': 'Imidazolium', 'smiles': '[n+]1cc[nH]c1', 'fragment_type': 'cation'},
    {'name': 'Pyridinium', 'smiles': '[n+]1ccccc1', 'fragment_type': 'cation'},
    {'name': 'Methylimidazolium', 'smiles': 'C[n+]1cc[nH]c1', 'fragment_type': 'cation'},
    {'name': 'Tetramethylammonium', 'smiles': '[N+](C)(C)(C)C', 'fragment_type': 'cation'},
    {'name': '1-Butyl-3-methylimidazolium', 'smiles': 'CCCC[n+]1cc[nH]c1C', 'fragment_type': 'cation'},
    {'name': '1-Octyl-3-methylimidazolium', 'smiles': 'CCCCCCCC[n+]1cc[nH]c1C', 'fragment_type': 'cation'},
    {'name': '1-Methylpyridinium', 'smiles': 'C[n+]1ccccc1', 'fragment_type': 'cation'},
    {'name': '1-Propylpyridinium', 'smiles': 'CCC[n+]1ccccc1', 'fragment_type': 'cation'},
    {'name': '1-Hexylpyridinium', 'smiles': 'CCCCCC[n+]1ccccc1', 'fragment_type': 'cation'},

    # Anions
    {'name': 'Chloride', 'smiles': '[Cl-]', 'fragment_type': 'anion'},
    {'name': 'Bromide', 'smiles': '[Br-]', 'fragment_type': 'anion'},
    {'name': 'Iodide', 'smiles': '[I-]', 'fragment_type': 'anion'},
    {'name': 'Dicyanamide', 'smiles': '[N-](C#N)C#N', 'fragment_type': 'anion'},
    {'name': 'Trifluoromethanesulfonate', 'smiles': 'FC(F)(F)S(=O)(=O)[O-]', 'fragment_type': 'anion'},
    {'name': 'Methanesulfonate', 'smiles': 'CS(=O)(=O)[O-]', 'fragment_type': 'anion'},
    {'name': 'Tosylate', 'smiles': 'Cc1ccc(cc1)S(=O)(=O)[O-]', 'fragment_type': 'anion'},
    {'name': 'Trifluoroacetate', 'smiles': 'FC(F)(F)C(=O)[O-]', 'fragment_type': 'anion'},
    {'name': 'Phosphate', 'smiles': '[O-]P(=O)([O-])[O-]', 'fragment_type': 'anion'},
    {'name': 'Bis(trifluoromethanesulfonyl)imide', 'smiles': 'FC(F)(F)S(=O)(=O)[N-]S(=O)(=O)C(F)(F)F', 'fragment_type': 'anion'},

    # Alkyl Chains
    {'name': 'Methyl', 'smiles': 'C', 'fragment_type': 'alkyl_chain'},
    {'name': 'Ethyl', 'smiles': 'CC', 'fragment_type': 'alkyl_chain'},
    {'name': 'Propyl', 'smiles': 'CCC', 'fragment_type': 'alkyl_chain'},
    {'name': 'Butyl', 'smiles': 'CCCC', 'fragment_type': 'alkyl_chain'},
    {'name': 'Pentyl', 'smiles': 'CCCCC', 'fragment_type': 'alkyl_chain'},
    {'name': 'Hexyl', 'smiles': 'CCCCCC', 'fragment_type': 'alkyl_chain'},

     # Functional Groups
    {'name': 'Hydroxyl', 'smiles': 'O', 'fragment_type': 'functional_group'},# -OH
    {'name': 'Amino', 'smiles': 'N', 'fragment_type': 'functional_group'}, # -NH2
]