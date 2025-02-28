fragments = [
    # Cations
    {'name': 'Ammonium', 'smiles': '[NH4+]', 'fragment_type': 'cation'},
    {'name': 'Phosphonium', 'smiles': '[PH4+]', 'fragment_type': 'cation'},
    {'name': 'Imidazolium', 'smiles': '[n+]1cc[nH]c1', 'fragment_type': 'cation'},
    {'name': 'Pyridinium', 'smiles': '[n+]1ccccc1', 'fragment_type': 'cation'},
    {'name': 'Methylimidazolium', 'smiles': 'C[n+]1cc[nH]c1', 'fragment_type': 'cation'},
    {'name': 'Tetramethylammonium', 'smiles': '[N+](C)(C)(C)C', 'fragment_type': 'cation'},
    {'name': '1-Butyl-3-methylimidazolium', 'smiles': 'CCCC[n+]1cc[nH]c1C', 'fragment_type': 'cation'},
    {'name': '1-Propylpyridinium', 'smiles': 'CCC[n+]1ccccc1', 'fragment_type': 'cation'},

    # Anions
    {'name': 'Chloride', 'smiles': '[Cl-]', 'fragment_type': 'anion'},
    {'name': 'Bromide', 'smiles': '[Br-]', 'fragment_type': 'anion'},
    {'name': 'Trifluoromethanesulfonate', 'smiles': 'FC(F)(F)S(=O)(=O)[O-]', 'fragment_type': 'anion'},
    {'name': 'Trifluoroacetate', 'smiles': 'FC(F)(F)C(=O)[O-]', 'fragment_type': 'anion'},
 
    # Alkyl Chains
    {'name': 'Methyl', 'smiles': 'C', 'fragment_type': 'alkyl_chain'},
    {'name': 'Ethyl', 'smiles': 'CC', 'fragment_type': 'alkyl_chain'},
    {'name': 'Propyl', 'smiles': 'CCC', 'fragment_type': 'alkyl_chain'},
    {'name': 'Butyl', 'smiles': 'CCCC', 'fragment_type': 'alkyl_chain'},

     # Functional Groups
    {'name': 'Hydroxyl', 'smiles': 'O', 'fragment_type': 'functional_group'},# -OH
    {'name': 'Amino', 'smiles': 'N', 'fragment_type': 'functional_group'}, # -NH2
]