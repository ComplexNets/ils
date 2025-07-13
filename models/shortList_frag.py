fragments = [

    # === Cation Cores ===
    {'name': 'Imidazolium', 'smiles': '[n+]1ccn(C)c1', 'fragment_type': 'cation'},
    {'name': 'Phosphonium', 'smiles': '[P+]', 'fragment_type': 'cation'},
    {'name': 'Piperidinium', 'smiles': '[N+]1(CCCCCC1)', 'fragment_type': 'cation'},

    # === Alkyl Chains ===
    {'name': 'Methyl', 'smiles': 'C', 'fragment_type': 'alkyl_chain'},
    {'name': 'Butyl', 'smiles': 'CCCC', 'fragment_type': 'alkyl_chain'},
    {'name': 'Pentyl', 'smiles': 'CCCCC', 'fragment_type': 'alkyl_chain'},
    {'name': 'Nonyl', 'smiles': 'CCCCCCCCC', 'fragment_type': 'alkyl_chain'},
    {'name': 'Octyl', 'smiles': 'CCCCCCCC', 'fragment_type': 'alkyl_chain'},
    {'name': 'Decyl', 'smiles': 'CCCCCCCCCC', 'fragment_type': 'alkyl_chain'},
    {'name': 'Oleyl', 'smiles': 'CCCCCCCC=CCCCCCCCCC', 'fragment_type': 'alkyl_chain'},  # cis-oleyl

    # === Functional Groups ===
    {'name': 'Ethoxyethyl', 'smiles': 'CCOCC', 'fragment_type': 'functional_group'},

    # === Anions ===
    {'name': 'Bis(trifluoromethanesulfonyl)imide', 'smiles': 'C(F)(F)(F)S(=O)(=O)[N-]S(=O)(=O)C(F)(F)F', 'fragment_type': 'anion'}
]