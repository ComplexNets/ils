fragments = [
    # Cations (Common imidazolium, pyridinium, and ammonium cations)
    {'name': 'Imidazolium', 'smiles': 'c1cc[nH]c1', 'fragment_type': 'cation'},  # Base imidazolium
    {'name': 'Pyridinium', 'smiles': 'c1ccncc1', 'fragment_type': 'cation'},  # Base pyridinium
    {'name': 'Ammonium', 'smiles': '[NH4+]', 'fragment_type': 'cation'},  # Ammonium
    {'name': 'Phosphonium', 'smiles': '[P+]', 'fragment_type': 'cation'},  # Phosphonium

    # Anions (Common anions in research and industry)
    {'name': 'Tetrafluoroborate', 'smiles': 'F[B-](F)(F)F', 'fragment_type': 'anion'},  # [BF4]
    {'name': 'Bis(trifluoromethanesulfonyl)imide', 'smiles': 'C(F)(F)(F)S(=O)(=O)[N-]S(=O)(=O)C(F)(F)F', 'fragment_type': 'anion'},  # [NTf2]
    {'name': 'Chloride', 'smiles': '[Cl-]', 'fragment_type': 'anion'},  # [Cl]
    {'name': 'Bromide', 'smiles': '[Br-]', 'fragment_type': 'anion'},  # [Br]
    {'name': 'Trifluoromethanesulfonate', 'smiles': 'C(F)(F)(F)S(=O)(=O)[O-]', 'fragment_type': 'anion'},  # [OTf]
    {'name': 'Dicyanamide', 'smiles': '[N-](C#N)C#N', 'fragment_type': 'anion'},  # [DCA]

    # Alkyl Chains (Common chain lengths)
    {'name': 'Methyl', 'smiles': 'C', 'fragment_type': 'alkyl_chain'},  # C1
    {'name': 'Ethyl', 'smiles': 'CC', 'fragment_type': 'alkyl_chain'},  # C2
    {'name': 'Propyl', 'smiles': 'CCC', 'fragment_type': 'alkyl_chain'},  # C3
    {'name': 'Butyl', 'smiles': 'CCCC', 'fragment_type': 'alkyl_chain'},  # C4
]