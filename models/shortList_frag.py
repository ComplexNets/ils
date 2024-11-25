fragments = [
    # Cations (Common imidazolium, pyridinium, and ammonium cations)
    {'name': 'Ethylimidazolium', 'smiles': 'CCn1cc[n+]c1', 'fragment_type': 'cation'},  # [EMIM]
    {'name': 'Methylimidazolium', 'smiles': 'Cn1cc[n+]c1', 'fragment_type': 'cation'},  # [MIM]
    {'name': 'Propylimidazolium', 'smiles': 'CCCn1cc[n+]c1', 'fragment_type': 'cation'},  # [PMIM]
    {'name': '1-Ethylpyridinium', 'smiles': 'CC[n+]1ccccc1', 'fragment_type': 'cation'},  # [EPy]

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