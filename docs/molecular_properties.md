# Molecular Properties Description

This document describes the molecular properties used in our ionic liquid fragment calculations.

## Property Definitions

### LogP (Partition coefficient)
- **Description**: Logarithm of the partition coefficient, which measures the molecule's lipophilicity (its ability to dissolve in fats, oils, and non-polar solvents)
- **Significance**: Indicates how well a molecule can permeate cell membranes and its distribution in biological systems
- **Range**: Typically between -3 and 7
- **Interpretation**: Higher values indicate greater lipophilicity

### Hydrogen Bond Properties
#### Hydrogen Bond Donor Count
- **Description**: Number of hydrogen atoms bonded to electronegative atoms (typically O, N, F)
- **Significance**: Important for predicting molecular interactions and solubility
- **Range**: Typically 0-5 for drug-like molecules
- **Interpretation**: Higher counts indicate stronger hydrogen bonding capability

#### Hydrogen Bond Acceptor Count
- **Description**: Number of electronegative atoms (typically O, N, F) with a lone pair of electrons
- **Significance**: Crucial for molecular recognition and solubility
- **Range**: Typically 0-10 for drug-like molecules
- **Interpretation**: Higher counts indicate stronger ability to accept hydrogen bonds

### TPSA (Topological Polar Surface Area)
- **Description**: Sum of surface areas of all polar atoms (primarily oxygen and nitrogen) in the molecule
- **Significance**: Predicts drug absorption, including intestinal absorption and blood-brain barrier penetration
- **Unit**: Square Angstroms (Å²)
- **Range**: Typically 0-140 Å² for drug-like molecules
- **Interpretation**: Higher values indicate greater polarity and typically lower membrane permeability

### Rotatable Bond Count
- **Description**: Number of single bonds that can rotate freely in the molecule
- **Significance**: Measures molecular flexibility, important for binding properties
- **Range**: Typically 0-15 for drug-like molecules
- **Interpretation**: Higher counts indicate greater molecular flexibility

### Complexity
- **Description**: A measure of the molecule's structural complexity based on bond types, symmetry, and atom counts
- **Significance**: Useful for predicting synthesis difficulty and stability
- **Range**: Typically 0-1000+ (arbitrary units)
- **Interpretation**: Higher values indicate more complex molecular structures

### Charge
- **Description**: The overall electrical charge of the molecule or fragment
- **Significance**: Critical for ionic liquid properties and interactions
- **Range**: Typically -1, 0, or +1 for ionic liquid components
- **Interpretation**: Determines ionic behavior and interactions

### Heavy Atom Count
- **Description**: Number of non-hydrogen atoms in the molecule
- **Significance**: Basic measure of molecular size and complexity
- **Range**: Typically 3-50 for ionic liquid fragments
- **Interpretation**: Higher counts indicate larger molecules

## Usage in Ionic Liquid Design
These properties are used to:
1. Screen potential ionic liquid combinations
2. Predict physicochemical properties
3. Assess molecular interactions
4. Evaluate potential applications
5. Guide the design of new ionic liquids with desired properties
