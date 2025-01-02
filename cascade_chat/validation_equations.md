# IL Structural Constraints from Paper

## Overview
The designed ILs need to satisfy certain rules to ensure chemical feasibility. These rules, termed as structural constraints, include feasibility rules such as octet rule, bonding rule, and complexity rules.

## Base Variable Definitions
- ci: vector of binary variables representing the cations
- aj: vector of binary variables representing the anions
- yi: vector of binary variables representing the alkyl chains l
- ngil: vector of integer variables representing the number of groups of type i in the alkyl chain l
- vci, vail: vectors of group valences of the cations and alkyl groups, respectively
- G: set of alkyl groups available for the cation side chains

## Equations 12-24

### Basic Structure Constraints (Eq. 12-13)
```
Σci = 1  (12)
i∈C
```
Variables:
- ci: binary variable (0 or 1) indicating if cation i is selected
- C: set of all possible cations
Purpose: Ensures exactly one cation is selected

```
Σaj = 1  (13)
j∈a
```
Variables:
- aj: binary variable (0 or 1) indicating if anion j is selected
- a: set of all possible anions
Purpose: Ensures exactly one anion is selected

### Valence and Bonding (Eq. 14-16)
```
Σyi = Σciνci  (14)
i=1   i∈C
```
Variables:
- yi: binary variable for alkyl chain i
- ci: binary variable for cation i
- vci: valence of cation i
Purpose: Fixes number of alkyl side chains based on cation valence

```

Explanation: 
Σ(2-νci)ci + ΣΣ(2-νail)yingil = 2  (15)
i∈C      l=1 i∈G
```
Variables:
- vci: valence of cation i
- ci: binary variable for cation i
- vail: valence of alkyl group i in chain l
- yingil: number of groups of type i in alkyl chain l
- G: set of alkyl groups
Purpose: Ensures valence balance in the molecule

```
Σyingil(2-νail) = 1  (16)
k∈G
```
Variables:
- yingil: number of groups of type i in alkyl chain l
- vail: valence of alkyl group i
- G: set of alkyl groups
Purpose: Ensures proper bonding in alkyl chains

### Size and Group Constraints (Eq. 17-18)
```
ΣΣyingil ≤ nGl  (17)
l=1 k∈G
```
Variables:
- yingil: number of groups of type i in alkyl chain l
- nGl: maximum number of groups allowed in chain l
- G: set of alkyl groups
Purpose: Controls total number of groups in each chain

```
Σyingil ≤ nGal  (18)
k∈G
```
Variables:
- yingil: number of groups of type i in alkyl chain l
- nGal: maximum number of groups of type i allowed
- G: set of alkyl groups
Purpose: Controls number of specific group types in chains

### Group Occurrence Constraints (Eq. 19-21)
```
Σyingil ≤ t1  (19)
k∈G
```
Variables:
- yingil: number of groups of type i in alkyl chain l
- t1: maximum occurrence threshold
- G: set of alkyl groups
Purpose: Sets upper limit on group occurrences in a chain

```
Σyingil ≥ t2  (20)
k∈G
```
Variables:
- yingil: number of groups of type i in alkyl chain l
- t2: minimum occurrence threshold
- G: set of alkyl groups
Purpose: Sets lower limit on group occurrences in a chain

```
Σyingil = t3  (21)
k∈G
```
Variables:
- yingil: number of groups of type i in alkyl chain l
- t3: exact occurrence requirement
- G: set of alkyl groups
Purpose: Requires exact number of group occurrences in a chain

### Cation-wide Group Constraints (Eq. 22-24)
```
ΣΣyingil ≤ t4  (22)
l=1 k∈G*
```
Variables:
- yingil: number of groups of type i in alkyl chain l
- t4: maximum total occurrence threshold
- G*: specific subset of alkyl groups
Purpose: Sets upper limit on total group occurrences across all chains

```
ΣΣyingil ≥ t5  (23)
l=1 k∈G*
```
Variables:
- yingil: number of groups of type i in alkyl chain l
- t5: minimum total occurrence threshold
- G*: specific subset of alkyl groups
Purpose: Sets lower limit on total group occurrences across all chains

```
ΣΣyingil = t6  (24)
l=1 k∈G*
```
Variables:
- yingil: number of groups of type i in alkyl chain l
- t6: exact total occurrence requirement
- G*: specific subset of alkyl groups
Purpose: Requires exact number of total group occurrences across all chains

## Notes
- Modified octet rule ensures structural feasibility and valence satisfaction
- Cation size controlled by upper bound on maximum number of groups
- Alkyl chain size controlled by upper bound on maximum groups per chain
- Equations 19-21 place restrictions on occurrences in side chains
- Equations 22-24 place restrictions on occurrences in whole cation
