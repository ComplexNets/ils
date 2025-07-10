# Ionic Liquid Screening (ILS) Project Summary

MINICONDA ENVIRONMENT: C:\Users\X1\miniconda3\envs\ils

## Overall Project Progress

### Core Features Implemented

1. **Fragment Management**

   - Short and long fragment lists defined in `models/shortList_frag.py`
   - Fragment selection system in `utils/fragment_selector.py`
   - Property calculation and storage for each fragment
2. **Validation System**

   - Robust validation rules in `utils/validation_rules.py`
   - Handles cation, anion, and alkyl chain combinations
   - Validates based on chemical rules and molecular structure
   - Uses RDKit for molecular property calculations
3. **Fragment Combination**

   - Main combination logic in `core/combine_fragments.py`
   - Parallel processing for efficient combination generation
   - Filters combinations based on property ranges
   - Integration with IL Thermo database
4. **Naming and Properties**

   - Standardized IL naming system in `utils/utils.py`
   - Property calculation using RDKit and database lookups
   - Integration with external databases

### Project Structure

```
ils/
├── core/
│   ├── combine_fragments.py  # Main fragment combination logic
│   ├── density.py            # Density calculation
│   ├── heat_capacity.py      # Heat capacity calculation
│   ├── hydrophobicity.py     # Hydrophobicity calculation
│   ├── pareto_optimizer.py   # Pareto optimization logic
│   ├── solubility.py         # Solubility calculation
│   ├── toxicity.py           # Toxicity prediction
│   └── viscosity.py          # Viscosity calculation
├── models/
│   └── shortList_frag.py     # Fragment definitions
├── utils/
│   ├── fragment_selector.py  # Fragment list management
│   ├── validation_rules.py   # Chemical validation rules
│   ├── utils.py             # General utilities
│   └── rdkit_utils.py       # RDKit integration
```

## Technical Details

### Validation Process

1. Fragment Count Validation

   - Ensures exactly one cation and one anion
   - Validates number of alkyl chains
2. Chemical Structure Validation

   - Uses RDKit for molecular analysis
   - Validates valence and bonding
   - Handles special cases (aromatic systems, charged centers)
3. Property Validation

   - Checks molecular properties
   - Validates bond capacities
   - Ensures chemical feasibility

### Performance Considerations

- Parallel processing using ThreadPoolExecutor
- Efficient fragment filtering
- Caching of molecular properties
- Integration with external databases for validation

## Primary File Explanations

This section provides a brief overview of the core scripts identified in the Project Structure.

### Core Modules (`core/`)

* **`combine_fragments.py`:** The main engine for generating potential ionic liquid (IL) candidates. It loads predefined fragment lists (cations, anions, alkyl chains, functional groups), fetches properties, generates all combinations, and uses `MolecularValidator` to check chemical validity in parallel. Calculates molecular weight and checks against the ILThermo database for valid combinations.
* **`density.py`:** Calculates IL density (kg/m³) using a modified Ye & Shreeve Group Contribution Method with temperature correction. Estimates molecular volume based on atomic volumes, ionic radii, ring corrections, and environmental factors. Combines fragment weights and volumes for the final calculation at standard conditions.
* **`heat_capacity.py`:** Estimates IL heat capacity (Cp in J/mol·K) using a mix of group contribution methods (UNIFAC-like parameters) and empirical correlations based on molecular weight, atom counts, bond flexibility, and hydrogen bonding. Sums fragment contributions and adds interaction terms.
* **`hydrophobicity.py`:** Calculates IL hydrophobicity (logP) using RDKit's Crippen method. Sums the logP contributions of individual fragments (cation, anion, alkyl, functional group) based on their SMILES strings and applies a temperature correction.
* **`pareto_optimizer.py`:** Implements multi-objective optimization using Pareto dominance. Allows defining property constraints (min/max, weight, direction). Identifies the Pareto front (non-dominated solutions) from a set of candidates and can rank all solutions based on a weighted score.
* **`solubility.py`:** Estimates water solubility (g/L) using a simplified Group Contribution Method. Calculates contributions for each fragment based on name/structure, sums them, applies a temperature correction, and validates the result.
* **`toxicity.py`:** Estimates IL toxicity (IC50 in mM, higher value = lower toxicity). Calculates RDKit molecular descriptors for the combined IL, normalizes them into a weighted toxicity score (0-1) adjusted by component types, and converts this score to an estimated IC50 value using an exponential scale.
* **`viscosity.py`:** Estimates IL viscosity (Pa·s) using a Quantitative Structure-Property Relationship (QSPR) model. Calculates RDKit descriptors for the combined IL, feeds them into a linear regression model to predict log(viscosity), applies temperature correction (Arrhenius equation), and uses an empirical calibration factor.

### Model Definitions (`models/`)

* **`shortList_frag.py`:** Defines a small, curated list of common molecular fragments (cations, anions, alkyl chains, functional groups) with their names, SMILES strings, and types. Used as input for generating a limited set of IL candidates. (Note: `mediumList_frag.py` and `longList_frag.py` likely exist for larger fragment sets).

### Utility Modules (`utils/`)

* **`fragment_selector.py`:** Utility class (`FragmentSelector`) to manage and select different predefined fragment lists (e.g., 'short', 'medium', 'long') imported from the `models` directory.
* **`validation_rules.py`:** Defines the `MolecularValidator` class for checking the chemical validity of IL combinations. Enforces rules like cation/anion counts, alkyl chain counts based on cation valence (using RDKit), and potentially constraints on functional groups per chain/molecule.

* **`utils.py`:** Provides common utility functions: standardizing and generating IUPAC-compliant IL names, loading/caching fragment properties from a CSV database (`fragment_data/autono17_ilselect_db.csv`), retrieving properties (falling back to RDKit), and checking IL existence in the ILThermo database (via local list and API query).
* **`rdkit_utils.py`:** Centralizes RDKit calculations. The `get_rdkit_properties` function takes a SMILES string and returns a dictionary of standard molecular properties and descriptors (MW, TPSA, LogP, H-bond counts, etc.).

### Description of ILSCOPE In current paper:

To evaluate the contribution of fragments to IL properties and
identify potential fine-tuning opportunities, IL-SCOPE was developed. This
Python-based software utilizes a reverse-forward design strategy in screening
ILs and calculating their properties to find optimal IL candidates based on
user-defined property ranges and mission objective weighting. The reverse
design phase incorporates the user defined ranges and works backwards to screen
fragments(cation, anion, alkyl chains and functional groups) that will likely
fall outside of these ranges.  This
approach avoids the need to presuppose which fragments are viable, allowing for
broader design space exploration without biasing the outcome toward
preconceived fragment choices.The forward phase then assembles these fragments
into IL candidates applies structural (AND OTHER) rules and calculates the
properties of the remaining feasible ILs. The feasible IL properties are then
weighted according to user input priorities yielding a subset of optimal ILs

## A. Reverse Design

![](file:///C:/Users/X1/AppData/Local/Temp/msohtmlclip1/01/clip_image002.gif)The reverse design
strategy works backward from target property values to loosely screen out
non-viable fragments, improving computational efficiency. IL-SCOPE uses
RDKit-derived descriptors (e.g., molecular weight for density) and UNIFAC-like
parameters (e.g., for heat capacity) to estimate each fragment’s contribution
using a scaled portion, plus margin, of the user-defined property range. If the
estimated value falls outside this range, the fragment is discarded. This
approach does not require the researcher to presuppose which fragments are
viable, but instead enables the import of a full complement of fragments for
initial, coarse-grained screening that reduces the solution set prior to IL
assembly.

## B. Forward Design

The forward design phase combines screened fragments into
candidate ILs, which are then evaluated for structural feasibility. These
candidates must satisfy a set of 13 constraint rules, based on the Octet Rule,
Bonding Logic, and Complexity Constraints (K_) that govern fragment
compatibility and molecular assembly. Only ILs that pass this validation step
are retained.

For each validated IL, properties are then calculated using the
models listed in Table X. While these models share a similar structure to those
used in the reverse phase, they operate on fully assembled ILs, accounting for
both individual fragment contributions and interactions between fragments
(e.g., hydrogen bonding, ion-pair interactions). The resulting properties are
then screened against user-defined target ranges, and any IL falling outside
those bounds is removed from the solution set.

| **Property**       | **Formula Type**      | **Key Methods**                              |
| ------------------------ | --------------------------- | -------------------------------------------------- |
| **Heat Capacity**  | Group Contribution          | UNIFAC, empirical scaling^20^                      |
| **Viscosity**      | QSPR Model                  | Molecular descriptors, Arrhenius correction^20,21^ |
| **Density**        | Molecular Volume Estimation | Ye & Shreeve + Gardas & Coutinho correction^22^    |
| **Hydrophobicity** | RDKit LogP Calculation      | Crippen LogP, temperature correction               |
| **Solubility**     | Group Contribution          | Modified UNIFAC-IL^23^                             |
| **Toxicity**       | QSAR-Based                  | IC50 from molecular descriptors^13,24^             |

## C. Pareto Optimization

The validated ionic liquid is then assigned a weighted score,
based on user-input, that ranges from 0 - 10, in increments of 1. These scores
are calculated and, using a multi-objective (Pareto) Optimization method, each
Il is ranked, ensuring no property can be changed without negatively impacting
other properties. While the Pareto optimization provides a ranked list of IL
candidates, the use of approximate models means that differences in score
should not be interpreted as substantive variability of performance. The
top-performing ILs represent a performance band, rather than a definitive
ranking. This band can be used as a basis for further analysis, including
identifying fragment–property correlations and exploring how specific fragments
contribute to high-performing profiles.

## D. Summary

IL-SCOPE enables effieicny early stage design fuidance and
fragment screening that lows researchers and engineer prform low to medium fidelity
property evalutions, based on mission requirements, of candidate ILs. This provides
users a ‘stepping off’ point for further research into high-fidelity models
like Aspen or Cosmos-rs to ensure property accuracy and Il stability and
feasibility. IL-SCOPE results represent a first step in fragment-level property
correlation analysis, intended to support early-stage investigations of
specific use cases for ILs. The tool provides mission objective, weighted
guidance constrained by user-defined property range requirements, enabling
informed screening and refinement of IL candidates.
