# Ionic Liquid Selector (ILS) Program Summary

## Overview
The Ionic Liquid Selector (ILS) is a sophisticated software tool designed to help researchers identify and optimize ionic liquids for specific applications. It combines fragment-based molecular design with multi-criteria optimization to predict and evaluate key properties of ionic liquids.

## Core Components

### 1. Fragment Combination (`combine_fragments.py`)
- Manages the assembly of ionic liquids from basic fragments (cations, anions, and alkyl chains)
- Validates molecular combinations using chemical rules
- Implements parallel processing for efficient combination generation
- Handles database integration for property lookup

### 2. Property Calculations

#### Density Calculation (`density.py`)
- Implements modified Ye & Shreeve method with temperature correction
- Features:
  - Atomic volume contributions
  - Ring corrections for cyclic structures
  - Environment-specific correction factors
  - Special handling for ionic species
  - Temperature and pressure dependencies

#### Heat Capacity Calculation (`heat_capacity.py`)
- Calculates molecular heat capacities using group contribution methods
- Features:
  - Fragment-specific contributions
  - Hydrogen bonding effects
  - Temperature-dependent corrections
  - UNIFAC method integration
  - Empirical correlation adjustments

#### Toxicity Prediction (`toxicity.py`)
- Predicts ionic liquid toxicity using structure-property relationships
- Features:
  - Molecular descriptor analysis
  - IC50 value estimation
  - Structure-based toxicity scoring
  - Component-wise toxicity evaluation

### 3. Optimization Framework

#### Pareto Optimization (`pareto_optimizer.py`)
- Implements multi-objective optimization for ionic liquid properties
- Features:
  - Non-dominated sorting
  - Property normalization
  - Weighted scoring
  - Constraint handling
  - Pareto front identification

#### Validation Rules (`validation_rules.py`)
- Enforces chemical and physical constraints
- Validates molecular structures
- Checks property ranges
- Ensures fragment compatibility

## User Interface

### Frontend Components
- Interactive property input
- Real-time visualization
- Results display with:
  - Pareto front plots
  - Property correlation analysis
  - Distribution visualization
  - Detailed solution rankings

### Key Features
- Property range specification
- Importance weighting
- Constraint definition
- Interactive visualization
- Export capabilities

## Technical Details

### Dependencies
- Python 3.7+
- Key packages:
  - streamlit (1.24.0)
  - pandas (2.0.3)
  - numpy (1.24.3)
  - plotly (5.15.0)
  - requests (2.31.0)

### Architecture
- Modular design with clear separation of concerns
- Efficient data structures for molecular representation
- Parallel processing capabilities
- Database integration for property lookup
- Extensible framework for adding new properties

## Future Enhancements
1. Additional property predictions
2. Machine learning integration
3. Expanded fragment database
4. Advanced visualization options
5. Batch processing capabilities

## Usage
1. Run the application using `run_ui.bat`
2. Input desired property ranges and constraints
3. Execute optimization
4. Review and analyze results
5. Export selected combinations

## Conclusion
The Ionic Liquid Selector provides a comprehensive platform for designing and optimizing ionic liquids, combining theoretical models with practical constraints to deliver actionable results for researchers and engineers.
