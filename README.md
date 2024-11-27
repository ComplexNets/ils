# Ionic Liquid Selector (ILS)

A tool for designing and optimizing ionic liquids based on desired properties.

## Setup

1. Install Python 3.8 or higher
2. Install dependencies:
```bash
pip install -r requirements.txt
```

## Running the Application

To run the application:
```bash
streamlit run frontend/ils_ui.py
```

## Data Sources

Fragment data is stored in `fragment_data/autono17_ilselect_db.csv`. This file contains:
- Molecular properties for cations, anions, and alkyl chains
- SMILES representations
- Physical properties like molecular weight, density, etc.

## Features

- Interactive property range selection
- Multi-objective optimization using Pareto front analysis
- Property visualization using parallel coordinates and radar plots
- Fragment screening based on desired properties
- Heat capacity, density, and toxicity predictions
