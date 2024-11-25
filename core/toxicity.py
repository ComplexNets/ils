import sys
import os
import math
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from typing import Dict, List, Optional
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
from models.shortList_frag import fragments
from utils.utils import get_fragment_properties
from utils.rdkit_utils import get_rdkit_properties

def get_fragment_smiles(fragment_name: str, fragment_type: str) -> str:
    """Get SMILES string for a given fragment name and type."""
    for fragment in fragments:
        if fragment['name'] == fragment_name and fragment['fragment_type'] == fragment_type:
            return fragment.get('smiles', '')
    return ''

def calculate_toxicity_score(mol_properties: Dict) -> float:
    """
    Calculate a normalized toxicity score based on molecular properties.
    Higher score indicates higher predicted toxicity.
    
    Args:
        mol_properties: Dictionary of molecular properties
    Returns:
        float: Normalized toxicity score (0-1 range)
    """
    # Key molecular descriptors that influence toxicity
    mw = mol_properties.get('molecular_weight', 0)
    logp = mol_properties.get('logp', 0)
    tpsa = mol_properties.get('topological_polar_surface_area', 0)
    rotatable_bonds = mol_properties.get('rotatable_bond_count', 0)
    hbd = mol_properties.get('hydrogen_bond_donor_count', 0)
    hba = mol_properties.get('hydrogen_bond_acceptor_count', 0)
    
    # Normalize each property to 0-1 range based on typical IL ranges
    norm_mw = min(1.0, mw / 500.0)  # Most ILs < 500 g/mol
    norm_logp = (logp + 2) / 8  # Typical range: -2 to 6
    norm_tpsa = min(1.0, tpsa / 150.0)  # Most ILs < 150 Å²
    norm_rb = min(1.0, rotatable_bonds / 10.0)  # Most ILs < 10 rotatable bonds
    norm_hb = min(1.0, (hbd + hba) / 10.0)  # Most ILs < 10 H-bond sites
    
    # Weight the contributions (these weights can be adjusted based on data)
    toxicity_score = (
        0.3 * norm_mw +      # Larger molecules tend to be more toxic
        0.3 * norm_logp +    # Higher logP often correlates with higher toxicity
        0.2 * norm_tpsa +    # TPSA correlates with membrane permeability
        0.1 * norm_rb +      # Flexibility can affect bioavailability
        0.1 * norm_hb        # H-bonding affects interaction with biomolecules
    )
    
    return min(1.0, max(0.0, toxicity_score))

def normalize_to_ic50(toxicity_score: float) -> float:
    """Convert normalized toxicity score to IC50 value in mM."""
    # Convert to IC50 using exponential scale (0.1 mM to 100 mM)
    min_ic50 = 0.1  # Very toxic
    max_ic50 = 100.0  # Less toxic
    
    # Invert score since higher IC50 means less toxic
    inverted_score = 1 - toxicity_score
    
    # Calculate IC50 using exponential scale
    log_min = math.log10(min_ic50)
    log_max = math.log10(max_ic50)
    log_ic50 = log_min + (log_max - log_min) * inverted_score
    
    return 10 ** log_ic50

def calculate_ionic_liquid_toxicity(combination: Dict) -> Optional[Dict]:
    """
    Calculate the toxicity of an ionic liquid using molecular descriptors.
    Returns IC50 value in mM, where higher values indicate lower toxicity.
    
    Args:
        combination: Dictionary containing cation, anion, and alkyl_chain components
    Returns:
        Optional[Dict]: Dictionary containing IC50 value and component properties
    """
    try:
        components = {}
        total_properties = {
            'molecular_weight': 0.0,
            'logp': 0.0,
            'topological_polar_surface_area': 0.0,
            'rotatable_bond_count': 0,
            'hydrogen_bond_donor_count': 0,
            'hydrogen_bond_acceptor_count': 0
        }
        
        # Calculate properties for each component
        for component_type in ['cation', 'anion', 'alkyl_chain']:
            component = combination[component_type]
            smiles = get_fragment_smiles(component['name'], component_type)
            
            if not smiles:
                print(f"No SMILES found for {component['name']}")
                continue
            
            # Try to get properties from database first
            props = get_fragment_properties(component['name'], component_type)
            
            # If not in database, calculate using RDKit
            if not props:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    props = {
                        'molecular_weight': Descriptors.ExactMolWt(mol),
                        'logp': Descriptors.MolLogP(mol),
                        'topological_polar_surface_area': Descriptors.TPSA(mol),
                        'rotatable_bond_count': Descriptors.NumRotatableBonds(mol),
                        'hydrogen_bond_donor_count': Descriptors.NumHDonors(mol),
                        'hydrogen_bond_acceptor_count': Descriptors.NumHAcceptors(mol)
                    }
            
            if props:
                components[component_type] = props
                for key in total_properties:
                    if key in props:
                        total_properties[key] += props[key]
        
        if not components:
            print("No valid components found for toxicity calculation")
            return None
        
        # Calculate overall toxicity score
        toxicity_score = calculate_toxicity_score(total_properties)
        
        # Convert to IC50
        ic50 = normalize_to_ic50(toxicity_score)
        
        return {
            'ic50_mm': ic50,
            'components': components,
            'total_properties': total_properties
        }
        
    except Exception as e:
        print(f"Error calculating toxicity: {str(e)}")
        return None
