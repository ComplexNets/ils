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

def calculate_toxicity_score(mol_properties: Dict, component_types: Dict) -> float:
    """
    Calculate a normalized toxicity score based on molecular properties.
    Higher score indicates higher predicted toxicity.
    
    Args:
        mol_properties: Dictionary of molecular properties
        component_types: Dictionary of component types (cation, anion, etc.)
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
    
    # Get cation and anion types for specific adjustments
    cation_name = component_types.get('cation', '').lower()
    anion_name = component_types.get('anion', '').lower()
    alkyl_name = component_types.get('alkyl_chain', '').lower()
    functional_group = component_types.get('functional_group', '').lower()
    
    # Adjust normalization parameters based on component types
    # Different cation families have different toxicity profiles
    mw_factor = 500.0  # Default
    if 'imidazolium' in cation_name:
        mw_factor = 450.0
    elif 'pyridinium' in cation_name:
        mw_factor = 475.0
    elif 'ammonium' in cation_name:
        mw_factor = 525.0
    elif 'phosphonium' in cation_name:
        mw_factor = 550.0
        
    # LogP normalization adjustments
    logp_min = -2.0  # Default minimum
    logp_range = 8.0  # Default range
    
    # Adjust based on anion type
    if 'bis(trifluoromethanesulfonyl)imide' in anion_name:
        logp_min = -1.0
        logp_range = 9.0
    elif 'tetrafluoroborate' in anion_name:
        logp_min = -2.5
        logp_range = 7.5
        
    # Adjust based on alkyl chain length
    alkyl_length_factor = 1.0
    if 'butyl' in alkyl_name:
        alkyl_length_factor = 1.2
    elif 'hexyl' in alkyl_name:
        alkyl_length_factor = 1.4
    elif 'octyl' in alkyl_name:
        alkyl_length_factor = 1.6
    elif 'decyl' in alkyl_name:
        alkyl_length_factor = 1.8
        
    # Functional group adjustments
    functional_group_factor = 1.0
    if functional_group:
        if 'hydroxyl' in functional_group:
            functional_group_factor = 0.8  # Hydroxyl groups generally reduce toxicity
        elif 'amino' in functional_group:
            functional_group_factor = 0.9
        elif 'carboxyl' in functional_group:
            functional_group_factor = 0.75
            
    # Normalize each property to 0-1 range based on adjusted parameters
    norm_mw = min(1.0, mw / mw_factor)
    norm_logp = (logp - logp_min) / logp_range
    norm_tpsa = min(1.0, tpsa / 150.0)
    norm_rb = min(1.0, rotatable_bonds / 10.0)
    norm_hb = min(1.0, (hbd + hba) / 10.0)
    
    # Ensure normalized values are in 0-1 range
    norm_mw = min(1.0, max(0.0, norm_mw))
    norm_logp = min(1.0, max(0.0, norm_logp))
    
    # Weight the contributions (these weights can be adjusted based on data)
    # Apply the functional group and alkyl chain length factors
    toxicity_score = (
        0.3 * norm_mw * alkyl_length_factor +
        0.3 * norm_logp * alkyl_length_factor +
        0.2 * norm_tpsa +
        0.1 * norm_rb +
        0.1 * norm_hb
    ) * functional_group_factor
    
    return min(1.0, max(0.0, toxicity_score))

def normalize_to_ic50(toxicity_score: float, component_types: Dict) -> float:
    """
    Convert normalized toxicity score to IC50 value in mM.
    Adjusts the scale based on component types.
    """
    # Get component types for specific adjustments
    cation_name = component_types.get('cation', '').lower()
    anion_name = component_types.get('anion', '').lower()
    
    # Adjust IC50 range based on component types
    min_ic50 = 0.1  # Default very toxic
    max_ic50 = 100.0  # Default less toxic
    
    # Different cation families have different toxicity ranges
    if 'imidazolium' in cation_name:
        min_ic50 = 0.08
        max_ic50 = 90.0
    elif 'pyridinium' in cation_name:
        min_ic50 = 0.05
        max_ic50 = 80.0
    elif 'ammonium' in cation_name:
        min_ic50 = 0.15
        max_ic50 = 120.0
    elif 'phosphonium' in cation_name:
        min_ic50 = 0.12
        max_ic50 = 110.0
        
    # Adjust based on anion type
    if 'bis(trifluoromethanesulfonyl)imide' in anion_name:
        min_ic50 *= 0.9  # More toxic
    elif 'tetrafluoroborate' in anion_name:
        max_ic50 *= 1.1  # Less toxic
        
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
        
        # Store component types for toxicity adjustments
        component_types = {
            'cation': combination['cation']['name'],
            'anion': combination['anion']['name'],
            'alkyl_chain': combination['alkyl_chain']['name']
        }
        
        # Add functional group if present
        if 'functional_group' in combination and combination['functional_group']:
            component_types['functional_group'] = combination['functional_group']['name']
        
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
        
        # Add functional group properties if present
        if 'functional_group' in combination and combination['functional_group']:
            func_group = combination['functional_group']
            smiles = get_fragment_smiles(func_group['name'], 'functional_group')
            
            if smiles:
                props = get_fragment_properties(func_group['name'], 'functional_group')
                
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
                    components['functional_group'] = props
                    for key in total_properties:
                        if key in props:
                            total_properties[key] += props[key]
        
        if not components:
            print("No valid components found for toxicity calculation")
            return None
        
        # Calculate overall toxicity score with component-specific adjustments
        toxicity_score = calculate_toxicity_score(total_properties, component_types)
        
        # Convert to IC50 with component-specific adjustments
        ic50 = normalize_to_ic50(toxicity_score, component_types)
        
        return {
            'ic50_mm': ic50,
            'components': components,
            'total_properties': total_properties,
            'toxicity_score': toxicity_score
        }
        
    except Exception as e:
        print(f"Error calculating toxicity: {str(e)}")
        return None
