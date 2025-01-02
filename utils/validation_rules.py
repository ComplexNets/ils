from typing import Dict, Tuple, List
import numpy as np
from itertools import product
import os
from utils.utils import get_fragment_properties
from rdkit import Chem

# Constants for chemical rules
MAX_BOND_CAPACITY = 4  # Maximum total bond capacity

class MolecularValidator:
    """Validates ionic liquid combinations based on chemical rules"""
    
    def __init__(self):
        """Initialize validator"""
        pass

# Cation and Anion Base Count
# Equations 12 and 13
# Purpose: Ensures a maximum of one cation base and one anion, respectively, for each IL candidate.

    def _get_anion_cation_count(self, fragments: List[Dict], fragment_type: str) -> int:
        """Get the count of fragments of a specific type"""
        return sum(1 for fragment in fragments if fragment.get('fragment_type', '').lower() == fragment_type)
    
    def _validate_anion_cation_counts(self, fragments: List[Dict]) -> Tuple[bool, str]:
        """
        Validate that there is exactly one cation and one anion
        Returns:
            Tuple[bool, str]: (is_valid, message)
        """
        cation_count = self._get_anion_cation_count(fragments, 'cation')
        anion_count = self._get_anion_cation_count(fragments, 'anion')
        
        if cation_count != 1:
            return False, f"Invalid number of cations: {cation_count} (must be 1)"
        if anion_count != 1:
            return False, f"Invalid number of anions: {anion_count} (must be 1)"

        return True, "Valid fragment counts"
# Alkyl Chain Count
# Equation 14
# Purpose: Fixes number of alkyl side chains based on cation valence

    def _get_valence_from_smiles(self, smiles: str) -> List[Dict]:
        """Get valence information for each atom in a SMILES string.
        
        Args:
            smiles: SMILES string to analyze
            
        Returns:
            List of dictionaries containing valence info for each atom
        """
        try:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                print(f"DEBUG: Failed to parse SMILES: {smiles}")
                return None
            
            valences = []
            for atom in mol.GetAtoms():
                # Calculate total bonds and max valence
                total_bonds = sum([bond.GetBondTypeAsDouble() for bond in atom.GetBonds()])
                max_valence = 4 if atom.GetSymbol() in ["N", "P"] else 3
                
                valence_info = {
                    "Atom": atom.GetSymbol(),
                    "Total Bonds": total_bonds,
                    "Max Valence": max_valence,
                    "Formal Charge": atom.GetFormalCharge(),
                    "Is Aromatic": atom.GetIsAromatic()
                }
                print(f"DEBUG: Atom valence info: {valence_info}")
                valences.append(valence_info)
            return valences
        except Exception as e:
            print(f"DEBUG: Error processing SMILES {smiles}: {str(e)}")
            return None

    def _alkyl_chain_count(self, fragments: List[Dict]) -> Tuple[bool, str]:
        """Eq. 14 fixes the number of alkyl side chains attached to the cation based on 
        the available free valence of the cation base.
        Returns:
            Tuple[bool, str]: (is_valid, message)
        """
        # Find the cation fragment
        cation = next((f for f in fragments if f['fragment_type'].lower() == 'cation'), None)
        if not cation:
            print("DEBUG: No cation found in fragments")
            return False, "No cation found"
        
        print(f"DEBUG: Processing cation: {cation['name']}")
        
        # Get cation SMILES and convert to RDKit mol
        cation_smiles = cation.get('smiles')
        if not cation_smiles:
            print("DEBUG: No SMILES found for cation")
            return False, "No SMILES found for cation"
            
        print(f"DEBUG: Cation SMILES: {cation_smiles}")
            
        # Get valence information
        valences = self._get_valence_from_smiles(cation_smiles)
        if not valences:
            print(f"DEBUG: Could not analyze valence for cation: {cation_smiles}")
            return False, f"Could not analyze valence for cation: {cation_smiles}"
            
        # Calculate available attachment points based on valence
        available_points = 0
        for v in valences:
            print(f"DEBUG: Processing atom: {v}")
            atom_symbol = v["Atom"]
            total_bonds = v["Total Bonds"]
            max_valence = v["Max Valence"]
            formal_charge = v["Formal Charge"]
            is_aromatic = v["Is Aromatic"]
            
            # For ammonium/phosphonium centers (NH4+ or PH4+)
            if atom_symbol in ["N", "P"] and formal_charge == 1 and not is_aromatic:
                available_points = 4  # Can replace all hydrogens
                print(f"DEBUG: Found tetravalent {atom_symbol}+ center with 4 available points")
                break
            
            # For nitrogen in aromatic rings (imidazolium, pyridinium)
            elif atom_symbol == "N" and is_aromatic:
                # For charged nitrogen (N+)
                if formal_charge == 1:
                    # If not fully substituted
                    if total_bonds < 3:
                        available_points += 1
                        print(f"DEBUG: Found available point on aromatic N+")
                # For uncharged NH in imidazolium
                elif total_bonds < max_valence:
                    available_points += 1
                    print(f"DEBUG: Found available point on aromatic NH")
        
        print(f"DEBUG: Total available points: {available_points}")
        
        # Count alkyl chains in the combination
        alkyl_chain_count = sum(1 for f in fragments if f['fragment_type'].lower() == 'alkyl_chain')
        print(f"DEBUG: Alkyl chain count: {alkyl_chain_count}")
        
        # Validate that we don't exceed available points
        if alkyl_chain_count > available_points:
            print(f"DEBUG: Too many alkyl chains ({alkyl_chain_count}) for available points ({available_points})")
            return False, f"Too many alkyl chains ({alkyl_chain_count}) for available attachment points ({available_points})"
        
        print(f"DEBUG: Valid alkyl chain count")
        return True, f"Valid number of alkyl chains ({alkyl_chain_count}) for available attachment points ({available_points})"
   
    # OCTET RULE
    # Equation 15
    # Ensures valence balance in the molecule




    def _get_fragment_bond_capacity(self, fragment: Dict) -> int:
        """Get fragment bond capacity based on properties"""
        props = get_fragment_properties(fragment['name'], fragment['fragment_type'])
        if not props:
            return 0
        return min(props.get('rotatable_bond_count', 0) + 1, MAX_BOND_CAPACITY)

    def validate(self, cation: Dict, anion: Dict, alkyl: Dict) -> Tuple[bool, str]:
        """
        Validate a combination of fragments
        Returns:
            Tuple[bool, str]: (is_valid, message)
        """
        try:
            # Check fragment counts
            fragments = [cation, anion, alkyl]
            is_valid, message = self._validate_anion_cation_counts(fragments)
            if not is_valid:
                return False, message

            # Check if we can get properties for all fragments
            fragments = [
                (cation, "Cation"),
                (anion, "Anion"),
                (alkyl, "Alkyl chain")
            ]
            
            for frag, frag_type in fragments:
                props = get_fragment_properties(frag['name'], frag['fragment_type'])
                if not props:
                    return False, f"{frag_type} {frag.get('name')} properties not found"
            
            # Validate alkyl chain count based on cation valence
            fragments = [cation, anion, alkyl]
            is_valid, message = self._alkyl_chain_count(fragments)
            if not is_valid:
                return False, message
            
            # Get and check bond capacities
            cat_bonds = self._get_fragment_bond_capacity(cation)
            an_bonds = self._get_fragment_bond_capacity(anion)
            alkyl_bonds = self._get_fragment_bond_capacity(alkyl)
            
            # Check total bonds
            total_bonds = cat_bonds + an_bonds + alkyl_bonds
            if total_bonds > MAX_BOND_CAPACITY:
                return False, f"Too many bonds: {total_bonds}"
            
            return True, "Valid combination"
            
        except Exception as e:
            return False, f"Validation error: {str(e)}"
