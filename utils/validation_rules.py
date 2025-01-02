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
    
    def __init__(self, validation_criteria=None):
        """Initialize validator with optional validation criteria"""
        self.validation_criteria = validation_criteria

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
    # Equation 15, 16
    # Ensures valence balance in the molecule


    # Upper Bound Cation and Alkyl Chain Size
    # Cation size: The size of the cation is controlled by introducing an upper bound on maximum number of groups nU  G 
    # that are allowed in the cation (Eq. 17).
    # Alkyl chain size: The size of the alkyl chain is controlled by introducing an upper bound on maximum number of groups nU  A 
    # that are allowed in the alkyl chain (Eq. 18).


    def _get_max_groups(self):
        """Get maximum groups allowed per chain"""
        return self.validation_criteria.max_groups_per_chain if self.validation_criteria else 6
        
    def _get_max_groups_per_type(self):
        """Get maximum groups allowed per type"""
        return self.validation_criteria.max_groups_per_type if self.validation_criteria else 2

    def _validate_total_groups_per_chain(self, alkyl_chain: Dict) -> Tuple[bool, str]:
        """Eq. 17: Validates that total number of functional groups on an alkyl chain doesn't exceed maximum allowed
        
        ΣΣyingil ≤ nGl : Sum of all functional groups on chain l must not exceed maximum allowed groups nGl
        
        Note: This counts modifications to the base alkyl chain, not the chain carbons themselves.
        For example, a butyl chain (CCCC) with one -OH group has:
        - Base chain length: 4 (not counted here)
        - Functional groups: 1 (-OH) (counted here)
        
        Returns:
            Tuple[bool, str]: (is_valid, message)
        """
        try:
            # Get alkyl chain properties
            props = get_fragment_properties(alkyl_chain['name'], alkyl_chain['fragment_type'])
            if not props:
                return False, f"Could not get properties for alkyl chain: {alkyl_chain['name']}"
            
            # Get total functional groups (modifications) on chain
            total_functional_groups = props.get('functional_group_count', 0)  # Total number of modifications
            max_allowed_groups = self._get_max_groups()  # Get max allowed functional groups from validation criteria
            
            print(f"DEBUG: Alkyl chain {alkyl_chain['name']} has {total_functional_groups} functional groups (max allowed: {max_allowed_groups})")
            
            if total_functional_groups > max_allowed_groups:
                return False, f"Too many functional groups on alkyl chain {alkyl_chain['name']}: {total_functional_groups} (max: {max_allowed_groups})"
            
            return True, f"Valid number of functional groups on chain ({total_functional_groups})"
            
        except Exception as e:
            return False, f"Error validating total functional groups on chain: {str(e)}"

    def _validate_specific_group_counts(self, alkyl_chain: Dict) -> Tuple[bool, str]:
        """Eq. 18: Validates that number of specific functional group types doesn't exceed their limits
        
        Σyingil ≤ nGal : Sum of functional groups of type i on chain l must not exceed maximum allowed nGal
        
        Note: This counts specific types of modifications to the base alkyl chain.
        For example, a butyl chain (CCCC) with two -OH groups has:
        - Base chain: CCCC (not counted here)
        - Functional groups: 2 (-OH) (counted here)
        
        Returns:
            Tuple[bool, str]: (is_valid, message)
        """
        try:
            # Get alkyl chain properties
            props = get_fragment_properties(alkyl_chain['name'], alkyl_chain['fragment_type'])
            if not props:
                return False, f"Could not get properties for alkyl chain: {alkyl_chain['name']}"
            
            # Get functional group counts
            functional_group_counts = props.get('functional_group_counts', {})  # Dictionary of {group_type: count}
            max_group_counts = self._get_max_groups_per_type()  # Get max allowed groups per type from validation criteria
            
            for group_type, count in functional_group_counts.items():
                print(f"DEBUG: Functional group {group_type} has count {count} (max allowed: {max_group_counts})")
                
                if count > max_group_counts:
                    return False, f"Too many groups of type {group_type} on chain: {count} (max: {max_group_counts})"
            
            return True, "Valid functional group type counts on chain"
            
        except Exception as e:
            return False, f"Error validating specific functional group counts: {str(e)}"

    def _validate_group_occurrences(self, alkyl_chain: Dict) -> Tuple[bool, str]:
        """Group Occurrence Constraints (Equations 19-21)
        
        Eq. 19: Σyingil ≤ t1 (upper limit)
        Eq. 20: Σyingil ≥ t2 (lower limit)
        Eq. 21: Σyingil = t3 (exact count)
        
        Returns:
            Tuple[bool, str]: (is_valid, message)
        """
        try:
            if not self.validation_criteria:
                return True, "No occurrence constraints set"
                
            # Get alkyl chain properties
            props = get_fragment_properties(alkyl_chain['name'], alkyl_chain['fragment_type'])
            if not props:
                return False, f"Could not get properties for alkyl chain: {alkyl_chain['name']}"
            
            # Get group counts
            group_counts = props.get('group_type_counts', {})
            
            for group_type, count in group_counts.items():
                # Check upper limit (t1) - Eq. 19
                if count > self.validation_criteria.group_occurrence_upper:
                    return False, f"Group {group_type} occurs {count} times, exceeding upper limit of {self.validation_criteria.group_occurrence_upper}"
                
                # Check lower limit (t2) - Eq. 20
                if count < self.validation_criteria.group_occurrence_lower:
                    return False, f"Group {group_type} occurs {count} times, below lower limit of {self.validation_criteria.group_occurrence_lower}"
                
                # Check exact count (t3) - Eq. 21
                if self.validation_criteria.group_occurrence_exact >= 0 and count != self.validation_criteria.group_occurrence_exact:
                    return False, f"Group {group_type} occurs {count} times, must occur exactly {self.validation_criteria.group_occurrence_exact} times"
                
                print(f"DEBUG: Group {group_type} occurs {count} times (limits: {self.validation_criteria.group_occurrence_lower}-{self.validation_criteria.group_occurrence_upper})")
            
            return True, "Valid group occurrences"
            
        except Exception as e:
            return False, f"Error validating group occurrences: {str(e)}"

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
            
            # Validate total groups per chain
            is_valid, message = self._validate_total_groups_per_chain(alkyl)
            if not is_valid:
                return False, message
            
            # Validate specific group counts
            is_valid, message = self._validate_specific_group_counts(alkyl)
            if not is_valid:
                return False, message
            
            # Validate group occurrences
            is_valid, message = self._validate_group_occurrences(alkyl)
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
