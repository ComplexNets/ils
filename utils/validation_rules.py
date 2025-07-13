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
            
        # Calculate available attachment points based on valence and implicit hydrogens
        mol = Chem.MolFromSmiles(cation_smiles) # Re-parse to get atom objects easily
        available_points = 0
        processed_central_atom = False
        for atom in mol.GetAtoms():
            atom_symbol = atom.GetSymbol()
            formal_charge = atom.GetFormalCharge()
            is_aromatic = atom.GetIsAromatic()
            implicit_hs = atom.GetNumImplicitHs()
            
            # print(f"DEBUG: Checking Atom {atom.GetIdx()}: Symbol={atom_symbol}, Aromatic={is_aromatic}, Charge={formal_charge}, ImplicitHs={implicit_hs}")

            # For non-aromatic ammonium/phosphonium centers (NH4+ or PH4+)
            if atom_symbol in ["N", "P"] and formal_charge == 1 and not is_aromatic:
                # Assume all 4 positions are potentially substitutable if represented like [NH4+] or [PH4+]
                # A more precise check might look at explicit Hs if the SMILES includes them.
                available_points = 4 
                print(f"DEBUG: Found non-aromatic {atom_symbol}+ center. Assuming 4 available points.")
                processed_central_atom = True
                break # Found the central atom, primary determinant of chain count

            # For aromatic rings (check C and N atoms with implicit Hs)
            elif is_aromatic:
                # Check aromatic carbons with implicit hydrogens
                if atom_symbol == "C" and implicit_hs > 0:
                     available_points += implicit_hs # Each implicit H represents a potential attachment point
                     print(f"DEBUG: Found {implicit_hs} available point(s) on aromatic C (Atom {atom.GetIdx()})")
                # Check aromatic NH groups (like in imidazolium)
                elif atom_symbol == "N" and formal_charge == 0 and implicit_hs > 0:
                     available_points += implicit_hs # The H on the N can be replaced
                     print(f"DEBUG: Found {implicit_hs} available point(s) on aromatic NH (Atom {atom.GetIdx()})")
                # Note: Aromatic N+ usually doesn't have implicit Hs for direct substitution

        # If we didn't find a central N+/P+ atom, the sum from ring atoms is the total.
        print(f"DEBUG: Calculated total available points: {available_points}")
        
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

    def _validate_total_groups_per_chain(self, alkyl_chain: Dict, functional_groups: List[Dict]) -> Tuple[bool, str]:
        """Eq. 17: Validates that total number of functional groups on an alkyl chain doesn't exceed maximum allowed
        
        ΣΣyingil ≤ nGl : Sum of all functional groups on chain l must not exceed maximum allowed groups nGl
        
        Args:
            alkyl_chain: The alkyl chain fragment
            functional_groups: List of functional group fragments to validate
            
        Returns:
            Tuple[bool, str]: (is_valid, message)
        """
        try:
            if not self.validation_criteria:
                return True, "No validation criteria set"
            
            # Count total functional groups
            total_groups = len(functional_groups)
            max_allowed_groups = self._get_max_groups()
            
            print(f"DEBUG: Alkyl chain {alkyl_chain['name']} has {total_groups} functional groups (max allowed: {max_allowed_groups})")
            
            if total_groups > max_allowed_groups:
                return False, f"Too many functional groups on alkyl chain {alkyl_chain['name']}: {total_groups} (max: {max_allowed_groups})"
            
            return True, f"Valid number of functional groups on chain ({total_groups})"
            
        except Exception as e:
            return False, f"Error validating total functional groups on chain: {str(e)}"

    def _validate_specific_group_counts(self, alkyl_chain: Dict, functional_groups: List[Dict]) -> Tuple[bool, str]:
        """Eq. 18: Validates that number of specific functional group types doesn't exceed their limits
        
        Σyingil ≤ nGal : Sum of functional groups of type i on chain l must not exceed maximum allowed nGal
        
        Args:
            alkyl_chain: The alkyl chain fragment
            functional_groups: List of functional group fragments to validate
            
        Returns:
            Tuple[bool, str]: (is_valid, message)
        """
        try:
            if not self.validation_criteria:
                return True, "No validation criteria set"
            
            # Count occurrences of each functional group type
            group_counts = {}
            for group in functional_groups:
                group_type = group['name']
                group_counts[group_type] = group_counts.get(group_type, 0) + 1
            
            max_group_counts = self._get_max_groups_per_type()
            
            for group_type, count in group_counts.items():
                print(f"DEBUG: Functional group {group_type} has count {count} (max allowed: {max_group_counts})")
                
                if count > max_group_counts:
                    return False, f"Too many groups of type {group_type} on chain: {count} (max: {max_group_counts})"
            
            return True, "Valid functional group type counts on chain"
            
        except Exception as e:
            return False, f"Error validating specific functional group counts: {str(e)}"

    def _validate_group_occurrences_per_chain(self, alkyl_chain: Dict, functional_groups: List[Dict]) -> Tuple[bool, str]:
        """Per-Chain Group Occurrence Constraints (Equations 19-21)
        
        These equations validate the number of occurrences of each group type WITHIN a single chain:
        Eq. 19: Σyingil ≤ t1 (upper limit of a specific group on chain i)
        Eq. 20: Σyingil ≥ t2 (lower limit of a specific group on chain i)
        Eq. 21: Σyingil = t3 (exact count of a specific group on chain i)
        
        Args:
            alkyl_chain: The alkyl chain fragment
            functional_groups: List of functional groups attached to this chain
            
        Returns:
            Tuple[bool, str]: (is_valid, message)
        """
        try:
            if not self.validation_criteria:
                return True, "No occurrence constraints set"
            
            # Count occurrences of each group type on this chain
            group_counts = {}
            for group in functional_groups:
                group_type = group['name']
                group_counts[group_type] = group_counts.get(group_type, 0) + 1
            
            for group_type, count in group_counts.items():
                # Check upper limit (t1) - Eq. 19
                if count > self.validation_criteria.group_occurrence_upper:
                    return False, f"Group {group_type} occurs {count} times on chain {alkyl_chain['name']}, exceeding chain limit of {self.validation_criteria.group_occurrence_upper}"
                
                # Check lower limit (t2) - Eq. 20
                if count < self.validation_criteria.group_occurrence_lower:
                    return False, f"Group {group_type} occurs {count} times on chain {alkyl_chain['name']}, below chain minimum of {self.validation_criteria.group_occurrence_lower}"
                
                # Check exact count (t3) - Eq. 21
                if self.validation_criteria.group_occurrence_exact >= 0 and count != self.validation_criteria.group_occurrence_exact:
                    return False, f"Group {group_type} occurs {count} times on chain {alkyl_chain['name']}, must be exactly {self.validation_criteria.group_occurrence_exact}"
                
                print(f"DEBUG: Group {group_type} occurs {count} times on chain {alkyl_chain['name']} (chain limits: {self.validation_criteria.group_occurrence_lower}-{self.validation_criteria.group_occurrence_upper})")
            
            return True, "Valid group occurrences on chain"
            
        except Exception as e:
            return False, f"Error validating group occurrences on chain: {str(e)}"

    def _validate_total_group_occurrences(self, alkyl_chains: List[Dict], all_functional_groups: List[List[Dict]]) -> Tuple[bool, str]:
        """Aggregate Group Occurrence Constraints (Equations 22-24)
        
        These equations validate the TOTAL number of occurrences of each group type ACROSS ALL chains:
        Eq. 22: ΣΣyingil ≤ t4 (upper limit on total occurrences across all chains)
        Eq. 23: ΣΣyingil ≥ t5 (lower limit on total occurrences across all chains)
        Eq. 24: ΣΣyingil = t6 (exact count of occurrences across all chains)
        
        Args:
            alkyl_chains: List of alkyl chain fragments
            all_functional_groups: List of lists of functional groups (one list per chain)
            
        Returns:
            Tuple[bool, str]: (is_valid, message)
        """
        try:
            if not self.validation_criteria:
                return True, "No total occurrence constraints set"
            
            # Count total occurrences of each group type across all chains
            total_group_counts = {}
            for chain_groups in all_functional_groups:
                for group in chain_groups:
                    group_type = group['name']
                    total_group_counts[group_type] = total_group_counts.get(group_type, 0) + 1
            
            for group_type, total_count in total_group_counts.items():
                # Check total upper limit (t4) - Eq. 22
                if total_count > self.validation_criteria.alkyl_total_upper:
                    return False, f"Group {group_type} occurs {total_count} times total, exceeding total limit of {self.validation_criteria.alkyl_total_upper}"
                
                # Check total lower limit (t5) - Eq. 23
                if total_count < self.validation_criteria.alkyl_total_lower:
                    return False, f"Group {group_type} occurs {total_count} times total, below total minimum of {self.validation_criteria.alkyl_total_lower}"
                
                # Check total exact count (t6) - Eq. 24
                if self.validation_criteria.alkyl_total_exact >= 0 and total_count != self.validation_criteria.alkyl_total_exact:
                    return False, f"Group {group_type} occurs {total_count} times total, must be exactly {self.validation_criteria.alkyl_total_exact}"
                
                print(f"DEBUG: Group {group_type} occurs {total_count} times across all chains (total limits: {self.validation_criteria.alkyl_total_lower}-{self.validation_criteria.alkyl_total_upper})")
            
            return True, "Valid total group occurrences"
            
        except Exception as e:
            return False, f"Error validating total group occurrences: {str(e)}"

    def _validate_cation_wide_constraints(self, cation: Dict, alkyl_chains: List[Dict]) -> Tuple[bool, str]:
        """Cation-wide Group Constraints (Equations 22-24)
        
        These equations control the total number of alkyl chains across the whole cation:
        Eq. 22: ΣΣyingil ≤ t4 (upper limit on total alkyl chains)
        Eq. 23: ΣΣyingil ≥ t5 (lower limit on total alkyl chains)
        Eq. 24: ΣΣyingil = t6 (exact total alkyl chain count)
        
        Args:
            cation: The cation fragment
            alkyl_chains: List of alkyl chain fragments attached to this cation
            
        Returns:
            Tuple[bool, str]: (is_valid, message)
        """
        try:
            if not self.validation_criteria:
                return True, "No cation-wide constraints set"
            
            # Get total alkyl chain count
            total_chains = len(alkyl_chains)

            # Check upper limit (t4) - Eq. 22
            # Note: Applying alkyl_total_upper constraint to the number of chains here
            if total_chains > self.alkyl_total_upper:
                return False, f"Total alkyl chains ({total_chains}) exceeds upper limit of {self.alkyl_total_upper}"

            # Check lower limit (t5) - Eq. 23
            # Note: Applying alkyl_total_lower constraint to the number of chains here
            if total_chains < self.alkyl_total_lower:
                return False, f"Total alkyl chains ({total_chains}) below lower limit of {self.alkyl_total_lower}"

            # Check exact count (t6) - Eq. 24
            # Note: Applying alkyl_total_exact constraint to the number of chains here
            if self.alkyl_total_exact >= 0 and total_chains != self.alkyl_total_exact:
                return False, f"Total alkyl chains ({total_chains}) must be exactly {self.alkyl_total_exact}"

            # print(f"DEBUG: Total alkyl chains: {total_chains} (limits: {self.alkyl_total_lower}-{self.alkyl_total_upper})")
            
            return True, "Valid total alkyl chain count"
            
        except Exception as e:
            return False, f"Error validating total alkyl chains: {str(e)}"

    def validate(self, cation: Dict, anion: Dict, alkyl_chains: List[Dict], functional_groups_per_chain: List[List[Dict]] = None) -> Tuple[bool, str]:
        """
        Validate a combination of fragments
        
        Args:
            cation: The cation fragment
            anion: The anion fragment
            alkyl_chains: List of alkyl chain fragments
            functional_groups_per_chain: List of lists of functional groups (one list per chain)
            
        Returns:
            Tuple[bool, str]: (is_valid, message)
        """
        try:
            # Define log_name early, before it's used in the loop
            log_name = f"{cation.get('name', 'C')}/{anion.get('name', 'A')}/{'/'.join(c.get('name', 'Alk') for c in alkyl_chains)}" # For logging

            # Initialize empty lists if no functional groups provided
            functional_groups_per_chain = functional_groups_per_chain or [[] for _ in alkyl_chains]
            
            # Ensure we have same number of functional group lists as chains
            if len(functional_groups_per_chain) != len(alkyl_chains):
                return False, f"Mismatch between number of chains ({len(alkyl_chains)}) and functional group lists ({len(functional_groups_per_chain)})"
            
            # Flatten all fragments for basic validation
            all_fragments = [cation, anion] + alkyl_chains + [g for chain_groups in functional_groups_per_chain for g in chain_groups]
            
            # Check fragment counts
            is_valid, message = self._validate_anion_cation_counts(all_fragments)
            if not is_valid:
                return False, message

            # Check if we can get properties for all fragments
            for frag in all_fragments:
                props = get_fragment_properties(frag['name'], frag['fragment_type'])
                if not props:
                    return False, f"{frag['fragment_type'].title()} {frag.get('name')} properties not found"
            
            # Validate alkyl chain count based on cation valence
            is_valid, message = self._alkyl_chain_count(all_fragments)
            if not is_valid:
                return False, message
            
            # Validate per-chain constraints
            for i, (chain, chain_groups) in enumerate(zip(alkyl_chains, functional_groups_per_chain)):
                chain_log_name = f"{log_name} Chain {i+1} ({chain.get('name', 'Unknown')})"
                print(f"DEBUG [{chain_log_name}]: Checking Total Groups (Eq 17)...")
                # Validate total groups on this chain
                is_valid, message = self._validate_total_groups_per_chain(chain, chain_groups)
                if not is_valid:
                    print(f"FAILED [{chain_log_name}]: Total Groups - {message}")
                    return False, message
                print(f"PASSED [{chain_log_name}]: Total Groups")

                print(f"DEBUG [{chain_log_name}]: Checking Specific Group Counts (Eq 18)...")
                # Validate specific group counts on this chain
                is_valid, message = self._validate_specific_group_counts(chain, chain_groups)
                if not is_valid:
                    print(f"FAILED [{chain_log_name}]: Specific Group Counts - {message}")
                    return False, message
                print(f"PASSED [{chain_log_name}]: Specific Group Counts")

                print(f"DEBUG [{chain_log_name}]: Checking Group Occurrences (Eq 19-21)...")
                # Validate group occurrences per chain (Eq. 19-21)
                is_valid, message = self._validate_group_occurrences_per_chain(chain, chain_groups)
                if not is_valid:
                    print(f"FAILED [{chain_log_name}]: Group Occurrences - {message}")
                    return False, message
                print(f"PASSED [{chain_log_name}]: Group Occurrences")

                # --- Skipping call to undefined _validate_intrachain_bonding_eq16 ---
                # print(f"DEBUG [{chain_log_name}]: Checking Intra-Chain Bonding (Eq 16)...")
                # # Validate intra-chain bonding (Eq. 16)
                # is_valid, message = self._validate_intrachain_bonding_eq16(chain, chain_groups)
                # if not is_valid:
                #     print(f"FAILED [{chain_log_name}]: Intra-Chain Bonding - {message}")
                #     return False, message
                # print(f"PASSED [{chain_log_name}]: Intra-Chain Bonding")

            # Validate total group occurrences across all chains (Eq. 22-24)
            print(f"DEBUG [{log_name}]: Checking Total Group Occurrences (Eq 22-24 for Func Groups)...")
            # Note: This checks functional groups, not alkyl groups as Eq 22-24 might imply
            is_valid, message = self._validate_total_group_occurrences(alkyl_chains, functional_groups_per_chain)
            if not is_valid:
                print(f"FAILED [{log_name}]: Total Group Occurrences - {message}")
                return False, message
            print(f"PASSED [{log_name}]: Total Group Occurrences")

            # Validate cation-wide constraints (Eq 22-24 applied to # chains - likely misinterpretation)
            print(f"DEBUG [{log_name}]: Checking Cation-Wide Constraints (Eq 22-24 for # Chains)...")
            is_valid, message = self._validate_cation_wide_constraints(cation, alkyl_chains)
            if not is_valid:
                print(f"FAILED [{log_name}]: Cation-Wide Constraints - {message}")
                return False, message
            print(f"PASSED [{log_name}]: Cation-Wide Constraints")

            # --- Missing Eq 15 & 16 validation ---
            # Note: The RDKit valence check added previously was removed as it was too strict for fragment SMILES.

            # log_name defined earlier now

            # --- Simple Charge Check (Proxy for Eq 15) ---
            total_charge = 0
            for frag in all_fragments:
                props = get_fragment_properties(frag['name'], frag['fragment_type'])
                # Use charge from properties if available, otherwise default based on type
                charge = props.get('charge')
                if charge is None:
                     if frag['fragment_type'].lower() == 'cation': charge = 1
                     elif frag['fragment_type'].lower() == 'anion': charge = -1
                     else: charge = 0
                total_charge += charge

            if total_charge != 0:
                 print(f"FAILED [{log_name}]: Total charge is not zero ({total_charge})")
                 return False, f"Invalid combination: Total charge is {total_charge}, expected 0."
            print(f"PASSED [{log_name}]: Total Charge Check (Result: {total_charge})")
            # --- End Charge Check ---

            print(f"PASSED [{log_name}]: All implemented validation checks.")
            return True, "Valid combination"
            
        except Exception as e:
            return False, f"Error validating combination: {str(e)}"
