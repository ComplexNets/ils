# RDKit Data Storage Strategy

## Overview
This document outlines different approaches for handling RDKit molecular property calculations in the ILS project, comparing database storage versus on-the-fly calculation strategies.

## Approach 1: Store RDKit Properties in Database

### Pros
1. Faster retrieval - no need to recalculate properties each time
2. Consistent values across runs - useful for reproducibility
3. Less computational overhead during runtime
4. Good for properties that rarely change or are expensive to calculate

### Cons
1. Data could become stale if RDKit algorithms update
2. Takes up database storage space
3. Need to maintain data consistency
4. Additional complexity in database schema management

## Approach 2: Calculate RDKit Properties On-the-fly

### Pros
1. Always using latest RDKit algorithms
2. No data synchronization needed
3. Simpler database schema
4. More flexible - can easily add new properties
5. Storage efficient

### Cons
1. Higher computational overhead
2. Could slow down runtime performance
3. May get slightly different results across runs (if RDKit updates)

## Recommended Solution: Hybrid Approach

### 1. Store Core Properties
Properties that are:
- Computationally expensive
- Used frequently
- Unlikely to change with RDKit updates
- Critical for reproducibility

Examples:
- molecular_weight
- heavy_atom_count

### 2. Calculate On-the-fly
Properties that are:
- Simple/fast calculations
- Might benefit from RDKit updates
- Less frequently used

Examples:
- rotatable_bond_count
- tpsa

### 3. Implement Caching Layer
- In-memory cache for frequently accessed properties
- Cache refreshes when application starts
- Provides best of both worlds:
  - Fast access
  - Up-to-date calculations

## Benefits of Hybrid Approach
1. Fast access to critical properties
2. Up-to-date calculations when needed
3. Reasonable storage requirements
4. Good balance of performance vs. flexibility

## Implementation Considerations
1. Use database for core molecular properties that are:
   - Expensive to calculate
   - Fundamental to molecule identity
   - Stable across RDKit versions

2. Implement caching system for:
   - Frequently accessed properties
   - Properties used in search/filtering
   - Properties needed for quick comparisons

3. Calculate on-the-fly for:
   - Properties that benefit from latest RDKit updates
   - Less frequently used properties
   - Properties where slight variations are acceptable