# Trajectory inference

## Overview

This workflow was used for the pseudotime ordering and downstream gene analysis of scRNA-seq data.

## Structure

1. **Pseudotime ordering** 
   - Using the `slingshot` package to infer pseudotime variables and construct smooth lineages

2. **Dynamic genes: linear model**
   - Significant difference testing with respect to pseudotime using `testPseudotime()`

3. **Dynamic genes: random forest**
   - Use gene expression to predict pseudotime value with random forests
