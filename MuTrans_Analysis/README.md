# MuTrans workflow

## Overview

This workflow was used for the data preparation and analysis of cell-state transitions using MuTrans.\
See the MuTrans GitHub page (https://github.com/cliffzhou92/MuTrans-release) for installation and usage instructions.

Also see the ICR RSE Group MuTrans GitHub (https://github.com/ICR-RSE-Group/MuTrans-release) for changes made to output plots.

## Structure

1. **Conversion to AnnData** 
   - Separate data by condition (untreated, cisplatin, JQ1) and convert to AnnData object: *`data_prep.R`*

2. **Scanpy pre-processing**
   - Pre-processing and clustering reproduced using `scanpy`: *`mutrans_(condition).ipynb`*

3. **Landscape construction**
   - Visualisation of the dynamical manifold, identification of attractor basins, and calculation of transition probabilities using MuTrans: *`mutrans_(condition).ipynb`*\
     N.B. this needs to be run from the ``./Example/`` folder - see MuTrans GitHub (https://github.com/cliffzhou92/MuTrans-release)
   - Gene set enrichment analysis to characterise the attractor basins carried out in R: *`gene_analysis.R`*
