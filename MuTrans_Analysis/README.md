# MuTrans workflow

## Overview

This workflow was used for the data preparation and analysis of cell-state transitions using MuTrans.\
See the MuTrans GitHub page (https://github.com/cliffzhou92/MuTrans-release) for installation and usage instructions.

Also see the ICR RSE Group MuTrans GitHub (https://github.com/ICR-RSE-Group/MuTrans-release) for changes made to output plots.

## Structure

1. **Conversion to AnnData** 
   - Separate data and convert to AnnData object at the end of trajectory R scripts

2. **Scanpy pre-processing**
   - Pre-processing and clustering reproduced using `scanpy` in prep scripts.

3. **Landscape construction**
   - Visualisation of the dynamical manifold, identification of attractor basins, and calculation of transition probabilities using MuTrans: *`mutrans_(condition)_dynam.py`*\
     N.B. this needs to be run from the ``./Example/`` folder - see MuTrans GitHub (https://github.com/cliffzhou92/MuTrans-release)
