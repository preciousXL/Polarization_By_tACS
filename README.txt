Code associated Huang XL, Wei XL, Wang J, Yi GS. (2024). "Frequency-dependent membrane polarization across neocortical cell types and subcellular elements by transcranial alternating current stimulation". (Under Review)

This repository contains the model set of cortical neurons and the Matlab/Python code for reproducing figures.

Steps for reproducing the paper results:

1. Run addpaths.m
2. Run RunSim_Figure1.m ~ RunSim_Figure5and6.m to produce data
3. Run snowp_figure/plotfigure1-7.m to preprocess simulation data for plotting figures
4. Run Plot_figures.ipynb in Jupyter Notebook to plot Figs. 1-9
Note: Selecting appropriate ion dynamics for Figs. 8 and 9 by revising the code in "nrn/cell_chooser.hoc" (lines 22-24)