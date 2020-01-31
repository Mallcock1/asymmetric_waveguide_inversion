# Applying asymmetric waveguide seismology to data

Uses the amplitude ratio method from [here](https://iopscience.iop.org/article/10.3847/1538-4357/aaad0c/meta).

Used for the data analysis of Section 4 of [this](https://www.frontiersin.org/articles/10.3389/fspas.2019.00048/full) paper.

analysis.py is the main file from which you run your inversion.

The process is made up of two parts:
- time_distance.py creates a class for your data, and has methods to draw out the boundary data.
- fibril_inversion.py creates a class for your boundary data, which can then be used for the inversion procedure.

### Prerequisites

astropy