# Applying asymmetric waveguide seismology to data

analysis.py is the main file from which you run your inversion.

The process is made up of two parts:
- time_distance.py creates a class for your data, and has methods to draw out the boundary data.
- fibril_inversion.py creates a class for your boundary data, which can then be used for the inversion procedure.

### Prerequisites

astropy