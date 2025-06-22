# AxionDCBH

Jupyter project for Direct Collapse Black Holes via Axion decay. Useful physics constants and definitions are given inside the physics/definitions.py, and the main project is run in the AxionDCBH.ipynb.

Variance2.0 holds the value for variance of cosmic density field as a function of mass. 

Packages needed to run the code: <br />
    1.) Python 3.12 <br />
    2.) SciPy 1.14 <br />
    3.) Numpy 2.1 <br />
    4.) Pandas 2.2 <br />
    5.) Matplotlib 3.9 <br />
    6.) h5py 3.12 <br />

    7.) We will also need to install HaloMod for a statistical analysis of halos. The library can be found at: https://github.com/halomod/halomod

The example.ipynb and physics/definitions.ipynb are just there to test out pieces of code. 'plots' directory contains all the necessary plots, and 'solutions' directory contains the solutions to the parameter space scans. These computations are extemely expensive so it is better to save them in a file and recall from there. 
