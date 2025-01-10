# Shape_Optimization

Integrated (MATLAB + Py) environment performing shape optimization using GA. Upon running the code, this initiates controlpoint coordinates of NACA0012 and then sends this as an input format to the trained Surrogate Model to get pressure and shear stress vectors. Further computes surface integral quantities and then the objective function as well. Moving along in a similar fashion, it generates a set of population size and then generates the function evaluation in a similar fashion eventually dropping the procedure of meshing and solving NS equations for CFD.  
