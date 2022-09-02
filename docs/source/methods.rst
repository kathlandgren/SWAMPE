Methods
===============

SWAMPE is based on a spectral method for solving shallow-water equations on a sphere
described in `Hack and Jakob (1992) <https://opensky.ucar.edu/islandora/object/technotes:112>`_. 
SWAMPE is designed to replace the previous Fortran implementations of this method.

Time stepping
----------------

While Hack and Jakob (1992) describe a semi-implicit time-stepping scheme, we have found
that modified Euler's method (which has been applied to atmospheric models of exoplanets by 
`Langton
<https://www.proquest.com/docview/304661389?pq-origsite=gscholar&fromopenview=true>`_
yields more stable results. 

Filters
----------------

To ensure numerical stability, SWAMPE applies the following filters:

* a modal-splitting filter as described in Hack and Jakob (1992)
* a sixth-degree hyperviscosity filter. We use the formulation based on Gelb and Gleeson (2001).


Testing
----------------

To ensure the correct operation of the spectral transforms, a series of unit tests are performed 
via continuous integration with Github Actions. 

SWAMPE has been benchmarked against end-to-end tests 1 and 2 from a standard test set for 
numerical shallow-water solvers (see Williamson et al. (199?)).



