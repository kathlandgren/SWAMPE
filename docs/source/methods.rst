Methods
============

SWAMPE is based on a spectral method for solving shallow-water equations on a sphere
described in `Hack and Jakob (1992)<https://opensky.ucar.edu/islandora/object/technotes:112>`_. 

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


Gelb and Gleeson
Williamson



