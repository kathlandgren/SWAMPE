Methods
===============

`SWAMPE` is based on a spectral method for solving shallow-water equations on a sphere
described in `Hack and Jakob (1992) <"https://opensky.ucar.edu/islandora/object/technotes:112">`_. 
`SWAMPE` is designed to replace the previous Fortran implementations of this method.

    "The basic idea behind the spectral transform method is to locally evaluate all nonlinear
    terms (including diabatic physical processes) in physical space on an associated
    finite-difference-like grid, most often referred to as the transform grid. 
    These terms are then transformed back into wavenumber space to calculate 
    linear terms and derivatives, and to obtain tendencies for the time-dependent state variables."

.. image:: https://github.com/kathlandgren/SWAMPE/blob/main/docs/_static/method_illustration.png?raw=true
    :width: 300
    :alt: Illustration of one time step using the spectral method employed by SWAMPE

Governing equations
-------------------

.. math::
    :nowrap:
    \begin{align}
    \frac{d\mathbf{V}}{dt}&=-f\mathbf{k}\times\mathbf V-\nabla\Phi\\
    \frac{d\Phi}{dt}&=-\Phi\nabla\cdot\mathbf{V}
    \end{align}

Time stepping
----------------

While `Hack and Jakob (1992) <"https://opensky.ucar.edu/islandora/object/technotes:112">`_ describe a semi-implicit 
time-stepping scheme, we follow `Langton (2008)
<https://www.proquest.com/docview/304661389?pq-origsite=gscholar&fromopenview=true>`_ and implement a 
modified Euler's method scheme. In practice, this method proves more stable, especially for 
strongly irradiated planets. 

Filters
----------------

To ensure numerical stability, SWAMPE applies the following filters:

* a modal-splitting filter as described in `Hack and Jakob (1992) <"https://opensky.ucar.edu/islandora/object/technotes:112">`_.
* a sixth-degree hyperviscosity filter. We use the formulation based on `Gelb and Gleeson (2001) <https://www.researchgate.net/publication/230675145_Spectral_Viscosity_for_Shallow_Water_Equations_in_Spherical_Geometry>`_.
Note that `SWAMPE`'s default hyperviscosity coefficient has been tested for hot Jupiter and sub-Neptune simulations but might require further tuning
for drastically different stellar forcings. The modal-splitting coefficient typically does not need to be adjusted from its default value.

Testing
----------------

To ensure the correct operation of the spectral transforms, a series of unit tests are performed 
via continuous integration with Github Actions. 

`SWAMPE` has been benchmarked against end-to-end tests 1 and 2 from a standard test set for 
numerical shallow-water solvers 
(see `Williamson and Drake (1992) <https://www.sciencedirect.com/science/article/pii/S0021999105800166>`_).
as well as strongly irradiated hot Jupiters described by `Perez-Becker and Showman (2013) <https://ui.adsabs.harvard.edu/abs/2013ApJ...776..134P/abstract>`_.


