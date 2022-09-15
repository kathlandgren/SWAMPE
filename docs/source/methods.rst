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
    :width: 600
    :alt: Illustration of one time step using the spectral method employed by SWAMPE


Governing equations
-------------------

The shallow-water equations on the sphere are given by: 

:math:`\frac{d\mathbf{V}}{dt}=-f\mathbf{k}\times\mathbf V-\nabla\Phi`

:math:`\frac{d\Phi}{dt}=-\Phi\nabla\cdot\mathbf{V}`

where :math:`\mathbf V=u\mathbf{i}+v\mathbf{j}` is the velocity vector on the surface of the sphere and
 :math:`\mathbf i` and :math:`\mathbf j` are the unit eastward and northward directions, respectively. 
  The free surface geopotential is given by :math:`\Phi\equiv gh`, where :math:`g` is the gravitational acceleration. 

The atmosphere is assumed to be a fluid that is incompressible and hydrostatically balanced.
For a derivation of the shallow water equations from the continuity equation and the equation of motion, see, e.g.
`Kaper and Engel (2013) <"https://epubs-siam-org.proxy.library.cornell.edu/doi/book/10.1137/1.9781611972610">`_.


For the spherical harmonic transform method that we use, it is convenient to multiply the velocity 
:math:`\mathbf V` by the cosine of latitude so that the velocity arguments are smooth at the poles. 
We will also rewrite the governing equations in terms of geopotential, the vertical component of relative
 vorticity :math:`\zeta=\mathbf k\cdot(\nabla\times\mathbf V)`, and horizontal divergence :math:`\delta=\nabla\cdot \mathbf V`.

Taking the curl (:math:`\mathbf k\cdot\nabla\times[ \ ]`) and divergence (:math:`\nabla\cdot[ \ ]`) of the momentum equation yields the following:

:math:`\frac{\partial\zeta}{\partial t}=-\nabla\cdot (\zeta+f)\mathbf V\`

:math:`\frac{\partial\delta}{\partial t}=\mathbf k\cdot\nabla\times((\zeta+f)\mathbf V)-\nabla^2\left(\Phi+\frac{\mathbf V\cdot\mathbf V}{2}\right).`

Writing the continuity equation so that the partial time derivative is the only term on the left hand side yields
:math:`\frac{\partial \Phi}{\partial t}=-(\mathbf V\cdot\nabla)\Phi-\Phi\nabla\cdot\mathbf V=-\nabla\cdot(\Phi\mathbf V).`

It is now convenient to use absolute vorticity :math:`\eta=\zeta+f` (as opposed to relative vorticity :math:`\eta=\zeta+f`), as one of the state variables
in addition to state variables divergence :math:`\delta`, and geopotential :math:`\Phi`.
The zonal winds :math:`U` and meridional winds :math:`V`, on the other hand, are diagnostic variables: they will not be time-stepped directly. 
The winds enter the time-stepping scheme via nonlinear components.

We follow `Hack and Jakob (1992) <"https://opensky.ucar.edu/islandora/object/technotes:112">`_ in our notation for the nonlinear terms, 
namely :math:`A=U\eta`, :math:`B=V\eta`, :math:`C=U\Phi`, :math:`D=V\Phi`, :math:`E=\frac{U^2+V^2}{2(1-\mu^2)}`.



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
**Note**
`SWAMPE`'s default hyperviscosity coefficient has been tested for hot Jupiter and sub-Neptune simulations but might require further tuning
for drastically different stellar forcings. The modal-splitting coefficient typically does not need to be adjusted from its default value.

Testing
----------------

To ensure the correct operation of the spectral transforms, a series of unit tests are performed 
via continuous integration with Github Actions. 

`SWAMPE` has been benchmarked against end-to-end tests 1 and 2 from a standard test set for 
numerical shallow-water solvers 
(see `Williamson and Drake (1992) <https://www.sciencedirect.com/science/article/pii/S0021999105800166>`_).
as well as strongly irradiated hot Jupiters described by `Perez-Becker and Showman (2013) <https://ui.adsabs.harvard.edu/abs/2013ApJ...776..134P/abstract>`_.


