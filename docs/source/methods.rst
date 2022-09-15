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


The modified Euler's method for differential equations has the form:

:math:`y_{n+1}=y_n+(\Delta t/2)[f(t_n,y_n)+f(t_{n+1}, y_n+\Delta t f(t_n,y_n)].`

Equivalently, we can write:

:math:`K_1=\Delta t f(t_n,y_n),`

:math:`K_2=\Delta t f(t_{n+1},y_n+K_1),`

:math:`y_{n+1}=y_n+\frac{K_1+K_2}{2}.`

We will now derive these coefficients for the state variables in our timestepping scheme.

Firstly, we can split the unforced shallow-water system into the linear and nonlinear components by rewriting it as follows:


:math:`\frac{d}{dt} \begin{bmatrix} \eta^m_n \\ \delta^m_n \\ \Phi^m_n
\end{bmatrix} =  \begin{bmatrix}
0 & 0 & 0\\
0 & 0 & \frac{n(n+1)}{a^2} \\
0 & -\overline{\Phi} & 0
\end{bmatrix}
\begin{bmatrix}
\eta^m_n \\
\delta^m_n \\
\Phi^m_n
\end{bmatrix}
+ \begin{bmatrix}
\mathscr{E} (t)\\
\mathscr{D} (t)\\
\mathscr{P} (t)
\end{bmatrix}`

where :math:`\begin{bmatrix}
\mathscr{E} (t)\\
\mathscr{D} (t)\\
\mathscr{P} (t)
\end{bmatrix}` represents the nonlinear time-dependent components.
We will evaluate the first component of the right hand side implicitly, while evaluating the second component explicitly.

The *unforced* nonlinear components can be expressed as follows:

:math:`\mathscr{E}(t)=-\frac{1}{a(1-\mu^2)}\frac{\partial A}{\partial \lambda}-\frac{1}{a}\frac{\partial B}{\partial \mu}`

:math:`\mathscr{D}(t)=\frac{1}{a(1-\mu^2)}\frac{\partial B}{\partial \lambda}-\frac{1}{a}\frac{\partial A}{\partial \mu}-\nabla^2E`

:math:`\mathscr{P}(t)=-\frac{1}{a(1-\mu^2)}\frac{\partial C}{\partial \lambda}-\frac{1}{a}\frac{\partial D}{\partial \mu}.`

Let :math:`F_{\Phi}`$` be the geopotential forcing (for `SWAMPE`, due to stellar irradiation, but more general in theory). 
Let :math:`F_{U}=F_{u}\cos \phi` and :math:`F_{V}=F_{v}\cos \phi` be momentum forcing. Then the *forced* nonlinear components are as follows:

:math:`\mathscr{E}(t)=-\frac{1}{a(1-\mu^2)}\frac{\partial} {\partial \lambda}(A-F_{V})-\frac{1}{a}\frac{\partial }{\partial \mu}(B+F_{U}),`

:math:`\mathscr{D}(t)=\frac{1}{a(1-\mu^2)}\frac{\partial }{\partial \lambda}(B+F_{U})-\frac{1}{a}\frac{\partial }{\partial \mu}(A-F_{V})-\nabla^2E,`

:math:`\mathscr{P}(t)=-\frac{1}{a(1-\mu^2)}\frac{\partial C}{\partial \lambda}-\frac{1}{a}\frac{\partial D}{\partial \mu}+ F_{\Phi}.`


Following the notation of the modified Euler's method, we write :math:`K^1=\Delta t f(t,y_t)`:


:math:`K^1_{\eta}=\Delta t (\mathscr{E} (t)),`

:math:`K^1_{\delta}=\Delta t \left(\dfrac{n(n+1)}{a^2}\Phi^{m(t)}_n+\mathscr{D} (t)\right),`

:math:`K^1_{\Phi}=\Delta t \left(-\overline{\Phi}\delta^{m(t)}_n+\mathscr{P} (t)\right).`

Then we can write the :math:`K^2=\Delta t (f(t+1,y_t+K^1))` coefficients. 

:math:`K^2_{\eta}=\Delta t (\mathscr{E} (t+1)),`

:math:`K^2_{\delta}=\Delta t \left(\mathscr{D} (t+1) +\dfrac{n(n+1)}{a^2}(\Phi^m_n+K^1_{\Phi})\right),`


:math:`K^2_{\Phi}=\Delta t \left(\mathscr{P} (t+1)-\overline{\Phi}(\delta^m_n+K^1_{\delta})\right).`


Expanding the equations for :math:`K^2_{\delta}` and :math:`K^2_{\Phi}`, we obtain:

:math:`K^2_{\delta}=\Delta t \left(\mathscr{D} (t+1) +\dfrac{n(n+1)}{a^2}(\mathscr{P}(t))+\dfrac{n(n+1)}{a^2}\Phi^m_n-\overline{\Phi}\dfrac{n(n+1)}{a^2}\delta^m_n \right),`

:math:`K^2_{\Phi}=\Delta t \left(\mathscr{P} (t+1)-\overline{\Phi}(\mathscr{D}(t))-\overline{\Phi}\delta^m_n-\overline{\Phi}\dfrac{n(n+1)}{a^2} \Phi^m_n\right).`


We evaluate the time-dependent terms explicitly, assuming
:math:`
    \begin{bmatrix}
\mathscr{E} (t)\\
\mathscr{D} (t)\\
\mathscr{P} (t)
\end{bmatrix}=
\begin{bmatrix}
\mathscr{E} (t+1)\\
\mathscr{D} (t+1)\\
\mathscr{P} (t+1)
\end{bmatrix}`
to first order. This is what is done in the semi-implicit method in \citet{hack1992description}. An alternative variant would be to approximate $\eta$, $\delta$, $\Phi$, $U$, and $V$ by a different method, such as forward Euler's method or a semi-implicit one. This would result in a higher computational cost and hopefully higher accuracy as well, while maintaining the stability properties of modified Euler's method. 

Note that in the current implementation, :math:`\eta` time-stepping is equivalent to forward Euler's method, since :math:`\eta` does not depend linearly on other state variables, only nonlinearly in the :math:`\mathscr{E}(t)` term.  
Writing :math:`(K^1+K^2)/2` in order to evaluate the modified Euler scheme, we can simplify:

:math:`\dfrac{K^1_{\delta}+K^2_{\delta}}{2}=\Delta t\left( \dfrac{n(n+1)}{a^2} \Phi^m_n +\mathscr{D}(t) + \dfrac{1}{2}\left(   \dfrac{n(n+1)}{a^2}(\mathscr{P}(t) -\overline{\Phi} \delta^m_n   \right)\right),`

and 

:math:`\dfrac{K^1_{\Phi}+K^2_{\Phi}}{2}=\Delta t\left( -\overline{\Phi}\delta^m_n +\mathscr{P}(t)\right)-\dfrac{\Delta t}{2}\overline{\Phi} \left( \mathscr{D}(t)+\dfrac{n(n+1)}{a^2} \right).`


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


