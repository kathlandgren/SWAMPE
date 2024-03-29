.. SWAMP-E documentation master file, created by
   sphinx-quickstart on Thu Apr 29 14:41:33 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

SWAMPE's documentation
===================================

SWAMPE is a 2D shallow-water general circulation model designed for simulating
exoplanet atmospheres. SWAMPE is a fully in Python implementation of a spectral
algorithm described in `Hack and Jakob (1992) <https://opensky.ucar.edu/islandora/object/technotes:112>`_
and previously available only in Fortran.

SWAMPE's features currently include:

* Exoplanet atmosphere simulations that can be run on your laptop.
* Tuned for synchronously rotating hot Jupiters and sub-Neptunes. 
* Can be easily adapted to model a variety of substellar objects.

This public release contains tutorials for running atmospheric simulations,
working with SWAMPE output, and plotting.

SWAMPE is available under the BSD 3-Clause License. 

.. image:: https://github.com/kathlandgren/SWAMPE/blob/main/docs/_static/testgif.gif?raw=true
   :width: 300
   :alt: Example geopotential and wind gif generated with SWAMPE

Contributing
==============

To report a bug or request a new feature, please open a `new issue <https://github.com/kathlandgren/SWAMPE/issues/new>`_.

If you would like to get more actively involved in SWAMPE development, you are invited to get in touch ek672-at-cornell-dot-edu,
and I will be happy to help you understand the codebase as you work on your contribution.

.. toctree::
   :maxdepth: 1
   :hidden:

   installation
   methods

.. toctree::
   :maxdepth: 2
   :caption: Tutorials:
   
   notebooks/QuickStart
   notebooks/Plots
   
.. toctree::
   :maxdepth: 2
   :caption: Code Documentation:

   continuation
   filters
   forcing
   initial_conditions
   model
   plotting
   spectral_transform
   time_stepping
   modEuler_tdiff



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
