.. SWAMP-E documentation master file, created by
   sphinx-quickstart on Thu Apr 29 14:41:33 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

SWAMPE's documentation
===================================

SWAMPE is a 2D shallow-water general circulation model designed for simulating
exoplanet atmospheres.

SWAMPE's features currently include:

* Simulations that can be run on your laptop.
* Can be easily adapted to model a variety of substellar objects.

This public release contains tutorials for running atmospheric simulations,
working with SWAMPE output, and plotting.

SWAMPE is available under the BSD 3-Clause License. 


.. toctree::
   :maxdepth: 1
   :hidden:

   installation

.. toctree::
   :maxdepth: 2
   :caption: Tutorials:
   
   notebooks/QuickStart
   notebooks/Plotting
   
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


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
