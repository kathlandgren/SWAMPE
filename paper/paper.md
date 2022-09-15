---
title: '`SWAMPE`: A Shallow-Water Atmospheric Model in Python for Exoplanets'
tags:
  - Python
  - astronomy
  - exoplanets
authors:
  - name: Ekaterina Landgren
    orcid: 0000-0001-6029-5216
    affiliation: 1
affiliations:
  - name: Center for Applied Mathematics, Cornell University, Rhodes Hall, Ithaca, NY 14853, USA
    index: 1
date: 02 Septempber 2022
bibliography: paper.bib
--- 


# Summary

In order to answer questions about potential habitability of exoplanets, it is important to develop a robust understanding of a variety of dynamic
processes that can take place in exoplanetary atmospheres. While many exoplanets are readily characterized with current facilities like Hubble and the recently-launched James Webb Space Telescope, exoplanet scientists work with indirect and limited observations of the planets that they study. In order to form hypotheses about their climate, weather, and atmospheric composition, astronomers need robust models that demonstrate how atmospheres act under different conditions. One-dimensional energy-balance models can capture complex mechanisms such as cloud formation and can rapidly explore the parameter ranges, but they fail to account for variations with longitude. In contrast, three-dimensional models capture the variation in latitude, longitude, and altitude, but they are computationally expensive, sometimes taking months to explore the parameter regimes. Their complexity can also obscure the mechanisms that govern atmospheric phenomena. This leaves a natural gap for two-dimensional models, which can capture the spatial variability as well as rapidly explore the parameter space and study the dynamical mechanisms.

`SWAMPE` is a Python package for modeling the dynamics of exoplanetary atmospheres. `SWAMPE` is an intermediate-complexity, two-dimensional shallow-water general circulation model. Benchmarked for synchronously rotating hot Jupiters and sub-Neptunes, 
the code is modular and could be easily modified to model dissimilar space objects, from Brown Dwarfs to terrestrial, potentially habitable exoplanets. 


# Modeling Exoplanet Atmospheres with SWAMPE

Exoplanets exist in a vast range of orbital and planetary parameters. `SWAMPE` is designed to be adaptable to a variety of possible regimes. The user can specify physical parameters such as radius, surface gravity, rotation rate, stellar radiation, and scale height. 

`SWAMPE` solves the shallow-water equations using the spectral method [@Hack:1992], with a modified Euler's method timestepping scheme [@Langton:2008]. To ensure numerical-stability, two filters are applied: the modal-splitting filter [@Hack:1992] and a sixth-degree hyperviscosity filter [@Gelb:2001]. `SWAMPE` has the capability to save simulation data at any user-specified frequency. The model outputs geopotential maps and the associated wind fields, which can be used to make inferences about the temperature profiles of exoplanet atmospheres and the dynamical mechanisms behind them.

![Sample `SWAMPE` output: geopotential maps for a hot Jupiter exoplanet at three values of radiative timescale $\tau_{\rm rad}$: 0.1 days, 1 day, and 10 days. \label{fig:SWAMPE-output}](timescale_example-1.png)

# Statement of need

Current efforts to model exoplanet atmospheres primarily focus on minimal-complexity
one-dimensional and high-complexity three-dimensional models. One-dimensional (1D)
energy-balance models can capture complex mechanisms [e.g., @Bell:2018]
and can rapidly explore the parameter space, but they fail to account for longitudinal variation.
Furthermore, recent observations of giant exoplanets have shown that one-dimensional models
cannot completely describe some of the key atmospheric processes [e.g., @Feng:2016].
On the other hand, complex three-dimensional (3D) models are able to capture variation in the physical space. 
They are frequently based on primitive equations [e.g., @Menou:2009; @Kataria:2016;
 @Parmentier:2013] or on the Navier-Stokes equations 
[e.g., @Cooper:2006; @Dobbs-Dixon:2013] and can be used to understand a variety of radiative,
chemical, and dynamical processes. 3D models such as `ROCKE-3D` [@Way:2017] can be tuned to
a variety of exoplanets. However, 3D models tend to be computationally expensive, sometimes taking months
to explore the parameter space. 

The difference in capability between 1D and 3D models leaves a natural gap for two-dimensional
shallow-water models, which can capture the spatial variability as well as run fast enough to
rapidly explore the parameter space and study the dynamical mechanisms. In particular, shallow-water models have been used to study solar system planets, including Earth
[e.g., @Ferrari:2011; @Brueshaber:2019]. Outside of the solar system,
shallow-water models have been used to understand a variety of atmospheric phenomena of hot Jupiters,
such as atmospheric variability [@Menou:2003] and superrotation [@Showman:2011].
They have also been used to make observational predictions for hot Jupiters [e.g., @Langton:2008b;
@Perez-Becker:2013]. However, many of these models are written in Fortran, which makes them difficult to adapt
for the varied needs of exoplanetary science.

`SWAMPE` offers a fully Python, open-source implementation of the 2D shallow-water system. This package does not require multiple cores, and is flexible and modular. `SWAMPE` is designed to be easily modified to model dissimilar space objects, from Brown Dwarfs to terrestrial, potentially habitable exoplanets. `SWAMPE` provides the capability to conduct
wide parameter sweeps and to produce maps of the thermal and wind properties of the planets in latitude and longitude, which can be used to help constrain and make predictions for observations of their atmospheres.

# Documentation

Documentation for `SWAMPE`, with step-by-step tutorials for research applications, is available at [https://swampe.readthedocs.io/en/latest/](https://swampe.readthedocs.io/en/latest/). 

# Similar tools

[Bell EBM](https://github.com/taylorbell57/Bell_EBM/tree/v1.3) [@Bell:2018] is an energy-balance model.
[MITgcm](https://github.com/MITgcm/MITgcm/) [@Marshall:1997] is an open-source, Fortran-based 3D global circulation model which includes a shallow-water mode. 

# Acknowledgements

EL expresses gratitude to the developers of many open-source Python packages used by `SWAMPE`, in particular `numpy` [@Harris:2020], `SciPy` [@Virtanen:2020], and `Matplotlib` [@Hunter:2007].

EL acknowledges financial support from 

EL thanks Nikole Lewis, Tiffany Kataria, Ryan J. MacDonald, Ishan Mishra, Trevor Foote, [more testers here] for helpful discussions.

# References