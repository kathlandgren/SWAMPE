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
processes that can take place in exoplanetary atmospheres. Exoplanet scientists work with indirect and limited observations of the planets that they study. In order to form hypotheses about their climate, weather, and atmospheric composition, astronomers need robust models that demonstrate how atmospheres act under different conditions. One-dimensional energy-balance models can capture complex mechanisms such as cloud formation and can rapidly explore the parameter ranges, but they fail to account for variations with longitude. In contrast, three-dimensional models capture the variation in the three-dimensional space, but they are computationally expensive, sometimes taking months to explore the parameter regimes. Their complexity can also obscure the mechanisms that govern atmospheric phenomena. This leaves a natural gap for two-dimensional models, which can capture the spatial variability as well as rapidly explore the parameter space and study the dynamical mechanisms.

`SWAMPE` is a Python package for modeling the dynamics of exoplanetary atmospheres. While tuned to synchronously rotating hot Jupiters and sub-Neptunes, 
the code is modular and could be easily modified to model dissimilar space objects, from Brown Dwarfs to terrestrial, potentially habitable exoplanets. 


# Modeling Exoplanet Atmospheres with SWAMPE

- input planetary parameters: radius, rotation rate, stellar radiation, scale height. 
- SWAMPE solves the shallow-water equations by employing the spectral method outlined in H&J
- the time-stepping scheme is modified Euler, following Langton thesis
- swampe has the capabilty to save data at any desired frequency
- outputs geopotential maps and wind fields, which can be used to make inferences about the temperature profile and the dynamical mechanisms in exoplanet atmospheres

# Statement of need

Current efforts to model exoplanet atmospheres primarily focus on minimal-complexity
one-dimensional and high-complexity three-dimensional models. One-dimensional (1D)
energy-balance models can capture complex mechanisms (e.g. @Bell:2018)
and can rapidly explore the parameter space, but they fail to account for longitudinal variation.
Furthermore, recent observations of giant exoplanets have shown that one-dimensional models
cannot completely describe some of the key atmospheric processes (e.g. \cite{feng2016impact}).
On the other hand, complex three-dimensional (3D) models can capture variation in the physical space. 
They are frequently based on primitive equations (e.g. \cite{menou2009atmospheric},
\cite{kataria2016atmospheric}, \cite{parmentier20133d}) or on the Navier-Stokes equations 
(e.g. \cite{cooper2006dynamics}, \cite{dobbs2013three}) and can be used to understand a variety of radiative,
chemical, and dynamical processes. 3D models such as \textit{ROCKE-3D} \cite{rocke3d} can be tuned to
a variety of exoplanets. However, 3D models tend to be computationally expensive, sometimes taking months
to explore the parameter space. 

The difference in capability between 1D and 3D models leaves a natural gap for two-dimensional
shallow-water models, which can capture the spatial variability as well as run fast enough to
rapidly explore the parameter space and study the dynamical mechanisms. 

In particular, shallow-water models have been used to study solar system planets (including Earth),
e.g. \cite{ferrari2011processes}, \cite{brueshaber2019dynamical}. Outside of the solar system,
shallow-water models have been used to understand a variety of atmospheric phenomena of hot Jupiters,
such as atmospheric variability \cite{menou2003weather} and superrotation \cite{showman2011equatorial}.
They have also been used to make observational predictions for hot Jupiters (e.g. \cite{langton2008hydrodynamic},
\cite{PBS}). However, many of these models are written in Fortran, which makes them difficult to adapt
for the varied needs of exoplanetary science.

`SWAMPE` offers a fully Python, open-source implementation of the spectral method, does not need multiple cores, and is flexible and modular. The code could be easily modified to model dissimilar
space objects, from Brown Dwarfs to terrestrial, potentially habitable exoplanets. 


# Documentation

Documentation for `SWAMPE`, with step-by-step tutorials for research applications, is available at [https://swampe.readthedocs.io/en/latest/](https://swampe.readthedocs.io/en/latest/). 

# Similar tools

[Bell EBM](https://github.com/taylorbell57/Bell_EBM/tree/v1.3) [@Bell:2018] is an energy-balance model.
[MITgcm](https://github.com/MITgcm/MITgcm/) [@Marshall:1997] is an open-source, Fortran-based 3D global circulation model which includes a shallow-water mode. 

# Acknowledgements

EL acknowledges financial support from 

EL thanks Nikole Lewis, Tiffany Kataria, Ryan J. MacDonald, Ishan Mishra, Trevor Foote, [more testers here] for helpful discussions.

# References