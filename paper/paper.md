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

# Statement of need
In order to answer questions about potential habitability of exoplanets,
it is important to develop a robust understanding of a variety of dynamic
processes that can take place in exoplanetary atmospheres. 

Current efforts to model exoplanet atmospheres primarily focus on minimal-complexity
one-dimensional and high-complexity three-dimensional models. One-dimensional (1D)
energy-balance models can capture complex mechanisms (e.g. \cite{bell2018increased})
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

There is a rich heritage of varying complexity atmospheric models that have been developed.
In particular, shallow-water models have been used to study solar system planets (including Earth),
e.g. \cite{ferrari2011processes}, \cite{brueshaber2019dynamical}. Outside of the solar system,
shallow-water models have been used to understand a variety of atmospheric phenomena of hot Jupiters,
such as atmospheric variability \cite{menou2003weather} and superrotation \cite{showman2011equatorial}.
They have also been used to make observational predictions for hot Jupiters (e.g. \cite{langton2008hydrodynamic},
\cite{PBS}). However, many of these models are written in Fortran, which makes them difficult to adapt
for the varied needs of exoplanetary science.

SWAMP-E is built fully in Python, and is flexible and modular. The code could be easily modified to model dissimilar
space objects, from Brown Dwarfs to terrestrial, potentially habitable exoplanets. 


# Documentation

Documentation for `SWAMPE`, with step-by-step tutorials illustrating research applications, is available at [https://swampe.readthedocs.io/en/latest/](https://swampe.readthedocs.io/en/latest/). 

# Similar tools

# Acknowledgements

# References