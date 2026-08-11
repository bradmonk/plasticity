---
layout: default
title: Plasticity — Multiplex Computational Model of Neural Dynamics
---

# Plasticity

**Multiplex Computational Model of Neural Dynamics**

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.19265.svg)](http://dx.doi.org/10.5281/zenodo.19265)
[![GitHub](https://img.shields.io/badge/GitHub-bradmonk%2Fplasticity-blue)](https://github.com/bradmonk/plasticity)

## About the Project

Dynamic changes to the signaling efficacy between neurons are thought to underlie many adaptive behaviors, including memory formation. Signaling efficacies (synaptic weights) are directly correlated with postsynaptic levels of AMPA-type glutamate receptors (AMPAR), the primary receptor mediating fast excitatory transmission.

This project explores factors involved in AMPAR synaptic regulation and presents a unified model for simulating AMPAR trafficking in and around synapses. The spatially-resolved stochastic model integrates experimental data on structural and molecular dynamics, simulating these processes in 3D space.

## Model Components

| Component | Description |
|-----------|-------------|
| **3D Surface Diffusion** | AMPARs diffuse along a 3D surface mesh representing a short dendritic segment with multiple spines, matched to empirical Brownian motion rates. |
| **Actin Scaffold Dynamics** | A dynamic actin filament network is simulated inside spine protrusions, capturing structural plasticity within the postsynaptic density. |
| **Multivalent SAPs** | Synaptic scaffold-associated proteins (SAPs) interact with AMPARs, actin filaments, and each other, forming the molecular machinery of plasticity. |

## Repository Structure

- **`MATmodel/`** — MATLAB implementation of the multiplex actin and receptor diffusion model (STARShiP, mesh-building utilities)
- **`PyDiffusion/`** — Python/FEniCS implementation of particle diffusion on 3D dendritic meshes with Cython-accelerated core routines

## Media

[![3D Mesh Building](https://bradleymonk.com/github/Fenics_Dolphin_3D_Mesh_Building.png)](https://www.youtube.com/watch?v=w8yt6XDgLMU)
[![Python Dendritic Mesh](https://bradleymonk.com/github/Python_Dendritic_Mesh_Building.png)](https://www.youtube.com/watch?v=tDKUU0SqbSA)
[![3D Surface Diffusion](https://bradleymonk.com/github/3D_Model_Surface_Diffusion.png)](https://www.youtube.com/watch?v=TGG_1ypuA4I)
[![Actin Dynamics](https://bradleymonk.com/github/Actin_Dynamics_Spine_Membrane.png)](https://www.youtube.com/watch?v=JH-hGjzhEFQ)

## Further Reading

- [DistMesh — A Simple Mesh Generator in MATLAB](http://persson.berkeley.edu/distmesh/)
- [NeuroML — Neuron Model & Mesh Library](http://www.neuroml.org/tool_support.php)
- [MMBioS — National Center for Multiscale Modeling of Biological Systems](http://mmbios.org/index.php/software)
- [Brownian Motion Simulation Basics](http://bradleymonk.com/Brownian_Motion)
- [Shouval Cluster Model Article (PNAS)](http://www.pnas.org/content/102/40/14440.full)
- [Computational Geometry Algorithm Library (CGAL)](http://www.cgal.org/)
- [Project Documentation Wiki](https://github.com/bradmonk/plasticity/wiki)
