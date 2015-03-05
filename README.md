# plasticity



## Project Summary

Dynamic changes to the signaling efficacy between neurons is thought to underlie many adaptive behaviors, including memory formation. Signaling efficacies (i.e. synaptic weights) are directly correlated with postsynaptic levels of AMPA-type glutamate receptors (AMPAR), the primary receptor used in fast excitatory transmissions. While synaptic regulation of AMPAR is considered fundamental to brain information storage, and despite the large number of studies done on this topic, there remains uncertainty regarding AMPAR-mediated synaptic plasticity. However, there is a wealth of information on AMPAR trafficking mechanisms from which to develop a unified model of synaptic plasticity. Such a model could be useful for addressing long-standing questions regarding how temporary signals induce persisting changes to signaling efficacy, and how these synaptic weights are maintained over long durations. This project explores factors involved in AMPAR synaptic regulation, and includes the synthesis of a unified model for simulating AMPAR trafficking in-and-around synapses. Primary components of the model include: (1) AMPARs that diffuse along a 3D dendritic surface, (2) a dynamic actin filament network, and (3) multivalent synaptic/scaffold-associates proteins (SAPs) that can interact with AMPARs, actin filaments, and other SAPs. This spatially-resolved stochastic model integrates experimental data on structural and molecular dynamics, and simulates these processes in 3D space.


## Graphical Abstract
<a href="http://bradleymonk.com/media2/figs/FIG1_ABC.png" target="_blank">
<img src="http://bradleymonk.com/media2/figs/FIG1_ABC.png" width="600" border="10" /></a>
* Overview of modeling environment. (A) Membrane diffusion of AMPAR was simulated on a 3D surface-mesh representing a short dendritic segment with several spines. Geometric values for surface components included dendrite segment length and diameter, density and distance between spines, and dimensions of spine head and neck. AMPAR diffused laterally along these surfaces with Brownian motion; rates were matched to empirical values. (B) Actin scaffolding dynamics were simulated inside the spine-protrusions. Components of the scaffold and their interactin with surface receptors is illustrated in (C).



## Media

<a href="http://youtu.be/9ipFHrxqLWc" target="_blank"><img src="http://camk2.com/pix/yt4.png" 
alt="3d actin network" width="300" border="10" /></a> <a href="http://youtu.be/QBx8F_5_y0g" target="_blank"><img src="http://camk2.com/pix/yt1.png" alt="3d actin network" width="300" border="10" /></a>


<a href="http://youtu.be/JH-hGjzhEFQ" target="_blank"><img src="http://camk2.com/pix/yt2.png" 
alt="3d actin network" width="300" border="10" /></a> <a href="http://youtu.be/t9Vzcvne40w" target="_blank"><img src="http://camk2.com/pix/yt3.png" alt="3d actin network" width="300" border="10" /></a> 


<a href="http://youtu.be/TK4iSQlOOHU" target="_blank"><img src="http://www.bradleymonk.com/w/images/7/73/Matdiffusion.png" 
alt="3d Diffusion Animation in matlab" width="300" border="10" /></a> <a href="http://bradleymonk.com/media2/vid9/Dendritic_Diffusion_3D.mp4" target="_blank"><img src="http://bradleymonk.com/media2/vid9/Dendritic_Diffusion_3D.png" 
alt="Dendritic surface diffusion" width="300" border="10" /></a>



### Analytical proof for steadystate distribution of particles
<a href="http://bradleymonk.com/media2/figs/SSproof1.png" target="_blank"><img src="http://bradleymonk.com/media2/figs/SSproof1.png" width="650" border="10" /></a>
<a href="http://bradleymonk.com/media2/figs/SSproof2.png" target="_blank"><img src="http://bradleymonk.com/media2/figs/SSproof2.png" width="650" border="10" /></a>



## Resources
* [djhbrm github.io site](http://subroutines.github.io/djhbrm/)
* [neuroml.org - neuron model/mesh library](http://www.neuroml.org/tool_support.php)
* [mcell.psc.edu - mcell diffusion software](http://www.mcell.psc.edu/tutorials/ficks_laws.html)
* [mcell on github](https://github.com/mcellteam/cellblender_tutorials)
* [mmbios.org - National Center for Multiscale Modeling of Biological Systems](http://mmbios.org/index.php/software)
* [ppm.org - particle mesh library](http://www.ppm-library.org/)
* [Max Plank Mosaic Group (creators of PPM)](http://mosaic.mpi-cbg.de/?q=research/gallery)
* [Brownian motion simulation basics](http://bradleymonk.com/Brownian_Motion)
* [Shouval cluster model article](http://www.pnas.org/content/102/40/14440.full)
* [Computational Geometry Algorithm Library](http://www.cgal.org/)
* [Kinetic Monte Carlo Simulations](http://www.roentzsch.org/SurfDiff/index.html)
* [Unified Form-assembly Code (UFC) User Manual][1]
* [Q&A about deterministic face indices][2]

[1]: http://fenicsproject.org/pub/documents/ufc/ufc-user-manual/ufc-user-manual.pdf
[2]: http://fenicsproject.org/qa/3233/are-face-and-other-indexes-deterministic

