

## TBD


### Things to implement/code
* Add scaffold cluster model at synaptic submembrane

======

### Things to analyze
* Time the mean spine escape latency for particles with a starting origin at the synapse (top of spine)


======
## Notes
Make sure Matplotlib, Numpy, etc. are up to date on Mac OS X
so we can properly do animations in 3D:

http://www.tapir.caltech.edu/~dtsang/python.html

Using `dolfin`:

http://fenicsproject.org/download/

Need to:
$ source /Applications/FEniCS.app/Contents/Resources/share/fenics/fenics.conf

```
pip install networkx
```


## Web Resources
* [djhbrm github.io site](http://subroutines.github.io/djhbrm/)
* [neuroml.org - neuron model/mesh library](http://www.neuroml.org/tool_support.php)
* [mcell.psc.edu - mcell diffusion software](http://www.mcell.psc.edu/tutorials/ficks_laws.html)
* [mcell on github](https://github.com/mcellteam/cellblender_tutorials)
* [mmbios.org - National Center for Multiscale Modeling of Biological Systems](http://mmbios.org/index.php/software)
* [ppm.org - particle mesh library](http://www.ppm-library.org/)
* [Max Plank Mosaic Group (creators of PPM)](http://mosaic.mpi-cbg.de/?q=research/gallery)
* [Brownian motion simulation basics](http://www.advancedlab.org/mediawiki/index.php/Simulating_Brownian_Motion)
* [Shouval cluster model article](http://www.pnas.org/content/102/40/14440.full)
* [Computational Geometry Algorithm Library](http://www.cgal.org/)
* [Kinetic Monte Carlo Simulations](http://www.roentzsch.org/SurfDiff/index.html)
* [Unified Form-assembly Code (UFC) User Manual][1]
* [Q&A about deterministic face indices][2]

[1]: http://fenicsproject.org/pub/documents/ufc/ufc-user-manual/ufc-user-manual.pdf
[2]: http://fenicsproject.org/qa/3233/are-face-and-other-indexes-deterministic


## Media
#### Diffusion to steadystate
<a href="http://bradleymonk.com/media/vid6/steadystate.mov" target="_blank"><img src="http://bradleymonk.com/media/vid6/steadystate.png" 
alt="Video on diffusion to steadystate" width="400" border="10" /></a>


#### Example animation rendered from "particle_diffusion_on_mesh.py"
<a href="http://bradleymonk.com/media/vid6/3dDiffusion.html" target="_blank"><img src="http://bradleymonk.com/media/vid6/3dDiffusion.png" 
alt="3d Diffusion Animation in dolfin" width="400" border="10" /></a>


#### Dendritic surface diffusion (on 3D mesh) from "particle_diffusion_on_mesh.py"
<a href="http://bradleymonk.com/media/vid9/Dendritic_Diffusion_3D.mp4" target="_blank"><img src="http://bradleymonk.com/media/vid9/Dendritic_Diffusion_3D.png" 
alt="Dendritic surface diffusion" width="400" border="10" /></a>


#### Example rendering dynamic actin network from "bran_actin_geo_xyz.spi"
<a href="http://bradleymonk.com/media/vid8/RenderActin.mp4" target="_blank"><img src="http://bradleymonk.com/media/vid8/SPiMactin.png" 
alt="3d actin network" width="400" border="10" /></a>



## Figures

#### Fig 1. Cluster model diagrams and synopsis
<a href="http://bradleymonk.com/media/figs/FIG1.png" target="_blank">
<img src="http://bradleymonk.com/media/figs/FIG1.png" width="600" border="10" /></a>
* Fig 1. Cluster model diagrams and synopsis. (A) In the Shouval model, particles form clusters through interactions with their nearest neighbors. Particles dissociate stochastically from the cluster at a constant rate, while the insertion probability for an unoccupied lattice space depends on the number of particles flanking that area. This is computed by convolving the cluster with a nearest-neighbor mask, and evaluating insertion probability based on these field weights. (B) The Czöndör model simulated the lateral 2D trajectories of AMPARs based on experimentally determined diffusion rates that were mapped to extrasynaptic (Dout), perisynaptic (Din), and synaptic (Dpsd), representations. (C) We incorporated key features from these two models into a single multilayer model where synaptic clustering of scaffold protein were simulated using the Shouval formulation and AMPAR surface diffusion was simulated using parameters from the Czöndör model. In addition, (a) we allowed scaffold particles to interact with nearby AMPARs effecting their diffusion rate; in turn, these interactions could impart stabilizing effects on the cluster by increasing lattice insertion propensity in the vicinity under a surface AMPAR. Geometric values for this model included: (b) dendritic segment length, (c) synapse spacing, (d) synaptic area, (e) spine area. 


#### Fig 2. Cluster math
<a href="http://bradleymonk.com/media/figs/FIG2.png" target="_blank">
<img src="http://bradleymonk.com/media/figs/FIG2.png" width="600" border="10" /></a>
* Fig 2. Simulation of particle clustering was conducted using the 2-step equation above [1]. (A) The first step of the on-rate equation (shown in the bottom right of the graph) yields conditional probabilities based on the number of neighbors surrounding a given lattice position (the off-rate equation works similarly). In this formulation each lattice position can have from zero to four neighbors, as displayed on the x-axis. The scalars L and β  are used to weight the cumulative effects of multiple neighbors on particle-insertion probability. When the number of neighboring particles exceed L, the conditional probability for insertion is greater than 0.5; β  adjusts the slope of this function such that higher values of β  result in steeper curves (changes in conditional probability) around L. (B) In the second part of the equation (full equation shown in top left), rate parameters are applied to the conditional probabilities. This figure is showing the various effects of R at the two values of β  indicated above (here L=2 and ∆t=0.1 across conditions). 


#### Fig 3. AMPARs can stabilize scaffold protein clusters
<a href="http://bradleymonk.com/media/figs/FIG3.png" target="_blank">
<img src="http://bradleymonk.com/media/figs/FIG3.png" width="700" border="10" /></a>
* Fig 3. AMPARs can stabilize scaffold protein clusters. (A) Scaffold protein clusters may dissociate over time unless stabilized by AMPA receptors. In this configuration, a starting cluster of 64 scaffold proteins will steadily degrade over the span of 1 hour [blue]; here however, cluster size is fully maintained through interactions with surface receptors [red]. (B) Scaffold cluster breakdown has a reciprocal effect on AMPAR synaptic numbers. As scaffold proteins become less available for tethering, synaptic expression of surface receptors declines. When cluster size nears the zero-point, surface receptor expression will parallel the perisynaptic steady-state. [values averaged over 10 simulations].


#### Fig 4. Scaffold protein cluster growth can support synaptic potentiation
<a href="http://bradleymonk.com/media/figs/FIG4.png" target="_blank">
<img src="http://bradleymonk.com/media/figs/FIG4.png" width="700" border="10" /></a>
* Fig 4. Scaffold protein cluster growth can support synaptic potentiation. (A) During a brief (5 min) time-window the tethering coefficient for one subtype of multivalent surface receptors was transiently increased at ‘synapse-1’, allowing it to assist in scaffold protein recruitment. This event induced scaffold protein accretion along the edges of the cluster, which were stabilized by neighboring scaffold particles. Interestingly, this resulted in a long-term increase in scaffold protein cluster size that persisted well after the transient event subsided. (B) The addition of new submembrane protein scaffolding was closely followed by increases in AMPAR synaptic expression. The synaptic quantity of surface particles tends to parallel scaffold protein cluster size, as such, the stable increase of cluster size resulted in a long-term increase of synaptic AMPAR levels.




#### Fig 5. Scaffold cluster model diagrams and synopsis
<a href="http://bradleymonk.com/media/figs/FIG5.png" target="_blank">
<img src="http://bradleymonk.com/media/figs/FIG5.png" width="600" border="10" /></a>
* Fig 5. Scaffold cluster model diagrams and synopsis. (A) In the scaffold cluster model, SAP particles (red circles) cluster around nodes of actin scaffolding (yellow stars) at the PSD. SAP particles can associated or dissociate stochastically from actin nodes or each other with neighbor-dependent probability rates (similar to the cluster model described above; see SI methods). Actin nodes in this 2D visualization represent the tips of (B) actin scaffolding located in close proximity to the postsynaptic submembrane. This lattice (blue box) is of special consideration because it’s near enough to the submembrane to interact with laterally-diffusing surface receptors. (C) In the full model, AMPAR particles diffused laterally on the surface of a 3D dendritic mesh; in parallel, scaffold-cluster dynamics were simulated in the dendritic spine compartments. SAP particles at the near-membrane lattice could interact with AMPARs diffusing in close proximity, effecting their diffusion rate. In some simulations, AMPAR-SAP interactions could briefly alter actin dynamics in the surrounding vicinity, increasing the probability of de novo actin node formation nearby. 


