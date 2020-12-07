# THERMINATOR 2
## THERMal heavy IoN generATOR version 2

##### For releases to 2.0.3 inclusive: 
### Mikolaj Chojnacki, Adam Kisiel, Wojciech Florkowski, Wojciech Broniowski

We present an extended version of THERMINATOR, a Monte Carlo event generator dedicated to studies of the statistical production of particles in relativistic heavy-ion collisions. The increased functionality of the code contains the following features: The input of any shape of the freeze-out hypersurface and the expansion velocity field, including the 3+1 dimensional profiles, in particular those generated externally with various hydrodynamic codes. The hypersufraces may have variable thermal parameters, which allows for studies departing significantly from the mid-rapidity region, where the baryon chemical potential becomes large. We include a library of standard sets of hypersurfaces and velocity profiles describing the RHIC Au+Au data at sqrt(s_(NN)) = 200 GeV for various centralities, as well as those anticipated for the LHC Pb+Pb collisions at sqrt(s_(NN)) = 5.5 TeV. A separate code, FEMTO-THERMINATOR, is provided to carry out the analysis of femtoscopic correlations which are an important source of information concerning the size and expansion of the system. We also include several useful scripts that carry out auxiliary tasks, such as obtaining an estimate of the number of elastic collisions after the freeze-out, counting of particles flowing back into the fireball and violating causality (typically very few), or visualizing various results: the particle p_T-spectra, the elliptic flow coefficients, and the HBT correlation radii. We also investigate the problem of the back-flow of particles into the hydrodynamic region, as well as estimate the elastic rescattering in terms of trajectory crossings. The package is written in C++ and uses the CERN ROOT environment. 

Manual can be found on arXiv: [1102.0273](https://arxiv.org/abs/1102.0273)

##### Extensions beyond release 2.0.3, conftributed by
### Malgorzata Gumberidze, Szymon Harabasz, Radoslaw Ryblewski
###### Spherically-symmetric blast-wave freeze-out model
A new freeze-out model of a spherically-symmetric Siemens-Rasumussen description of freeze-out is implemented, which is better applicable to heavy-ion collisions in the few-GeV energy regime. The details and results of the model are described in:
Szymon Harabasz, Wojciech Florkowski, Tetyana Galatyuk, Malgorzata Gumberidze, Radoslaw Ryblewski, Piotr Salabura, and Joachim Stroth, [Phys. Rev. C 102, 054903](https://journals.aps.org/prc/abstract/10.1103/PhysRevC.102.054903).

The model is switched on by selecting `SR` (standing for Siemens-Rasmussen) in the file `events.ini`, the freeze-out parameters are set in the file `fomodel/SR.ini`.

###### New output format compatible with UNIGEN
It is now possible to save the generated events in the [UNIGEN](https://www.gsi.de/work/wissenschaftliche_netzwerke/helmholtz_virtuelle_institute/unigen.htm) format. This is a ROOT-based binary format, widely used for lower-energy heavy-ion collisions. The advantages:
* Thanks to the common format, the very same analysis scripts can be used as for the outputs from other models
* In constrast to the standard THERMINATOR output format, in UNIGEN each event has direct access to the list of produced particles, which makes analysis of correlations easier.

To use the UNIGEN output format one needs:
* UNIGEN installed in the system
* Build therminator providing the path to the UNIGEN installation, e.g.:
`make DIR_UNIGEN=/home/myuser/unigen/v2.1/`
* In the `events.ini` file set the output format to `unigen` instead of `root` or `text`
