
Intro:
------------------
#### data

* Generate simulation data under folder "simdata"

* Build Matrices containing Nearest Neighbor informations for NNGP under folder "sorted_x_order"

#### projects

* This folder provides code for response GP and NNGP, latent GP and NNGP. The required stan code are under folder "src"

#### conj_models

* This folder provides code for sampling latent process by conjugate latent NNGP models

#### compare

* Provide code for calculating KL-D and RMSPE


Notes
---------
* This folder contains code files for the simulations in section 4.2. For the sake of simplicity, the replicated code files are cleaned. One may need to modify the code, say "M = 5" to "M = 10", to run all simulations.
* The recommeded Stan code for NNGP is available through [Stan case study](http://mc-stan.org/users/documentation/case-studies/nngp.html)






