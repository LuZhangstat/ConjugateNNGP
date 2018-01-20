
Structure
------------------
+-- data<br />
|   &nbsp;&nbsp;        +-- simdata_num:
|   &nbsp;&nbsp;   &nbsp;&nbsp;  &nbsp;&nbsp;  &nbsp;&nbsp;  +-- simdata_num.R:
Generate simulation data<br />
|   &nbsp;&nbsp;   &nbsp;&nbsp;  &nbsp;&nbsp;  &nbsp;&nbsp;  +-- buildNN_num.R:
Build Matrices containing Nearest Neighbor informations for NNGP<br />
|<br />
+-- projects: R code files for conducting response and latent NNGP <br />
|    &nbsp; &nbsp; &nbsp;  &nbsp;      The corresponding stan codes are under folder "src"<br />
|<br />
+-- src: stan code files for response and latent NNGP<br />

Notes
---------
* This folder contains code files for the simulations in section 4.2. For the sake of simplicity, the replicated code files are cleaned. One may need to modify the code, say "M = 5" to "M = 10", to run all simulations.
* The recommeded Stan code for NNGP is available through [Stan case study](http://mc-stan.org/users/documentation/case-studies/nngp.html)






