
Structure
------------------
+-- data<br />
|   &nbsp;&nbsp;        +-- (step 1)simdata: Generate simulation data by "simdata.R"<br />
|   &nbsp;&nbsp;       +-- (step 2)buildNN: Build Matrices containing Nearest Neighbor informations for NNGP through<br />
|    &nbsp; &nbsp;  &nbsp;   &nbsp;         "buildNNMatrics.R"<br />
|<br />
+-- projects: R code files for full GP, response NNGP and latent NNGP.<br />
|  
|<br />
+-- src: code files for building nearest-neighbor matrices, constructing NNGP approximation and calculaing RMSPE <br />
|<br />
+-- compare:  Code for calculating KL-D and RMSPE for all models<br />


Notes
---------
* The recommeded Stan code for NNGP is available through [Stan case study](http://mc-stan.org/users/documentation/case-studies/nngp.html)





