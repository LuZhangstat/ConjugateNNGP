
Structure
------------------
+-- data<br />
|   &nbsp;&nbsp;&nbsp;&nbsp;        +-- rawdata: run "CallFunc_ReadNcdf" to import SST<br />
|   &nbsp;&nbsp;&nbsp;&nbsp;        +-- data: run "datapre" to do simple EDA and save useful data<br />
|   &nbsp;&nbsp;&nbsp;&nbsp;       +-- buildNN: Build Matrices containing Nearest Neighbor informations for NNGP through<br />
|<br />
+-- projects: R code files for conjugate response NNGP and conjugate latent NNGP.<br />
|   &nbsp;&nbsp;&nbsp;&nbsp;     +-- Conj_RNNGP.R  <br />
|   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;   code for conjugate response NNGP model
|   &nbsp;&nbsp;&nbsp;&nbsp;     +-- CVLNNGP.R  <br />
|   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;   code for cross-validation for conjugate latent NNGP model<br />
|   &nbsp;&nbsp;&nbsp;&nbsp;       +-- Conj_LNNGP_paral.R: code for generating posterior samples of conjugate latent NNGP model <br />
|   &nbsp;&nbsp;&nbsp;&nbsp;       +-- RMSPE_conj_LNNGP.R: code for calculating RMSPE for conjugate latent NNGP model<br />
|   &nbsp;&nbsp;&nbsp;&nbsp;       +-- pic.R: code for drawing pictures in Section V<br />
|<br />
+-- src: code files for building nearest-neighbor matrices, constructing NNGP approximation and calculaing RMSPE <br />


Notes
---------
* The recommeded Stan code for NNGP is available through [Stan case study](http://mc-stan.org/users/documentation/case-studies/nngp.html)





