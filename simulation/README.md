
Structure
------------------
+-- data<br />
|   &nbsp;&nbsp;        +-- simdata: Generate simulation data by "simdata.R"<br />
|   &nbsp; &nbsp;       +-- buildNN: Build Matrices containing Nearest Neighbor informations for NNGP through<br />
|    &nbsp; &nbsp;  &nbsp;   &nbsp;         "buildNNMatrics.R"<br />
|<br />
+-- projects: R code files for response full GP and NNGP, latent full GP and NNGP.<br />
|    &nbsp; &nbsp; &nbsp;  &nbsp;      The corresponding stan codes are under folder "src"<br />
|<br />
+-- src: stan code files for conducting response full GP and NNGP, latent full GP and NNGP<br />
|<br />
+-- conj_models:<br />
|       &nbsp;&nbsp;      +-- functions.R: functions needed in "samplew_model1&2.R"<br />
|       &nbsp;&nbsp;      +-- spNNGPCV.R: code for obtaining point estimates of fixed parameters using cross-validation<br />
|       &nbsp; &nbsp;     +-- samplew_model1.R & samplew_model1.R. <br />
|       &nbsp;    &nbsp;  &nbsp; &nbsp;    code files for sampling latent process by conjugate latent NNGP models<br />
|<br />
+-- compare:  Code for calculating KL-D and RMSPE for all models<br />


Notes
---------
* The recommeded Stan code for NNGP is available through [Stan case study](http://mc-stan.org/users/documentation/case-studies/nngp.html)





