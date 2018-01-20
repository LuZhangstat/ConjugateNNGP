
Structure
------------------
+-- data<br />
|   &nbsp;&nbsp;        +-- rawdata: run "CallFunc_ReadNcdf" to import SST<br />
|   &nbsp;&nbsp;        +-- data: run "datapre" to do simple EDA and save useful data<br />
|   &nbsp; &nbsp;       +-- buildNN: Build Matrices containing Nearest Neighbor informations for NNGP through<br />
|<br />
+-- projects: R code files for response NNGP, model1 and model2.<br />
|       &nbsp; &nbsp;     +-- samplew_model1.R & samplew_model1.R. <br />
|       &nbsp;    &nbsp;  &nbsp; &nbsp;    code files for sampling latent process by conjugate latent NNGP models<br />
|       &nbsp;&nbsp;      +-- spNNGPCV.R: code for obtaining point estimates of fixed parameters using cross-validation<br />
|       &nbsp;&nbsp;      +-- functions.R: functions needed in "samplew_model1&2.R"<br />
|       &nbsp;&nbsp;      +-- nngp. R: R code for response NNGP, the corresponding stan codes are under folder "src"<br />
|<br />
+-- src: stan code files for conducting response NNGP, model1 and model2<br />
|<br />
+-- compare:  Code for calculating KL-D and RMSPE for all models<br />


Notes
---------
* The recommeded Stan code for NNGP is available through [Stan case study](http://mc-stan.org/users/documentation/case-studies/nngp.html)





