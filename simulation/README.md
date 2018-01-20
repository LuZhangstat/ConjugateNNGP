
Structure
------------------

|Folder Name |   subfolder    |          |
|:------ |:----------- |:----------- |
|data|simdata: |Generate simulation data by "simdata.R"|
|       |buildNN: |Build Matrices containing Nearest Neighbor informations for NNGP through "buildNNMatrics.R"|
|projects:|<td colspan=2> R code files for response full GP and NNGP, latent full GP and NNGP. |
|               |<td colspan=2>The corresponding stan codes are under folder "src"|
|asym_sim|  <td colspan=2> Simulations in section 4.3 "Comparison of response and latent process models"|
|SST_study| <td colspan=2>Sea Surface Temperature (SST) data Analysis |




#### data
|       +--     simdata: Generate simulation data by "simdata.R"
|       +--     buildNN: Build Matrices containing Nearest Neighbor informations for NNGP through "buildNNMatrics.R"
|
#### projects: R code files for response full GP and NNGP, latent full GP and NNGP.
|                   The corresponding stan codes are under folder "src"
|
+-- src: stan code files for response full GP and NNGP, latent full GP and NNGP
|
+-- conj_models:
|       +-- samplew_model1.R & samplew_model1.R:
|                code files for sampling latent process by conjugate latent NNGP models
|       +-- functions.R: functions needed in "samplew_model1&2.R"
|       +-- spNNGPCV.R: code for obtaining point estimates of fixed parameters using cross-validation
|
+-- compare:  Code files for calculating KL-D and RMSPE of all models


#### data

* Generate simulation data under folder "simdata"

* Build Matrices containing Nearest Neighbor informations for NNGP under folder "sorted_x_order"

#### projects

* This folder provides code for response full GP and NNGP, latent GP and NNGP. The required stan code are under folder "src"

#### conj_models

* This folder provides code for sampling latent process by conjugate latent NNGP models

#### compare

* Provide code for calculating KL-D and RMSPE


Notes
---------






