# Effective Connectivity calculation

Package for Transfer Entropy, Cross-correlation and Cross-covariance Calculation and analysis

Author: Elisa Tentori (University of Padua)

Date:   May 2025

## Contents:

<ul>
  <li>
    <a href="https://github.com/elisatentori/EC_calculation/tree/main/SC_TE_packageP" target="_blank">  SC_TE_package</a>: Package for Transfer Entropy, Cross-correlation and Cross-covariance Calculation. <br>See below for a detailed description.
  </li>
  <li>
    <a href="https://github.com/elisatentori/EC_calculation/tree/main/demo_TE_SC_calculation" target="_blank">  demo_TE_SC_calculation</a>: Package application on spiking data from neural cultures. <br>See below for a detailed description.
  </li>
</ul>


## SC_TE_package

<ul>
<li> Codes are adapted to handle data from the <a href="http://doi.org/10.6080/K0PC308P" target="_blank">  Collaborative Research in Computational Neuroscience (CRCNS)</a> data sharing initiative. </li>
 <li> Some functions are extentions of the ones belonging of freely available TE toolbox developed by John Beggs’ group: <a href="http://code.google.com/p/transfer-entropy-toolbox/" target="_blank"> transfer-entropy-toolbox</a>; Ito et al., 2011. In the list below, these function are targeted as "ITO_TE".
</ul>


### 1. Preparation

<ul>
<li>Compile the mex file (C program) with gcc or lcc (Matlab default C compiler).</li>
</ul>

### 2. Another Spiking Data Format (ASDF)

<ul>
<li>
Another Spike Data Format is basically cell array of spike timing of each neuron. <br>In order to calculate TE correctly, we recommend to use only integer for the timing. <br>To ensure that, you can do >> asdf = ChangeBinning(asdf, 1);
<br>

Last two cells contains special information of the data.
  <ul>
     <li>asdf{end-1}: Binning size in unit of ms (milisecond). (e.g. 1.2 -> 1.2ms/bin, 10 -> 10ms/bin etc...)</li>
     <li>asdf{end}: Array of number of neurons, number of total bins of the data<br>
       (e.g. [10 300000] -> 10 neurons, 300000 time bins)</li>
  </ul>

</li>
</ul>

### 3. Functions

<ul>
  <li><b>TE, SC and XCov Calculation:</b><br>
  <a href="https://github.com/elisatentori/EC_calculation/blob/main/SC_TE_package/ASDFTE_CN.m" target="_blank"> ASDFTE_CN.m</a>:      (requires transentmex) delayed higher order TE calculator for ASDF.<br>
  <a href="https://github.com/elisatentori/EC_calculation/blob/main/SC_TE_package/ASDFTE_CN_perm.m" target="_blank"> ASDFTE_CN_perm.m</a>: (requires transentmex) delayed higher order TE calculator for ASDF and ASDF2.<br>
   <a href="https://github.com/elisatentori/EC_calculation/blob/main/SC_TE_package/Calculate_TE_CN.m" target="_blank"> Calculate_TE_CN.m</a>: TE, SC, XCov (and related Z-scored measures) + CI Calculation<br>
   <a href="https://github.com/elisatentori/EC_calculation/blob/main/SC_TE_package/CIReduce.m" target="_blank"> CIReduce.m</a>:        Computes CI for any EC metric<br>
   <a href="https://github.com/elisatentori/EC_calculation/blob/main/SC_TE_package/save_all_measures.m" target="_blank"> save_all_measures.m</a>: saves data in .mat files 
   </li> <br>
    
   <li><b>Significance test:</b><br>
   <a href="https://github.com/elisatentori/EC_calculation/blob/main/SC_TE_package/NullModel_SignTest.m" target="_blank"> NullModel_SignTest.m</a>: jitter test
   </li><br>
   
  <li><b>Changing Data Format:</b> <br>
   <a href="https://github.com/elisatentori/EC_calculation/blob/main/SC_TE_package/SparseToASDF.m" target="_blank"> SparseToASDF.m</a>: Convert matrix form of raster to ASDF.
  </li><br>
  
   <li> <b>ASDF utilities:</b> <br>
   <a href="https://github.com/elisatentori/EC_calculation/blob/main/SC_TE_package/ASDFSubsample.m" target="_blank"> ASDFSubsample.m</a>:     Subsample specified neurons from ASDF.<br>
   <a href="https://github.com/elisatentori/EC_calculation/blob/main/SC_TE_package/ASDFChooseTime.m" target="_blank"> ASDFChooseTime.m</a>:    Crop a time segment from larger ASDF.<br>
   <a href="https://github.com/elisatentori/EC_calculation/blob/main/SC_TE_package/ASDFGetfrate.m" target="_blank"> ASDFGetfrate.m</a>:      Get firing rate (per bin) of all the neurons.<br>
   <a href="https://github.com/elisatentori/EC_calculation/blob/main/SC_TE_package/ASDFChangeBinning.m" target="_blank"> ASDFChangeBinning.m</a>: Change binning size of ASDF.
   </li><br>

 <li>
   <b>Supporting functions:</b><br>
   <a href="https://github.com/elisatentori/EC_calculation/blob/main/SC_TE_package/transent_CN.c" target="_blank"> transent_CN.c</a>:      Mex file for rapid calculation of TE. <br>
   <a href="https://github.com/elisatentori/EC_calculation/blob/main/SC_TE_package/transent_CN_perm.c" target="_blank"> transent_CN_perm.c</a>: Mex file for rapid calculation of TE between jittered senders and original receivers.
 </li> <br>
</ul>

   See <a href="https://github.com/elisatentori/EC_calculation/blob/main/SC_TE_package/README_FUNCTIONS_MANUAL.txt" target="_blank"> README_FUNCTIONS_MANUAL.txt</a> for functions description.

__

### Index of functions:

-) ASDFTE_CN.m:        (adapted from ITO_TE) <br>
-) ASDFTE_CN_perm.m    (adapted from ITO_TE) <br>
-) Calculate_TE_CN.m <br>
-) CIReduce.m          (from ITO_TE) <br>
-) ElecDistance.m <br>
-) NullModel_SignTest.m <br>
-) save_all_measures.m <br>
-) SparseToASDF.m      (from ITO_TE) <br>
-) transent.c          (from ITO_TE) <br>
-) transent_perm.c     (adapted from ITO_TE) <br>


## Demo

<ul>
  <li>Datasets used:<br>
   - <a href="https://github.com/elisatentori/EC_calculation/tree/main/demo_TE_SC_calculation/Data_CRCNS" target="_blank"> Data_CRCNS</a>: 5 minutes recordings (Spontaneous activity, 2 samples) from the <a href="http://doi.org/10.6080/K0PC308P" target="_blank"> CRCNS</a> data sharing initiative
  <br>
   - <a href="https://github.com/elisatentori/EC_calculation/tree/main/demo_TE_SC_calculation/Data_MaxOne" target="_blank"> Data_MaxOne</a>: our lab acquisitions (2 samples). Spontaneous activity, from rat embryo hippocampal cultures coupled with <a href="https://www.mxwbio.com/products/maxone-mea-system-microelectrode-array/" target="_blank"> MaxOne Single-Well HD-MEA System</a>.
  </li>

  <li> <a href="https://github.com/elisatentori/EC_calculation/tree/main/demo_TE_SC_calculation/scripts_demo" target="_blank"> scripts_demo</a>:<br>
    - <a href="https://github.com/elisatentori/EC_calculation/blob/main/demo_TE_SC_calculation/scripts_demo/TECN_Calculation.m"> TECN_Calculation.m</a>: computes EC (all metrics)<br>
    - <a href="https://github.com/elisatentori/EC_calculation/blob/main/demo_TE_SC_calculation/scripts_demo/TECN_Sign.m"> TECN_Sign.m</a>: significance test<br>
    - <a href="https://github.com/elisatentori/EC_calculation/blob/main/demo_TE_SC_calculation/scripts_demo/job.sh"> job.sh</a>: to launch the scripts with bash/SLURM
  </li>

  <li> <a href="https://github.com/elisatentori/EC_calculation/tree/main/demo_TE_SC_calculation/Notebooks_ECoutputs" target="_blank"> Notebooks_ECoutputs</a>:<br>
    -  <a href="https://github.com/elisatentori/EC_calculation/blob/main/demo_TE_SC_calculation/Notebooks_ECoutputs/ongoing_activity_and_EC_visualization.ipynb" target="_blank"> ongoing_activity_and_EC_visualization.ipynb</a>: to manage spike-trains and visualize EC
    -  <a href="https://github.com/elisatentori/EC_calculation/blob/main/demo_TE_SC_calculation/Notebooks_ECoutputs/utils" target="_blank"> utils</a>: functions useful for anaysis
  </li>
</ul>

%=======================================================================%
% Copyright (c) 2025, University of Padua, Italy                        %
% All rights reserved.                                                  %
%                                                                       %
% Authors: Elisa Tentori (elisa.tentori@phd.unipd.it)                   %
%          LiPh Lab - NeuroChip Lab, University of Padua, Italy         %
%=======================================================================%
