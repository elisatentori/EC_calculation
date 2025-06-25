Package for Transfer Entropy, Cross-correlation and Cross-covariance
Calculation and analysis

Author: Elisa Tentori <tentorielisa@gmail.com> (University of Padua)
Date:   May 2025


_________________________________________________________________________
_________________________________________________________________________



- Codes are adapted to handle data from the Collaborative Research in 
  Computational Neuroscience (CRCNS) data sharing
  initiative at http://doi.org/10.6080/K0PC308P .

- Some functions are updates of analogues belonging of freely available 
  TE toolbox developed by John Beggs’ group 
  (posted at: http://code.google.com/p/transfer-entropy-toolbox/ ; 
  Ito et al., 2011). These function are targeted as "ITO_TE".


_________________________________________________________________________
_________________________________________________________________________



1. Preparation
   –––––––––––

   Compile the mex file (C program) with gcc or lcc (Matlab default C compiler).



2. Another Spiking Data Format (ASDF)
   ––––––––––––––––––––––––––––––––––

   Another Spike Data Format is basically cell array of spike timing of each neuron.
   In order to calculate TE correctly, I recommend to use only integer for the timing.
   To ensure that, you can do
   >> asdf = ChangeBinning(asdf, 1);

   Last two cells contains special information of the data.

   asdf{end-1}: Binning size in unit of ms (milisecond). (e.g. 1.2 -> 1.2ms/bin, 10 -> 10ms/bin etc...)
   asdf{end}: Array of number of neurons, number of total bins of the data
     (e.g. [10 300000] -> 10 neurons, 300000 time bins)



3. Functions (Please refer to help of each function for details)
   –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

   * TE, SC and XCov Calculation
   ASDFTE_CN.m:      (requires transentmex) delayed higher order TE calculator for ASDF.
   ASDFTE_CN_perm.m: (requires transentmex) delayed higher order TE calculator for ASDF and ASDF2.
   Calculate_TE_CN.m: TE, SC, XCov (and related Z-scored measures) + CI Calculation
   CIReduce.m:        Computes CI for any EC metric
   save_all_measures.m: saves data in .mat files
    
   * Significance test
   NullModel_SignTest.m: jitter test
   
   
   * Changing Data Format
   SparseToASDF.m: Convert matrix form of raster to ASDF.

   * ASDF utilities
   ASDFSubsample.m:     Subsample specified neurons from ASDF.
   ASDFChooseTime.m:    Crop a time segment from larger ASDF.
   ASDFGetfrate.m:      Get firing rate (per bin) of all the neurons.
   ASDFChangeBinning.m: Change binning size of ASDF.


   * Supporting functions (not to be excuted directly)
   transent_CN.c:      Mex file for rapid calculation of TE.
   transent_CN_perm.c: Mex file for rapid calculation of TE between jittered senders and original receivers.
   

   See README_FUNCTIONS_MANUAL.txt for functions description.

_________________________________________________________________________
_________________________________________________________________________


Index of functions:

-) ASDFTE_CN.m:        (adapted from ITO_TE)
-) ASDFTE_CN_perm.m    (adapted from ITO_TE)
-) Calculate_TE_CN.m
-) CIReduce.m          (from ITO_TE)
-) ElecDistance.m
-) NullModel_SignTest.m
-) save_all_measures.m
-) SparseToASDF.m      (from ITO_TE)
-) transent.c          (from ITO_TE)
-) transent_perm.c     (adapted from ITO_TE)

_________________________________________________________________________


Other scripts:

-) demos/                  in which is presented a simple use of some
                           functions in the package
                       
-) how_to_create_list.m    shows how to format data struct for usage

_________________________________________________________________________
_________________________________________________________________________



%=======================================================================%
% Copyright (c) 2025, University of Padua, Italy                        %
% All rights reserved.                                                  %
%                                                                       %
% Authors: Elisa Tentori (elisa.tentori@phd.unipd.it)                   %
%          LiPh Lab - NeuroChip Lab, University of Padua, Italy         %
%=======================================================================%