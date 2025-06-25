Package for Transfer Entropy Calculation and analysis

Author: Elisa Tentori <tentorielisa@gmail.com> (University of Padua)
Date: January 25, 2023


_________________________________________________________________________
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
_________________________________________________________________________



- Codes are adapted to handle data from the Collaborative Research in 
  Computational Neuroscience (CRCNS) data sharing
  initiative at http://doi.org/10.6080/K0PC308P .

- Some functions are updates of analogues belonging of freely available 
  TE toolbox developed by John Beggs’ group 
  (posted at: http://code.google.com/p/transfer-entropy-toolbox/ ; 
  Ito et al., 2011). These function are targeted as "ITO_TE".



Functions are presented in progressive order.



______________________________________________________________________________________________________________________________

         ||||||||||||||||||||||||||||||           ||||||||||||||||||||||||||||||           ||||||||||||||||||||||||||||||
______________________________________________________________________________________________________________________________



Index of functions:
–––––––––––––––––––

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


-) ASDFSubsample.m:
-) ASDFChooseTime.m
-) ASDFGetfrate.m
-) ASDFChangeBinning.m




Other scripts:
––––––––––––––

-) demos/                  in which is presented a simple use of some
                           functions in the package


______________________________________________________________________________________________________________________________

     |||||||||||||||||||||||||||||||||               F U N C T I O N S               ||||||||||||||||||||||||||||||||||||
______________________________________________________________________________________________________________________________




-) ElecDistance.m          [r] = ElecDistance(path_data, b_elecs, h_elecs)
   

   Description:

   Given a MEA, with electrodes disposed as a grid, the function creates a list of possible distances between electrodes couples.


   Parameters:

   path_data   - path to which save the matrix

   b_elecs     - number of electrodes of the base of the rectangle, assuming a rectangular-based MEA

   h_elecs     - number of electrodes of the heigh of the rectangle (see b_elecs)
  

   Returns:

   r           - array with all possible distances in a MEA

   [ABOUT DISTANCES:] writes the list of possible distances between neuron couples in a file at path_data


______________________________________________________________________________________________________________________________


-) Calculate_TE_CN.m          Calculate_TE_CN(path_data, list_cultures, path_results, binsize, maxdelay)


   Description:

   0. Makes a binning of spike times series for each neuron, giving 1/0 series if spike happens/not happens.
   1. Calculates Transfer Entropy (check which type, given order and delay) for a neuron culture recording.
      There's the possibility to calculate TE for a list of recordings related to one or more cultures.
   2. Saves Neurons position coordinates (via channel number).
   3. Calculates matrix of distances between neurons' channels. 


   Parameters:

   path_data          - path to which take the recordings

   list_cultures      - list of cultures to analyze (name_files without format)

   path_results       - path to which save results

   binsize            - (opt) size of bins to binarize spike times series
                        size is intended as number of measure time-steps (0.05ms)
                        [default 20] (20 time-steps=1ms)
   
   maxdelay           - (opt) max time-lag on pre-synaptic neuron to which calculate TE
                        If maxdelay>1, for each couple TE is peakTE: the max[TE(d)], with d in [1,maxdelay]
                        [default 1]


   Returns:

   void funtion

   [ABOUT TE:]                 - saves TE matrices and, if maxdelay>1, delays matrices (with delays that maximize TE for each couple)

   [ABOUT NEURON POSITIONS:]   - saves coordinates of each neuron

   [ABOUT DISTANCES:]          - saves NxN matrix with distance between each neuron couple



______________________________________________________________________________________________________________________________


-) SparseToASDF.m          asdf = SparseToASDF(raster, binunit)


   Description :

   This function converts sparse matrix version of the data to ASDF version.


   Parameters:

   raster    - (N_neu, duration) (Sparse) time raster expressed as a sparse matrix.

   binunit   - (scalar, double) the unit of time in the data in ms. (length of a bin in real scale)


   Returns:

   asdf      - {n_neu + 2, 1} ASDF version of the data



______________________________________________________________________________________________________________________________


-) ASDFTE.m

   [te_result, ci_te_result, te_delay, all_delayed_te, sc_result, ci_sc_result, sc_delay, centered_all_coincs, xcov_result, ci_xcov_result, xcov_delay, all_XCov_d, all_coincidences] = ASDFTE_CN(asdf, j_delay, i_order, j_order, windowsize)


   Description :

   This function manages Transfer Entropy, Signed Cross-Correlation and Cross-Cvariance calculation via transent_CN.c.


 Parameters:

   asdf           - Time series in Another Spike Data Format (ASDF)
   
   j_delay        - Number of bins to lag sender (j series) or a vector [default 1]
   
   i_order        - Order of receiver [default 1]
   
   j_order        - Order of sender [default 1]
   
   windowsize     - window size used for Coincidence Index calculation (odd number only)


 Returns:

   te_result      - (nNeu, nNeu) NxN matrix where N(i, j) is the transfer entropy from i->j
                    If multiple j_delay is given, this is a peak value of TE over delay.
   ci_te_result   - (nNeu, nNeu) NxN matrix where N(i, j) is the Coincidence Index from i->j for TE
                     Multiple delays are necessary to calculate it.
   te_delay       - (nNeu,nNeu) NxN matrix where D(i,j) is the delay that maximizes the TE.
   all_te         - (nNeu, nNeu, delays) NxNxd matrix where N(i, j, k) is the transfer entropy
                    from i->j at delay of jdelay(k). For those who need all the delays.


   sc_result      - (nNeu, nNeu) NxN matrix where N(i, j) is the transfer entropy from i->j
                    If multiple j_delay is given, this is a peak value of SC over delay.
   ci_sc_result   - (nNeu, nNeu) NxN matrix where N(i, j) is the Coincidence Index from i->j for SC
                     Multiple delays are necessary to calculate it.
   sc_delay       - (nNeu,nNeu) NxN matrix where D(i,j) is the delay that maximizes the SC.
   all_sc         - (nNeu, nNeu, delays) NxNxd matrix where N(i, j, k) is the SC
                    from i->j at delay of jdelay(k). For those who need all the delays.

   xc_result      - (nNeu, nNeu) NxN matrix where N(i, j) is the Cross-covariance from i->j
                    If multiple j_delay is given, this is a peak value of XCov over delay.
   ci_xc_result   - (nNeu, nNeu) NxN matrix where N(i, j) is the Coincidence Index from i->j for XCov
                     Multiple delays are necessary to calculate it.
   xc_delay       - (nNeu,nNeu) NxN matrix where D(i,j) is the delay that maximizes the XCov.
   all_xc         - (nNeu, nNeu, delays) NxNxd matrix where N(i, j, k) is the cross-covariance
                    from i->j at delay of jdelay(k). For those who need all the delays.

   all_CN         - (nNeu, nNeu, delays) NxNxd matrix with the delayed coincidences for each jdelay(k).



______________________________________________________________________________________________________________________________


-) transent.c          transent(asdf, pre_delay, post_order, pre_order)


   Description:
   
   C code that effectively does Transfer Entropy calculation for asdf matrix of spike-times (binarized).
   This C code is called by ASDFTE function (should be in the same folder).
 
   NOTES:
   (1) transent.c contains 2 mutually exclusive functions: transent_1(...), transent_ho(...).
   Up to how we are using this library, normally transent_1(...) is called
   For correct usage of transent_ho(...), the entire library TE_matlab has to be appropriately expanded (not huge modifications are requested)


   Parameters:

   asdf.         - Time series in Another Spike Data Format (ASDF)

   pre_delay     - Number of bins to lag sender (asdf series) or a vector [default 1]   [PRE-SYNAPTIC]

   post_order    - Order of receiver [default 1]   [POST-SYNAPTIC]

   pre_order     - Order of sender [default 1]     [PRE-SYNAPTIC]


   Command to compile in Matlab: 
   
   mex transent.c


______________________________________________________________________________________________________________________________


-) NullModel.m          NullModel(path_data, path_nullmodel, list_cultures, time_jitter, binsize, n_permutations, maxdelay)


   Description:

   Directed connections are expected to have both higher TE and CI values than by chance.
   As a null model, we jitter the spike times only from the sender neuron to conserve the auto-prediction of the receiver neuron.

   1. For each culture considerated in list_cultures, the function jitters (randomly) each spike in the source series by using a Gaussian 
   distribution centered in the actual spike time and with a standard deviation of (time_jitter)ms [default time_jitter=10ms]. The Gaussian 
   distribution and the short standard deviation make the analysis stringent.

   2. For each culture considerated in list_cultures, having an asdf matrix jittered (representing neurons time series of sender neurons) and 
   an asdf matrix with original series (representing the receiver neurons), the function calculates TE between sender and receiver neurons.
   
   Points 1. and 2. are iterated n_permutations times [default n_permutations=1000].
   After iteration the TE matrices are calculated.



   Parameters:

   path_data          - path to which take the recordings

   path_nullmodel     - path to which save TE matrices for jittered series

   list_cultures      - list of culture registration filenames

   time_jitter        - (opt) jittering time for each spike time of presynaptic neurons
                        [default 10ms]

   binsize            - (opt) size of bins to binarize spike times series
                        size is intended as number of measure time-steps (0.05ms)
                        [default 20] (20 time-steps=1ms)

   n_permutations     - (opt) number of null-models to generate per each culture
                        [default 1000]

   maxdelay           - (opt) max delay of the pre-synaptic neuron to which calculate TE
                        [default 1]



   Returns:

   void funtion

   [ABOUT NULL MODEL:] saves TE matrix for each permutation (for each culture in list_cultures)



______________________________________________________________________________________________________________________________


-) ASDFTE_perm.m    [te_result, ci_result, all_te, delay_result] = ASDFTE_perm(asdf, asdf2, j_delay > 1, i_order, j_order, windowsize)

   [te_result] = ASDFTE_perm(asdf, asdf2, j_delay = 1, i_order, j_order, windowsize)

   Parameters:

   asdf        - Jittered time series in Another Spike Data Format (ASDF)
                 [asdf = is the jittered matrix referred to neurons as sender neurons]

   asdf2       - Time series in Another Spike Data Format (ASDF)
                 [asdf2 = is the NON-jittered matrix, referred to neurons as receiver neurons]

   j_delay     - Number of bins to lag sender (asdf series) or a vector
                 [default 1]   [PRE-SYNAPTIC or sender]

   i_order     - Order of receiver
                 [default 1]   [POST-SYNAPTIC or receiver]

   j_order     - Order of sender 
                 [default 1]   [PRE-SYNAPTIC or sender]

   windowsize  - window size used for Coincidence Index calculation (odd number only)


   Returns:

   te_result      - (nNeu, nNeu) NxN matrix where N(i, j) is the transfer entropy from i->j
                    If multiple j_delay is given, this is a peak value of TE over delay.
   ci_result      - (nNeu, nNeu) NxN matrix where N(i, j) is the Coincidence Index from i->j
                    Multiple delays are necessary to calculate it. (OPT - it is returned only when j_delay>1  
                                                                          –> see ASDFTE_perm.m header for clarifications)
   all_te         - (nNeu, nNeu, delays) NxNxd matrix where N(i, j, k) is the transfer entropy
                    from i->j at delay of jdelay(k). For those who need all the delays. 
                    (OPT - it is returned only when j_delay>1 –> see ASDFTE_perm.m header for clarifications)
                    
   delay_result   - (nNeu,nNeu) NxN matrix where D(i,j) is the delay that maximizes the TE. 
                    (OPT - it is returned only when j_delay>1 –> see ASDFTE_perm.m header for clarifications)



______________________________________________________________________________________________________________________________


-) transent_perm.c         transent_perm(asdf_pre, asdf_post, pre_delay, post_order, pre_order)


   Description:
   
   (Before, see descriptions of NullModel.m, ASDFTE_perm.m
   Same of transent.c modified in order to calculate TE between pre-synaptic jittered timeseries
   and post-synaptic not-jittered timeseries.
   
   This C code is called by ASDFTE_perm function (should be in the same folder).
 
  
   Parameters:

   asdf_pre      - Jittered time series in Another Spike Data Format (ASDF)    
                   [asdf = is the jittered matrix referred to neurons as sender neurons]

   asdf_post     - Time series in Another Spike Data Format (ASDF)     
                   [asdf2 = is the NON-jittered matrix, referred to neurons as receiver neurons]

   pre_delay     - Number of bins to lag sender (asdf series) or a vector 
                   [default 1]   [PRE-SYNAPTIC]

   post_order    - Order of receiver 
                   [default 1]   [POST-SYNAPTIC]

   pre_order     - Order of sender
                   [default 1]     [PRE-SYNAPTIC]


   
   NOTES:
   (1) In this code x=i is referred to the POST synaptic neurons series of spike-times (receivers),
       While y=j is the PRE synaptic neuron series of spike-times (senders)
   (2) The original transent_1(..) function calculates the TE as
       TE(x(t),x(t-1),y(t-d)). The choice, motivated by the autors in Ito et al.-2011,
       is based on the assumption that the only past of the post-synaptic neuron
       that counts is the one related to its last 1 ms time window.
       Here the function has been modified to calculate TE as
       TE(x(t),x(t-d),y(t-d)). Indeed, we prefer to evaluate also x at time d.
   (3) transent_perm.c contains 2 mutually exclusive functions: transent_1(...), transent_ho(...).
       Up to how we are using this library, normally transent_1(...) is called
       For correct usage of transent_ho(...), the entire library TE_matlab has to be appropriately expanded (not huge modifications are requested)


   Returns:

   TE matrix


   Command to compile in Matlab: 

   mex transent_perm.c

_________________________________________________________________________________________________________________________________


-) ASDFSubsample.m: 

   Subsample specified neurons from ASDF.


_________________________________________________________________________________________________________________________________


-) ASDFChooseTime.m: 

   Crop a time segment from larger ASDF.


_________________________________________________________________________________________________________________________________


-) ASDFGetfrate.m: 

   Get firing rate (per bin) of all the neurons.


_________________________________________________________________________________________________________________________________


-) ASDFChangeBinning.m: 

   Change binning size of ASDF.


_________________________________________________________________________________________________________________________________







%=======================================================================%
% Copyright (c) 2025, University of Padua, Italy                        %
% All rights reserved.                                                  %
%                                                                       %
% Authors: Elisa Tentori (elisa.tentori@phd.unipd.it)                   %
%          LiPh Lab - NeuroChip Lab, University of Padua, Italy         %
%=======================================================================%





