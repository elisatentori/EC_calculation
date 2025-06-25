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


-) ASDFSubsample.m     (from ITO_TE)
-) ASDFChooseTime.m    (from ITO_TE)
-) ASDFGetfrate.m      (from ITO_TE)
-) ASDFChangeBinning.m (from ITO_TE)




Other scripts:
––––––––––––––

-) demos/                  in which is presented a simple use of some
                           functions in the package


______________________________________________________________________________________________________________________________

     |||||||||||||||||||||||||||||||||               F U N C T I O N S               ||||||||||||||||||||||||||||||||||||
______________________________________________________________________________________________________________________________



______________________________________________________________________________________________________________________________


-) ASDFTE.m

   [te_result, ci_te_result, te_delay, all_delayed_te, sc_result, ci_sc_result, sc_delay, centered_all_coincs, xcov_result, ci_xcov_result, xcov_delay, all_XCov_d, all_coincidences] = ASDFTE_CN(asdf, j_delay, i_order, j_order, windowsize)


    Description :

       This function manages Transfer Entropy, Signed Cross-Correlation and Cross-Cvariance calculation via transent_CN.c.
    ASDFTE_CN(asdf, j_delay, i_order, j_order, windowsize)


    Args:

        asdf             – Time series in Another Spike Data Format (ASDF)
        j_delay          – Number of bins to lag sender (j series) or a vector [default 1]
        i_order          – Order of receiver [default 1]
        j_order          – Order of sender [default 1]
        windowsize       – window size used for Coincidence Index calculation (odd number only)


    Returns:

        te_result        – (nNeurons×nNeurons) matrix of peak transfer entropy (TE) values for each i→j.
        ci_te_result     – (nNeurons×nNeurons) coincidence‐index on all_delayed_te for TE (empty for single‐delay).
        te_delay         – (nNeurons×nNeurons) delay values that maximize TE for each i→j.
        all_delayed_te   – (nNeurons×nNeurons×D) TE at each tested delay (D = num_delays, or num_delays–1 if j_delay(1)==0).

        sc_result        – (nNeurons×nNeurons) signed cross‐correlation peak values for each i→j.
        ci_sc_result     – (nNeurons×nNeurons) coincidence‐index on centered_all_coincs for SC (empty for single‐delay).
        sc_delay         – (nNeurons×nNeurons) delay values that maximize SC for each i→j.
        centered_all_coincs – (nNeurons×nNeurons×(D-1)) SC values at each delay, mean‐centered (zero‐lag excluded).

        xcov_result      – (nNeurons×nNeurons) cross‐covariance peak values for each i→j.
        ci_xcov_result   – (nNeurons×nNeurons) coincidence‐index on XCov_d for XCov (empty for single‐delay).
        xcov_delay       – (nNeurons×nNeurons) delay values that maximize XCov for each i→j.
        XCov_d           – (nNeurons×nNeurons×(D-1)) normalized XCov at each delay (zero‐lag excluded).

        all_coincs       – (nNeurons×nNeurons×D) raw coincidence counts for each delay.

______________________________________________________________________________________________________________________________


-) ASDFTE_perm.m    
    [te_result, ci_te_result, sc_result, ci_sc_result, xcov_result, ci_xcov_result] =
    
    ASDFTE_CN_perm(asdf, asdf2, j_delay, i_order, j_order, windowsize)

   
    Description:
    
         Compute transfer entropy (TE), signed cross-correlation (SC), and cross-covariance (XCov)
         between original and jittered spike-time series via permutation control.
        
    Args:
   
      asdf        - ASDF cell array of original spike-time data
      asdf2       - ASDF cell array of jittered spike-time data
      j_delay     - Vector of sender delays to test (default = 1)
      i_order     - Markov order of receiver (default = 1)
      j_order     - Markov order of sender   (default = 1)
      windowsize  - Odd integer window for coincidence index (default = 5)
   
    Returns:
   
      te_result      - (nNeu×nNeu) peak transfer entropy values for each i→j
      ci_te_result   - (nNeu×nNeu) coincidence index for TE (empty for single-delay)
      sc_result      - (nNeu×nNeu) peak signed cross-correlation values for each i→j
      ci_sc_result   - (nNeu×nNeu) coincidence index for SC (empty for single-delay)
      xcov_result    - (nNeu×nNeu) peak cross-covariance values for each i→j
      ci_xcov_result - (nNeu×nNeu) coincidence index for XCov (empty for single-delay)

______________________________________________________________________________________________________________________________


-) Calculate_TE_CN.m          Calculate_TE_CN(path_data, list_cultures, path_results_TE, path_results_SC, path_results_XCov, binsize, maxdelay, CI_tau_window)


   Description:

       0. Makes a binning of spike times series for each neuron, giving 1/0 series if spike happens/not happens.
       1. Calculates Transfer Entropy (check which type, given order and delay) and the related CI for a neuron 
          culture recording.
          Calculates also the z-scored versions of the metrics.
       2. Calculates Signed Cross-Correlation and Cross-covariance with the related CI and z-scored versions.

        There's the possibility to calculate TE, SC, XCov for a list of recordings related to one or more cultures.

       4. Saves Neurons position coordinates (via channel number).
       5. Calculates matrix of distances between neurons' channels. 


    Args:

     path_data           - path to which take the recordings
     
     list_cultures       - list of cultures to analyze (name_files without format)
     
     path_results_TE     - path to which save results for TE
     
     path_results_SC     - path to which save results for SC
     
     path_results_XCov   - path to which save results for XCov
     
     binsize             - (opt) size of bins to binarize spike times series
                                 size is intended as number of measure time-steps (0.05ms)
                                 [default 20] (20 time-steps=1ms)
                                 
     maxdelay            - (opt) max delay (Unit: number of time-bins) 
                                 of the pre-synaptic neuron at which to calculate TE, SC and XCov
                                 [default 1]
                                 
     CI_tau_window       - (opt) window length for CI (Unit: number of time-bins)


    Returns:

     [] - void funtion

    Saved files:
   
      ABOUT TE:          saves TE matrix, CI matrix, delay matrix (delay that maximizes TE),
                         all_TE matrices (i.e. TE matrices each for fixed delay)
                         
      ABOUT SC:          saves SC matrix, CI matrix, delay matrix (delay that maximizes SC),
                         all_SC matrices (i.e. SC matrices each for fixed delay)
                         
      ABOUT XCov:        saves XCov matrix, CI matrix, delay matrix (delay that maximizes XCov),
                         all_XCov matrices (i.e. XCov matrices each for fixed delay)
                         
      ABOUT DISTANCES:   saves NxN matrix with distance between each neuron couple

______________________________________________________________________________________________________________________________


-) CIReduce.m          result = CIReduce(cells, args)

    Description:

        This function calculates the Coincidence Interval (CI) for a given Effective Connectivity (EC) link i–>j 
        profile as a function of delay. It computes the ration between the sum of the EC values within a time window 
        centered on the global peak of the profile and the sum of EC contained in the whole delays domain.

    Args:

        cells    - 1D Effective Connectivity profile of a link i–>j across time delays (e.g., from a 
                   delay-dependent measure like Transfer Entropy or Cross-Correlation).
   
        args     - Window size (ci_window) in number of bins around the EC peak, to compute the CI. 
                    >> Interval around the peak: (-ci_window/2,ci_window/2)

    Returns:

        result   - Ratio between the EC within the coincidence interval and the total EC across all delays.

______________________________________________________________________________________________________________________________


-) ElecDistance.m          [r] = ElecDistance(path_data, b_elecs, h_elecs)
   

   Description:

       Given a MEA, with electrodes disposed as a 8x8 grid, the function creates a list of possible distances 
       between electrodes couples.
       To be used only with CRCNS data.


   Args:

       path_data   - path to which save the matrix

       b_elecs     - number of electrodes of the base of the rectangle, assuming a rectangular-based MEA (8x8)
    
       h_elecs     - number of electrodes of the heigh of the rectangle (see b_elecs)
  

   Returns:

       r           - array with all possible distances in a MEA

       [ABOUT DISTANCES:] writes the list of possible distances between neuron couples in a file at path_data

______________________________________________________________________________________________________________________________


-) NullModel_SignTest.m

 [ ] = NullModel_SignTest(path_data, list_cultures, time_jitter, binsize, n_permutations, maxdelay,CI_tau_window)
   
   
    Parameters:
   
        path_data          - path to which take the recordings
        
        list_cultures      - list of culture registration filenames
        
        time_jitter        - (opt) jittering time for each spike time of presynaptic neurons
                                   [default 10ms]
                                   
        binsize            - (opt) size of bins to binarize spike times series
                                   size is intended as number of measure time-steps (Unit: 0.05ms)
                                   [default 20] (binsize=20 <==> 20 time-steps=1ms)
                                   
        n_permutations     - (opt) number of null-models to generate per each culture
                                   [default 1000]
                                   
        maxdelay           - (opt) max delay (Unit: number of time-bins) 
                                   of the pre-synaptic neuron at which to calculate TE, SC and XCov
                                   [default 1]
                                   
        CI_tau_window      - (opt) window length for CI (Unit: number of time-bins)
   
   
    Returns:
   
      void funtion
   
   
    Stores:
   
      for each EC metric: P-value matrix (and related counts matrix), Zscored_EC matrix (based on jittered test)
   
______________________________________________________________________________________________________________________________
 
-) save_all_measures.m          save_all_measures(outdir, basename, suffixes, mats, txt_delim)

    Description:
    
        Stores a collection of matrices both .mat (decomment for .txt).
    
    Args:
    
       outdir            -  directory in cui salvare i file
       basename          -  prefisso comune (es. list_cultures{num})
       suffixes          -  cell array di stringhe, es:
                            {'TE_delay1','TE_delay2',...,'SC_delay1',...,'CN_delay0',...}
       mats              -  cell array di matrici, in corrispondenza a suffixes
       txt_delim         -  delimiter per writematrix (es. '\t')

______________________________________________________________________________________________________________________________


-) SparseToASDF.m          asdf = SparseToASDF(raster, binunit)


   Description :

   This function converts sparse matrix version of the data to ASDF version.


   Args:

   raster    - (N_neu, duration) (Sparse) time raster expressed as a sparse matrix.

   binunit   - (scalar, double) the unit of time in the data in ms. (length of a bin in real scale)


   Returns:

   asdf      - {n_neu + 2, 1} ASDF version of the data



______________________________________________________________________________________________________________________________


-) transent.c          transent(asdf, pre_delay, post_order, pre_order)


   Description:
   
       C code that effectively computes Transfer Entropy and the delayed coincidences (for SC and XCov) 
       for asdf matrix of spike-times (binarized).
       This C code is called by ASDFTE function (should be in the same folder).
 
   NOTES:
   (1) transent.c contains 2 mutually exclusive functions: transent_1(...), transent_ho(...).
       Up to how we are using this library, normally transent_1(...) is called
       For correct usage of transent_ho(...), the entire library SC_TE_matlab has to be appropriately expanded
       (not huge modifications are requested)


   Args:

       asdf          - Time series in Another Spike Data Format (ASDF)

       pre_delay     - Number of bins to lag sender (asdf series) or a vector [default 1]   [PRE-SYNAPTIC]

       post_order    - Order of receiver [default 1]   [POST-SYNAPTIC]

       pre_order     - Order of sender [default 1]     [PRE-SYNAPTIC]

   Returns:

       TE matrix and spike coincidence matrix
       
   Command to compile in Matlab:      >> mex transent.c

______________________________________________________________________________________________________________________________


-) transent_perm.c         transent_perm(asdf_pre, asdf_post, pre_delay, post_order, pre_order)


   Description:
   
   (Before, see descriptions of NullModel.m, ASDFTE_perm.m
   Same of transent.c modified in order to calculate TE and delayed coincidences between pre-synaptic jittered 
   timeseries and post-synaptic not-jittered timeseries.
   
   This C code is called by ASDFTE_perm function (should be in the same folder).
 
  
   Args:

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
           For correct usage of transent_ho(...), the entire library TE_matlab has to be appropriately expanded 
           (not huge modifications are requested)


   Returns:

       TE matrix and spike coincidence matrix


   Command to compile in Matlab:      >> mex transent_perm.c
______________________________________________________________________________________________________________________________



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





