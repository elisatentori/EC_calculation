% [ ] = NullModel_SignTest(path_data, list_cultures, time_jitter, binsize, n_permutations, maxdelay,CI_tau_window)
%
%
% Parameters:
%
%   path_data          - path to which take the recordings
%   list_cultures      - list of culture registration filenames
%   time_jitter        - (opt) jittering time for each spike time of presynaptic neurons
%                        [default 10ms]
%   binsize            - (opt) size of bins to binarize spike times series
%                        size is intended as number of measure time-steps (Unit: 0.05ms)
%                        [default 20] (binsize=20 <==> 20 time-steps=1ms)
%   n_permutations     - (opt) number of null-models to generate per each culture
%                             [default 1000]
%   maxdelay           - (opt) max delay (Unit: number of time-bins) 
%                              of the pre-synaptic neuron at which to calculate TE, SC and XCov
%                              [default 1]
%   CI_tau_window      - (opt) window length for CI (Unit: number of time-bins)
%
%
% Returns:
%
%   void funtion
%   [ABOUT NULL MODEL:] saves TE matrix for each permutation (for each culture in list_cultures)


%==============================================================================%
% Copyright (c) 2024, University of Padua, Italy                               %
% All rights reserved.                                                         %
%                                                                              %
% Authors: Elisa Tentori (elisa.tentori@phd.unipd.it)                          %
%          LiPh Lab - NeuroChip Lab, University of Padua, Italy                %
%==============================================================================%



function NullModel_SignTest(path_data,path_results_TE,path_results_SC,path_results_XCov,list_cultures,time_jitter,binsize,n_permutations,maxdelay,CI_tau_window)


% ========= Set default arguments ========= %

if nargin<6
    time_jitter=10;
    n_permutations=1000;
    binsize=20;
    maxdelay=1;
    CI_tau_window=5;
end

if nargin<7
    n_permutations=1000;
    binsize=20;
    maxdelay=1;
    CI_tau_window=5;
end

if nargin<8
    binsize=20;
    maxdelay=1;
    CI_tau_window=5;
end

if nargin<9
    maxdelay=1;
    CI_tau_window=5;
end

if nargin<10
    CI_tau_window=5;
end




list_to_load=path_data+list_cultures+".mat";


for num=1:length(list_cultures)
    
    data=load(list_to_load(num));
    
    
    %============= binarizing original spike-trains =============%
    
    binIdx = cellfun(@(sp) reshape(floor(sp / binsize) + 1, [], 1),  data.spikes, 'UniformOutput', false);

    rows = repelem((1:data.nNeurons)', cellfun(@numel, binIdx));
    cols = vertcat(binIdx{:});
    nbins = max([1; cols]);

    binarized = sparse(rows, cols, true, data.nNeurons, nbins);

    % spike-trains of all neurons (seen as post-synaptic)
    asdf_post = SparseToASDF(binarized, 1);  %%= turns into cell array
    
    clear binarized;
    
    %=============        EC measures loading        =============%
    
    t    = load(path_results_TE+list_cultures(num)+"_TEPk.mat");
    TE   = t.peakEC; clear t;
    t    = load(path_results_SC+list_cultures(num)+"_SCPk.mat");
    SC   = t.peakEC; clear t;
    t    = load(path_results_XCov+list_cultures(num)+"_XCovPk.mat");
    XCov = t.peakEC; clear t;

    t       = load(path_results_TE+list_cultures(num)+"_TE_CI.mat");
    ci_TE   = t.ci; clear t;
    t       = load(path_results_SC+list_cultures(num)+"_SC_CI.mat");
    ci_SC   = t.ci; clear t;
    t       = load(path_results_XCov+list_cultures(num)+"_XCov_CI.mat");
    ci_XCov = t.ci; clear t;
    
    
    %=============       initializing count variables        =============%
    
    count_TE     = zeros(size(TE));
    count_SC     = zeros(size(SC));
    count_XCov   = zeros(size(XCov));
    count_ciTE   = zeros(size(ci_TE));
    count_ciSC   = zeros(size(ci_SC));
    count_ciXCov = zeros(size(ci_XCov));
    
    
    %=============      initializing z-scoring counters      =============%
    
    sum_TE     = zeros(size(TE));      sum_TE_sq     = zeros(size(TE));
    sum_SC     = zeros(size(SC));      sum_SC_sq     = zeros(size(SC));
    sum_XCov   = zeros(size(XCov));    sum_XCov_sq   = zeros(size(XCov));
    sum_ciTE   = zeros(size(ci_TE));   sum_ciTE_sq   = zeros(size(ci_TE));
    sum_ciSC   = zeros(size(ci_SC));   sum_ciSC_sq   = zeros(size(ci_SC));
    sum_ciXCov = zeros(size(ci_XCov)); sum_ciXCov_sq = zeros(size(ci_XCov));
    
    %============================== Null Model ==============================%
    
    for m=1:n_permutations
        m
                
        %============= Jittering spike-trains with Norm(spike-time,time_jitter) =============%
        
        J = cell(size(data.spikes));

        for i = 1:numel(data.spikes)
            rand_vector = normrnd(0, time_jitter * binsize, size(data.spikes{i}));
            J{i} = ceil(data.spikes{i} + rand_vector);
            J{i} = max(J{i}, 1);
            J{i} = unique(J{i});
        end
        
        %============= binarizing jittered spike-trains =============%
        
        nonempty_spikes = data.spikes(~cellfun('isempty', data.spikes));
        if isempty(nonempty_spikes)
            tmax = 0;
        else
            tmax = max(cellfun(@max, nonempty_spikes));
        end
        tmax = tmax + time_jitter / binsize;
        binarized = false(data.nNeurons,ceil(tmax/binsize)+1);
        
        for i=1:length(J)
            a=J{i};
            spiking_bins = floor(a/binsize)+1;
            binarized(i,spiking_bins)=true;
        end
        
        % jittered spike-trains of all neurons (seen as pre-synaptic)
        asdf_pre = SparseToASDF(binarized, 1);     %== turns into cell array
        
        clear binarized J;

        
        
        %============= computing EC for asdf_pre –> asdf_post =============%
        
        if size(asdf_pre)==size(asdf_post)
            %ASDFTE_CN_perm(timeseries_jit_pre, timeseries_post, x_delay(s), x_order, y_order, CI_tau_window);
            if maxdelay==1
                [peakTE] = ASDFTE_CN_perm(asdf_pre, asdf_post, 1,1,1)
            else
                [peakTE, ciTE, peakSC, ciSC, peakXCov, ciXCov] = ASDFTE_CN_perm(asdf_pre, asdf_post, 0:maxdelay,1,1,CI_tau_window);
            end
            
        else
            disp("problem!!");
            disp(size(asdf_pre),' ',size(asdf_pre));
        end
        
        
        %=============     updating count variables     =============%

        count_TE      = count_TE      + (peakTE   >= TE);
        count_SC      = count_SC      + (peakSC   >= SC);
        count_XCov    = count_XCov    + (peakXCov >= XCov);
        count_ciTE    = count_ciTE    + (ciTE     >= ci_TE);
        count_ciSC    = count_ciSC    + (ciSC     >= ci_SC);
        count_ciXCov  = count_ciXCov  + (ciXCov   >= ci_XCov);

        %=============    updating z-scoring counters   =============%

        sum_TE        = sum_TE        + peakTE;
        sum_TE_sq     = sum_TE_sq     + peakTE.^2;
        
        sum_SC        = sum_SC        + peakSC;
        sum_SC_sq     = sum_SC_sq     + peakSC.^2;

        sum_XCov      = sum_XCov      + peakXCov;
        sum_XCov_sq   = sum_XCov_sq   + peakXCov.^2;

        sum_ciTE      = sum_ciTE      + ciTE;
        sum_ciTE_sq   = sum_ciTE_sq   + ciTE.^2;

        sum_ciSC      = sum_ciSC      + ciSC;
        sum_ciSC_sq   = sum_ciSC_sq   + ciSC.^2;

        sum_ciXCov    = sum_ciXCov    + ciXCov;
        sum_ciXCov_sq = sum_ciXCov_sq + ciXCov.^2;


        clear peakTE ciTE peakSC ciSC peakXCov ciXCov;
        
    end


    %=============   dediagonalizing count matrices   =============%
     
    count_TE     = count_TE     - diag(diag(count_TE));
    count_SC     = count_SC     - diag(diag(count_SC));
    count_XCov   = count_XCov   - diag(diag(count_XCov));
    count_ciTE   = count_ciTE   - diag(diag(count_ciTE));
    count_ciSC   = count_ciSC   - diag(diag(count_ciSC));
    count_ciXCov = count_ciXCov - diag(diag(count_ciXCov));

    %=======     Saving p-value matrices     =======%
        
    P_TE     = count_TE     / n_permutations;
    P_SC     = count_SC     / n_permutations;
    P_XCov   = count_XCov   / n_permutations;
    P_ciTE   = count_ciTE   / n_permutations;
    P_ciSC   = count_ciSC   / n_permutations;
    P_ciXCov = count_ciXCov / n_permutations;
    
    %=======      Calculate means over permutations      =======%
    
    mean_TE     = sum_TE     / n_permutations;
    mean_SC     = sum_SC     / n_permutations;
    mean_XCov   = sum_XCov   / n_permutations;
    mean_ciTE   = sum_ciTE   / n_permutations;
    mean_ciSC   = sum_ciSC   / n_permutations;
    mean_ciXCov = sum_ciXCov / n_permutations;
    
    %=======       Calculate std devs over permutations       =======%
    
    % std = sqrt((sum_sq - (sum^2 / N)) / (N - 1))
    std_TE     = sqrt(abs(sum_TE_sq     - (sum_TE.^2)     / n_permutations) / (n_permutations - 1));
    std_SC     = sqrt(abs(sum_SC_sq     - (sum_SC.^2)     / n_permutations) / (n_permutations - 1));
    std_XCov   = sqrt(abs(sum_XCov_sq   - (sum_XCov.^2)   / n_permutations) / (n_permutations - 1));
    std_ciTE   = sqrt(abs(sum_ciTE_sq   - (sum_ciTE.^2)   / n_permutations) / (n_permutations - 1));
    std_ciSC   = sqrt(abs(sum_ciSC_sq   - (sum_ciSC.^2)   / n_permutations) / (n_permutations - 1));
    std_ciXCov = sqrt(abs(sum_ciXCov_sq - (sum_ciXCov.^2) / n_permutations) / (n_permutations - 1));
    
    %=======            Z-scoring             =======%
    
    Zscored_TE     = (TE      - mean_TE)     ./ std_TE;
    Zscored_SC     = (SC      - mean_SC)     ./ std_SC;
    Zscored_XCov   = (XCov    - mean_XCov)   ./ std_XCov;

    Zscored_ciTE   = (ci_TE   - mean_ciTE)   ./ std_ciTE;
    Zscored_ciSC   = (ci_SC   - mean_ciSC)   ./ std_ciSC;
    Zscored_ciXCov = (ci_XCov - mean_ciXCov) ./ std_ciXCov;



    %=======       Saving sum matrices       =======%
    
    sum = sum_TE;
    filename = fullfile(path_results_TE, [list_cultures{num}, '_TEpk_sum.mat']);
    save(filename, 'sum'); clear sum;

    sum = sum_SC;
    filename = fullfile(path_results_SC, [list_cultures{num}, '_SCpk_sum.mat']);
    save(filename, 'sum'); clear sum;

    sum = sum_XCov;
    filename = fullfile(path_results_XCov, [list_cultures{num}, '_XCovpk_sum.mat']);
    save(filename, 'sum'); clear sum;

    sum = sum_ciTE;
    filename = fullfile(path_results_TE, [list_cultures{num}, '_TE_CI_sum.mat']);
    save(filename, 'sum'); clear sum;

    sum = sum_ciSC;
    filename = fullfile(path_results_SC, [list_cultures{num}, '_SC_CI_sum.mat']);
    save(filename, 'sum'); clear sum;

    sum = sum_ciXCov;
    filename = fullfile(path_results_XCov, [list_cultures{num}, '_XCov_CI_sum.mat']);
    save(filename, 'sum'); clear sum;
    
    
    %=======       Saving sum_sq matrices       =======%
    
    sum_sq = sum_TE_sq;
    filename = fullfile(path_results_TE, [list_cultures{num}, '_TEpk_sum_sq.mat']);
    save(filename, 'sum_sq'); clear sum_sq;

    sum_sq = sum_SC_sq;
    filename = fullfile(path_results_SC, [list_cultures{num}, '_SCpk_sum_sq.mat']);
    save(filename, 'sum_sq'); clear sum_sq;

    sum_sq = sum_XCov_sq;
    filename = fullfile(path_results_XCov, [list_cultures{num}, '_XCovpk_sum_sq.mat']);
    save(filename, 'sum_sq'); clear sum_sq;

    sum_sq = sum_ciTE_sq;
    filename = fullfile(path_results_TE, [list_cultures{num}, '_TE_CI_sum_sq.mat']);
    save(filename, 'sum_sq'); clear sum_sq;

    sum_sq = sum_ciSC_sq;
    filename = fullfile(path_results_SC, [list_cultures{num}, '_SC_CI_sum_sq.mat']);
    save(filename, 'sum_sq'); clear sum_sq;

    sum_sq = sum_ciXCov_sq;
    filename = fullfile(path_results_XCov, [list_cultures{num}, '_XCov_CI_sum_sq.mat']);
    save(filename, 'sum_sq'); clear sum_sq;
    
    
    %=======       Saving count matrices       =======%
    
    count = count_TE;
    filename = fullfile(path_results_TE, [list_cultures{num}, '_TEpk_Pcount.mat']);
    save(filename, 'count'); clear count;

    count = count_SC;
    filename = fullfile(path_results_SC, [list_cultures{num}, '_SCpk_Pcount.mat']);
    save(filename, 'count'); clear count;

    count = count_XCov;
    filename = fullfile(path_results_XCov, [list_cultures{num}, '_XCovpk_Pcount.mat']);
    save(filename, 'count'); clear count;

    count = count_ciTE;
    filename = fullfile(path_results_TE, [list_cultures{num}, '_TE_CI_Pcount.mat']);
    save(filename, 'count'); clear count;

    count = count_ciSC;
    filename = fullfile(path_results_SC, [list_cultures{num}, '_SC_CI_Pcount.mat']);
    save(filename, 'count'); clear count;

    count = count_ciXCov;
    filename = fullfile(path_results_XCov, [list_cultures{num}, '_XCov_CI_Pcount.mat']);
    save(filename, 'count'); clear count;
    
    
    %=======     Saving p-value matrices     =======%
    
    P = P_TE;
    filename = fullfile(path_results_TE, [list_cultures{num}, '_TEpk_Pval.mat']);
    save(filename, 'P'); clear P;
    
    P = P_SC;
    filename = fullfile(path_results_SC, [list_cultures{num}, '_SCpk_Pval.mat']);
    save(filename, 'P'); clear P;
    
    P = P_XCov;
    filename = fullfile(path_results_XCov, [list_cultures{num}, '_XCovpk_Pval.mat']);
    save(filename, 'P'); clear P;
    
    P = P_ciTE;
    filename = fullfile(path_results_TE, [list_cultures{num}, '_TE_CI_Pval.mat']);
    save(filename, 'P'); clear P;
    
    P = P_ciSC;
    filename = fullfile(path_results_SC, [list_cultures{num}, '_SC_CI_Pval.mat']);
    save(filename, 'P'); clear P;
    
    P = P_ciXCov;   
    filename = fullfile(path_results_XCov, [list_cultures{num}, '_XCov_CI_Pval.mat']);
    save(filename, 'P'); clear P;
    
    
    
    %=======     Saving Z-scored matrices     =======%
    
    Zscored_EC = Zscored_TE;
    filename = fullfile(path_results_TE, [list_cultures{num}, '_TEpk_Zscored.mat']);
    save(filename, 'Zscored_EC'); clear Zscored_EC;
    
    Zscored_EC = Zscored_SC;
    filename = fullfile(path_results_SC, [list_cultures{num}, '_SCpk_Zscored.mat']);
    save(filename, 'Zscored_EC'); clear Zscored_EC;
    
    Zscored_EC = Zscored_XCov;
    filename = fullfile(path_results_XCov, [list_cultures{num}, '_XCovpk_Zscored.mat']);
    save(filename, 'Zscored_EC'); clear Zscored_EC;
    
    Zscored_EC = Zscored_ciTE;
    filename = fullfile(path_results_TE, [list_cultures{num}, '_TE_CI_Zscored.mat']);
    save(filename, 'Zscored_EC'); clear Zscored_EC;
    
    Zscored_EC = Zscored_ciSC;
    filename = fullfile(path_results_SC, [list_cultures{num}, '_SC_CI_Zscored.mat']);
    save(filename, 'Zscored_EC'); clear Zscored_EC;
    
    Zscored_EC = Zscored_ciXCov;   
    filename = fullfile(path_results_XCov, [list_cultures{num}, '_XCov_CI_Zscored.mat']);
    save(filename, 'Zscored_EC'); clear Zscored_EC;
    
end
