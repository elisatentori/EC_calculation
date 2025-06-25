% [ ] = Calculate_TE_CN(path_data, list_cultures, path_results_TE, path_results_SC, path_results_XCov, binsize, maxdelay, CI_tau_window)
%
% Parameters:
%   path_data           - path to which take the recordings
%   list_cultures       - list of cultures to analyze (name_files without format)
%   path_results_TE     - path to which save results for TE
%   path_results_SC     - path to which save results for SC
%   path_results_XCov   - path to which save results for XCov
%   binsize             - (opt) size of bins to binarize spike times series
%                               size is intended as number of measure time-steps (0.05ms)
%                               [default 20] (20 time-steps=1ms)
%   maxdelay            - (opt) max delay (Unit: number of time-bins) 
%                               of the pre-synaptic neuron at which to calculate TE, SC and XCov
%                               [default 1]
%   CI_tau_window       - (opt) window length for CI (Unit: number of time-bins)
%
% Returns:
%   void funtion
%
% Saved files:
%    ABOUT TE:    saves TE matrix, CI matrix, delay matrix (delay that maximizes TE),
%                       all_TE matrices (i.e. TE matrices each for fixed delay)
%    ABOUT SC:    saves SC matrix, CI matrix, delay matrix (delay that maximizes SC),
%                       all_SC matrices (i.e. SC matrices each for fixed delay)
%    ABOUT XCov:  saves XCov matrix, CI matrix, delay matrix (delay that maximizes XCov),
%                       all_XCov matrices (i.e. XCov matrices each for fixed delay)
%    ABOUT DISTANCES:  saves NxN matrix with distance between each neuron couple


%===============================================================================%
% Copyright (c) 2025, University of Padua, Italy                                %
% All rights reserved.                                                          %
%                                                                               %
% Authors: Elisa Tentori (elisa.tentori@phd.unipd.it)                           %
%          LiPh Lab - NeuroChip Lab, University of Padua, Italy                 %
%===============================================================================%


function Calculate_TE_CN(path_data, list_cultures, path_results_TE, path_results_SC, path_results_XCov, path_results_CN, binsize, maxdelay, CI_tau_window)
% Calculate TE, SC, XCov and save results for each culture

    if nargin<7, binsize = 20; end
    if nargin<8, maxdelay = 1; end
    if nargin<9, CI_tau_window = 5; end

    list_to_load = strcat(path_data, list_cultures, ".mat");
    nCult        = numel(list_cultures);

    for idx = 1:nCult
        data = load(list_to_load(idx));
        fprintf('\n[%d/%d] Processing %s\n', idx, nCult, list_cultures{idx});

        % Binarize spike trains, handling empty entries
        binIdx = cellfun(@(sp) reshape(floor(sp / binsize) + 1, [], 1),  data.spikes, 'UniformOutput', false);

        rows  = repelem((1:data.nNeurons)', cellfun(@numel, binIdx));
        cols  = vertcat(binIdx{:});
        nbins = max([1; cols]);

        binarized = sparse(rows, cols, true, data.nNeurons, nbins);
        fprintf('  Binarized (nNeurons = %d, nbins = %d)\n', data.nNeurons, nbins - 1);


        % Convert to ASDF
        asdf = SparseToASDF(binarized, 1);

        % Compute rate
        rate  = full(sum(binarized,2)'/nbins);
        fname = fullfile(path_results_TE, list_cultures{idx} + "_rate.txt");
        writematrix(rate, fname, 'Delimiter','\t');
        fprintf('  Rate saved: %s\n', fname);

        % Compute EC measures
        fprintf('  Computing EC measures (maxdelay=%d)...\n', maxdelay);
        if maxdelay==1
            [peakTE,CN] = ASDFTE_CN(asdf,1,1,1);
            peakSC      = CN;  SCdelays   = 1;
            peakXCov    = CN;  XCovdelays = 1;
            allTE       = [];  allSC      = [];
            allXCov     = [];  allCN      = CN;
            ciTE        = [];  ciSC       = []; ciXCov = [];
        else
            [peakTE,ciTE,TEdelays,allTE,peakSC,ciSC,SCdelays,allSC,peakXCov,ciXCov,XCovdelays,allXCov,allCN] = ASDFTE_CN(asdf,0:maxdelay,1,1,CI_tau_window);
        end

        % Save EC matrices per delay
        maxd = maxdelay;
        suffixes = {}; mats = {};
        for d = 1:maxd
            suffixes{end+1} = sprintf('TE_delay%d',d);    mats{end+1} = allTE(:,:,d);
            suffixes{end+1} = sprintf('SC_delay%d',d);    mats{end+1} = allSC(:,:,d);
            suffixes{end+1} = sprintf('XCov_delay%d',d);  mats{end+1} = allXCov(:,:,d);
            suffixes{end+1} = sprintf('CN_delay%d',d-1);  mats{end+1} = allCN(:,:,d);
        end
        suffixes{end+1} = sprintf('CN_delay%d',maxd); mats{end+1} = allCN(:,:,maxd+1);
        save_all_measures(path_results_TE, list_cultures{idx}, suffixes(1:2*maxd),   mats(1:2*maxd), '\t');
        save_all_measures(path_results_CN, list_cultures{idx}, suffixes(2*maxd+1:end), mats(2*maxd+1:end), '\t');

        % Save communication delays
        save_all_measures(path_results_TE,   list_cultures{idx}, {'TE_comm_delays'},   {int64(TEdelays)},   '\t');
        save_all_measures(path_results_SC,   list_cultures{idx}, {'SC_comm_delays'},   {int64(SCdelays)},   '\t');
        save_all_measures(path_results_XCov, list_cultures{idx}, {'XCov_comm_delays'}, {int64(XCovdelays)}, '\t');

        % Save peak EC and CI
        save_all_measures(path_results_SC,   list_cultures{idx}, {'SCPk','SC_CI'},   {peakSC, ciSC},   '\t');
        save_all_measures(path_results_XCov, list_cultures{idx}, {'XCovPk','XCov_CI'},{peakXCov, ciXCov},'\t');
        save_all_measures(path_results_TE,   list_cultures{idx}, {'TEPk','TE_CI'},   {peakTE, ciTE},   '\t');

        % Distance matrix
        coords = data.pos(:,1:2);
        mat_d  = squareform(pdist(coords));
        fname  = fullfile(path_results_TE, list_cultures{idx} + "_DistanceMat.mat");
        save(fname, 'mat_d');
        fprintf('  Distance matrix saved: %s\n', fname);
    end

end
        


