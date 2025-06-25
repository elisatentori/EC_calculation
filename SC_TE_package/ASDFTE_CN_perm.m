% [te_result, ci_te_result, sc_result, ci_sc_result, xcov_result, ci_xcov_result] = ASDFTE_CN_perm(asdf, asdf2, j_delay, i_order, j_order, windowsize)
%
% Compute transfer entropy (TE), signed cross-correlation (SC), and cross-covariance (XCov)
% between original and jittered spike-time series via permutation control.
%
% Parameters:
%
%   asdf        - ASDF cell array of original spike-time data
%   asdf2       - ASDF cell array of jittered spike-time data
%   j_delay     - Vector of sender delays to test (default = 1)
%   i_order     - Markov order of receiver (default = 1)
%   j_order     - Markov order of sender   (default = 1)
%   windowsize  - Odd integer window for coincidence index (default = 5)
%
% Returns:
%
%   te_result      - (nNeu×nNeu) peak transfer entropy values for each i→j
%   ci_te_result   - (nNeu×nNeu) coincidence index for TE (empty for single-delay)
%   sc_result      - (nNeu×nNeu) peak signed cross-correlation values for each i→j
%   ci_sc_result   - (nNeu×nNeu) coincidence index for SC (empty for single-delay)
%   xcov_result    - (nNeu×nNeu) peak cross-covariance values for each i→j
%   ci_xcov_result - (nNeu×nNeu) coincidence index for XCov (empty for single-delay)
%
% Examples:
% >> [te, ci_te, sc, ci_sc, xcov, ci_xcov] = ASDFTE_CN_perm(asdf, asdf2);
% >> [te, ci_te, sc, ci_sc, xcov, ci_xcov] = ASDFTE_CN_perm(asdf, asdf2, 1:30);
% >> [te, ci_te, sc, ci_sc, xcov, ci_xcov] = ASDFTE_CN_perm(asdf, asdf2, 0:10, 2, 3, 7);
%



%===============================================================================%
% Code update                                                                   %
% Copyright (c) 2025, University of Padua, Italy                                %
% All rights reserved.                                                          %
%                                                                               %
% Authors: Elisa Tentori (elisa.tentori@phd.unipd.it)                           %
%          LiPh Lab - NeuroChip Lab, University of Padua, Italy                 %
%                                                                               %
% Function written to calculate TE,SC and SCov between pre-synaptic             %
% jittered timeseries and post-synaptic not-jittered timeseries.                %
%===============================================================================%


%===============================================================================%
% Original toolbox:                                                             %
% Copyright (c) 2011, The Trustees of Indiana University                        %
% All rights reserved.                                                          %
%                                                                               %
% Authors: Michael Hansen (mihansen@indiana.edu), Shinya Ito (itos@indiana.edu) %
%===============================================================================%


function [te_result, ci_te_result, sc_result, ci_sc_result, xcov_result, ci_xcov_result] = ASDFTE_CN_perm(asdf, asdf2, j_delay, i_order, j_order, windowsize)


% ========= Set default arguments ========= %

    % default args
    if nargin<3, j_delay    = 1; end
    if nargin<4, i_order    = 1; end
    if nargin<5, j_order    = 1; end
    if nargin<6, windowsize = 5; end

    % extract sizes
    info        = asdf{end};
    num_neurons = info(1);
    n_bins      = info(2);
    num_delays  = length(j_delay);

   % Single-delay branch: TE, SC and XCov for one delay
   
    if num_delays == 1
        [te_result, coincidence] = transent_CN_perm(asdf, asdf2, j_delay, i_order, j_order);
        ci_te_result = [];                    % no CI with single delay
        
        sc_result    = coincidence;           % SC = raw coincidences
        ci_sc_result = [];                    % no CI with single delay

        % compute XCov for single delay
        spike_prob     = diag(coincidence) / n_bins;
        prob_product   = spike_prob * spike_prob';
        xcov_mat       = coincidence / n_bins - prob_product;
        xcov_result    = xcov_mat;
        ci_xcov_result = [];                  % no CI with single delay
    
    else

        % multi‐delay: preallocate storage
        all_te     = zeros(num_neurons, num_neurons, num_delays);
        all_coincs = zeros(num_neurons, num_neurons, num_delays);
        for d = 1:num_delays
            % compute TE and coincidences at each delay
            [all_te(:,:,d), all_coincs(:,:,d)] = transent_CN_perm(asdf, asdf2, j_delay(d), i_order, j_order);
        end

        % TE: optionally drop zero‐lag then find peak & CI
        if j_delay(1)==0
            te_data = all_te(:,:,2:end);
        else
            te_data = all_te;
        end
        [te_result, ~]    = max(te_data, [], 3);   % peak TE
        ci_te_result      = zeros(num_neurons);
        for i = 1:num_neurons
            for j = 1:num_neurons
                ci_te_result(i,j) = CIReduce(te_data(i,j,:), windowsize);
            end
        end

        % signed cross‐correlation: center around mean, then peak & CI
        flipped = flip(all_coincs(:,:,2:end), 3);
        means   = mean(cat(3, flipped, all_coincs(:,:,1:end-1)), 3);
        sc_data = all_coincs(:,:,2:end) - means;
        sc_result    = zeros(num_neurons);
        ci_sc_result = zeros(num_neurons);
        for i = 1:num_neurons
            for j = 1:num_neurons
                tmp = squeeze(sc_data(i,j,:));
                [~, idx]        = max(abs(tmp));
                sc_result(i,j)  = tmp(idx);
                ci_sc_result(i,j) = CIReduce(abs(tmp), windowsize);
            end
        end

        % cross‐covariance: normalize, then peak & CI
        spike_prob     = diag(all_coincs(:,:,1)) / n_bins;
        prob_product   = spike_prob * spike_prob';
        xcov_data      = all_coincs(:,:,2:end) / n_bins - prob_product;
        xcov_result    = zeros(num_neurons);
        ci_xcov_result = zeros(num_neurons);
        for i = 1:num_neurons
            for j = 1:num_neurons
                tmp = squeeze(xcov_data(i,j,:));
                [~, idx]          = max(abs(tmp));
                xcov_result(i,j)  = tmp(idx);
                ci_xcov_result(i,j) = CIReduce(abs(tmp), windowsize);
            end
        end
    end
end



