% [te_result, ci_te_result, te_delay, all_delayed_te, sc_result, ci_sc_result, sc_delay, centered_all_coincs, xcov_result, ci_xcov_result, xcov_delay, XCov_d, all_coincs] = ASDFTE_CN(asdf, j_delay, i_order, j_order, windowsize)
%
%
% Parameters:
%
%   asdf             – Time series in Another Spike Data Format (ASDF)
%   j_delay          – Number of bins to lag sender (j series) or a vector [default 1]
%   i_order          – Order of receiver [default 1]
%   j_order          – Order of sender [default 1]
%   windowsize       – window size used for Coincidence Index calculation (odd number only)
%
%
% Returns:
%
%
%   te_result        – (nNeurons×nNeurons) matrix of peak transfer entropy (TE) values for each i→j.
%   ci_te_result     – (nNeurons×nNeurons) coincidence‐index on all_delayed_te for TE (empty for single‐delay).
%   te_delay         – (nNeurons×nNeurons) delay values that maximize TE for each i→j.
%   all_delayed_te   – (nNeurons×nNeurons×D) TE at each tested delay (D = num_delays, or num_delays–1 if j_delay(1)==0).
%
%   sc_result        – (nNeurons×nNeurons) signed cross‐correlation peak values for each i→j.
%   ci_sc_result     – (nNeurons×nNeurons) coincidence‐index on centered_all_coincs for SC (empty for single‐delay).
%   sc_delay         – (nNeurons×nNeurons) delay values that maximize SC for each i→j.
%   centered_all_coincs – (nNeurons×nNeurons×(D-1)) SC values at each delay, mean‐centered (zero‐lag excluded).
%
%   xcov_result      – (nNeurons×nNeurons) cross‐covariance peak values for each i→j.
%   ci_xcov_result   – (nNeurons×nNeurons) coincidence‐index on XCov_d for XCov (empty for single‐delay).
%   xcov_delay       – (nNeurons×nNeurons) delay values that maximize XCov for each i→j.
%   XCov_d           – (nNeurons×nNeurons×(D-1)) normalized XCov at each delay (zero‐lag excluded).
%
%   all_coincs       – (nNeurons×nNeurons×D) raw coincidence counts for each delay.
%
%
% Examples:
%
%   % Single‐delay default (j_delay = 1): returns TE, SC and XCov for delay=1
%   [te,~,~,~,sc,~,~,~,xcov,~,~,~,all_coincs] = ASDFTE_CN(asdf);
%
%   % Explicit single‐delay call: returns TE, SC and XCov for j_delay=5
%   [te,ci_te,te_delay,all_te,sc,ci_sc,sc_delay,sc_all,xcov,ci_xcov,xcov_delay,xcov_all,all_coincs] = ...
%       ASDFTE_CN(asdf, 5);
%
%   % Multi‐delay TE only (ignore SC/XCov outputs if you want):
%   [te,ci_te,te_delay,all_te] = ASDFTE_CN(asdf, 1:10);
%
%   % Full multi‐delay with SC and XCov, custom orders and CI window:
%   [ te,ci_te,te_delay,all_te, sc,ci_sc,sc_delay,sc_all, xcov,ci_xcov,xcov_delay,xcov_all, all_coincs ] = ...
%       ASDFTE_CN(asdf, 0:10, 2, 3, 7);
%


%===============================================================================%
% Code update                                                                   %
% Copyright (c) 2025, University of Padua, Italy                                %
% All rights reserved.                                                          %
%                                                                               %
% Authors: Elisa Tentori (elisa.tentori@phd.unipd.it)                           %
%          LiPh Lab - NeuroChip Lab, University of Padua, Italy                 %
%                                                                               %
% Script modified in order to calculate TE, SC and XCov                         %
%===============================================================================%

%===============================================================================%
% Original Script: ASDFTE.m                                                     %
% Copyright (c) 2011, The Trustees of Indiana University                        %
% All rights reserved.                                                          %
%                                                                               %
% Authors: Michael Hansen (mihansen@indiana.edu), Shinya Ito (itos@indiana.edu) %
%===============================================================================%


function [te_result, ci_te_result, te_delay, all_delayed_te, sc_result, ci_sc_result, sc_delay, centered_all_coincs, xcov_result, ci_xcov_result, xcov_delay, XCov_d, all_coincs] = ASDFTE_CN(asdf, j_delay, i_order, j_order, windowsize)


    % ========= Set default arguments ========= %

    if nargin<2, j_delay = 1; end
    if nargin<3, i_order = 1; end
    if nargin<4, j_order = 1; end
    if nargin<5, windowsize = 5; end


    % =========     info     ========= %

    num_delays  = length(j_delay); % how many different time-lags
    info        = asdf{end};       % last 2 cells of cell array contain info about the time-series
    num_neurons = info(1);         % info(1) contains the dim of the matrix (number of different time-series)
    n_bins      = info(2);         % info(2) contains the number of bins


    % Single-delay branch: TE, SC and XCov for one delay
    if num_delays == 1
        [te_result, coincidence] = transent_CN(asdf, j_delay, i_order, j_order); % Single delay
        te_delay         = ones(num_neurons) * j_delay;
        all_delayed_te   = te_result;
        ci_te_result     = [];

        centered_all_coincs = coincidence;
        sc_delay         = te_delay;
        sc_result        = coincidence;
        ci_sc_result     = [];

        spike_prob     = diag(coincidence) / n_bins;
        prob_product   = spike_prob * spike_prob';
        XCov_d         = coincidence / n_bins - prob_product;
        xcov_delay     = te_delay;
        xcov_result    = XCov_d;
        ci_xcov_result = [];

        all_coincs = coincidence; % for consistency with multi-delay output
        return;


    % ========= delayed TE, SC, XCov ========= %

    else

        % Multi-delay: preallocate
        all_te     = zeros(num_neurons, num_neurons, num_delays);
        all_coincs = zeros(num_neurons, num_neurons, num_delays);
        for d = 1:num_delays % Change this for to parfor for parallelization.
            [all_te(:,:,d), all_coincs(:,:,d)] = transent_CN(asdf, j_delay(d), i_order, j_order);
        end

        % TE reduction (exclude zero-lag if first delay is 0)
        if j_delay(1) == 0
            data_te = all_te(:,:,2:end);
        else
            data_te = all_te;
        end
        [te_result, te_delay] = max(data_te, [], 3);
        all_delayed_te        = data_te;

        % CI for TE
        ci_te_result = zeros(num_neurons);
        for i = 1:num_neurons
            for j = 1:num_neurons
                ci_te_result(i,j) = CIReduce(all_delayed_te(i,j,:), windowsize);
            end
        end


        % Signed cross-correlation (SC)
        flipped             = flip(all_coincs(:,:,2:end), 3);
        means               = mean(cat(3, flipped, all_coincs(:,:,1:end-1)), 3);
        centered_all_coincs = all_coincs(:,:,2:end) - means;
        [~, sc_delay]      = max(abs(centered_all_coincs), [], 3);
        sc_result          = zeros(num_neurons);
        ci_sc_result       = zeros(num_neurons);
        for i = 1:num_neurons
            for j = 1:num_neurons
                idx                = sc_delay(i,j);
                sc_result(i,j)    = centered_all_coincs(i,j,idx);
                ci_sc_result(i,j) = CIReduce(abs(centered_all_coincs(i,j,:)), windowsize);
            end
        end


        % Cross-covariance (XCov)
        spike_prob     = diag(all_coincs(:,:,1)) / n_bins;
        prob_product   = spike_prob * spike_prob';
        XCov_d         = all_coincs(:,:,2:end) / n_bins - prob_product;
        [~, xcov_delay] = max(abs(XCov_d), [], 3);
        xcov_result    = zeros(num_neurons);
        ci_xcov_result = zeros(num_neurons);
        for i = 1:num_neurons
            for j = 1:num_neurons
                idx                 = xcov_delay(i,j);
                xcov_result(i,j)    = XCov_d(i,j,idx);
                ci_xcov_result(i,j) = CIReduce(abs(XCov_d(i,j,:)), windowsize);
            end
        end

    end
end


