function TECN_Sign(cultureID,bs,delay_ms,CI_tau_window,time_jitter,n_permutations,path,path_out)

    args.cultureID     = char(cultureID);                 % Culture folder
    args.bs            = str2double(bs);                  % chosen bin size for data (spike times binning)
    args.delay_ms      = str2double(delay_ms);            % max delay for TE-SC-XCov computation
    args.CI_tau_window = str2double(CI_tau_window);       % CI time-window length

    args.n_permutations = str2double(n_permutations);     % number of null model generations for significant test
    args.time_jitter    = str2double(time_jitter);        % (ms) time sigma for jittering 
                                                          % the sender neuron spike-times with N(spike_time,sigma)
    args.path           = char(path);
    args.path_out      = char(path_out);

    
    %=== Add path to TE-SC-XCov package ===%
    addpath("../../SC_TE_package/");


    % . . . Log directory
    logs_path = fullfile('./', '_logs_jobs2');
    mkdir(logs_path);

    % . . . Log file
    log_file = fullfile(logs_path, [args.cultureID, '_output.log']);
    diary(log_file);
    diary('on');
    fprintf('Start time: %s\n', datestr(datetime('now')));


    %==============================%==============================%==============================
    %                                        P A T H S
    %==============================%==============================%==============================


    path_data = [args.path, args.cultureID,'/'];

    suffix = ['_binsize', num2str(args.bs),'/'];
    
    mkdir(path_out);
    path_results      = fullfile(path_out,args.cultureID,'/')
    path_results_TE   = fullfile(path_results, ['TECN_TE', suffix]);
    path_results_SC   = fullfile(path_results, ['TECN_SC', suffix]);
    path_results_XCov = fullfile(path_results, ['TECN_XCov', suffix]);
    path_results_CN   = fullfile(path_results, ['TECN_CN', suffix]);

    mkdir(path_results);
    mkdir(path_results_TE);
    mkdir(path_results_SC);
    mkdir(path_results_XCov);
    mkdir(path_results_CN);


    %==============================%==============================%==============================
    %                                      F e a t u r e s
    %==============================%==============================%==============================


    time_step    = 0.05;                          % (ms) sampling time-step
    bin_size     = fix(args.bs / time_step);      % steps per bin
    max_delay_TE = fix(args.delay_ms / args.bs);  % max delay in bins

    disp(['Bin size (steps): ', num2str(bin_size)]);


    %==============================%==============================%==============================
    %                         D E L A Y E D   T E ,   S C ,   X C o v 
    %                              S I G N I F I C A N T   T E S T
    %==============================%==============================%==============================

    %===    List of cultures we want to analyze
    %  (modify the script args and the list below to compute TE,SC,XCov for more than one culture)
    
    cultures_DIV = ["Cult"];

    %===    TE significant test
    NullModel_SignTest(path_data, path_results_TE, path_results_SC, path_results_XCov, cultures_DIV, args.time_jitter, bin_size, args.n_permutations, max_delay_TE, args.CI_tau_window)
    
    disp('.....done');

%================================================================================================


fprintf('End time: %s\n', datestr(datetime('now')));
diary('off');

