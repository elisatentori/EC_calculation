%
%    Save a collection of matrices both .mat (decomment for .txt)
%
%   outdir   : directory in cui salvare i file
%   basename : prefisso comune (es. list_cultures{num})
%   suffixes : cell array di stringhe, es:
%              {'TE_delay1','TE_delay2',...,'SC_delay1',...,'CN_delay0',...}
%   mats     : cell array di matrici, in corrispondenza a suffixes
%   txt_delim: delimiter per writematrix (es. '\t')
%
% Example:
%   suffixes = {'TE_delay1','TE_delay2',...};
%   mats     = {allTE(:,:,1), allTE(:,:,2), ...};
%   save_all_measures(path_results_TE, list_cultures{num}, suffixes, mats, '\t');
%

function save_all_measures(outdir, basename, suffixes, mats, txt_delim)

    if nargin < 5, txt_delim = '\t'; end

    for k = 1:numel(suffixes)
        suf = suffixes{k};
        mat = mats{k};

        fmat = fullfile(outdir, sprintf('%s_%s.mat', basename, suf));
        %ftxt = fullfile(outdir, sprintf('%s_%s.txt', basename, suf));

        % name var based on suffix
        if endsWith(suf, 'Pk')
            peakEC = mat;
            save(fmat, 'peakEC');
        elseif endsWith(suf, '_CI')
            ci = mat;
            save(fmat, 'ci');
        elseif endsWith(suf, 'comm_delays')
            comm_delays = mat;
            save(fmat, 'comm_delays');
        else
            save(fmat, 'mat');
        end

        %writematrix(mat, ftxt, 'Delimiter', txt_delim);
    end
end
