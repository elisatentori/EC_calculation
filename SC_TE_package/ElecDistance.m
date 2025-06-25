% [r] = ElecDistance(path_data, b_elecs, h_elecs)
%
% Parameters:
%   path_data  - path to which save the matrix
%   b_elecs    - number of electrodes of the basis
%                assuming a rectangular-based MEA
%   h_elecs    - number of electrodes of the heigh
%
% Returns:
%   array with all possible distances in a MEA
%   (saves the list of distances in a file at path_data)

%==============================================================================%
% Copyright (c) 2022, University of Padua, Italy							   %
% All rights reserved.														   %
%																			   %
% Authors: Elisa Tentori (elisa.tentori@phd.unipd.it)						   %
%          LiPh Lab - NeuroChip Lab, University of Padua, Italy				   %
%==============================================================================%

function [r]=ElecDistance(path_data,b_elecs,h_elecs)

tot_elecs = b_elecs*h_elecs;  % tot electrodes of the grid

c_x_ = zeros(1,tot_elecs);    % x_coordinate of electrode
c_y_ = zeros(1,tot_elecs);    % y_coordinate of electrode

idx = 1;
for i=1:h_elecs
    for j=1:b_elecs
        c_x_(idx)=i;
        c_y_(idx)=j;
        idx=idx+1;
    end
end

% distance between the electrode in one corner of the grid and all the others + itself
d = zeros(1,tot_elecs);       % array of distances
for i=1:length(c_x_)
    d(i) = sqrt((c_x_(1)-c_x_(i))^2+(c_y_(1)-c_y_(i))^2);
end

%saving the list of possible distances
r = sort(unique(d));

filename=path_data+"r_8x8.txt";
dlmwrite(filename,r);
