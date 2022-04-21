%% DESCRIPTION
% 1. INTRODUCTION
%
% 2. PARAMETERS
%
% 3. EXAMPLES
%       
% 4. COPYRIGHT
% Copyright 2021 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal

function [lmx, lmy] = LocalMaxima2D(XY, wdx, wdy, thrI)

cshiftXY = zeros([size(XY) ((2*wdx+1)*(2*wdy+1))+1]);

ct = 0;
for cx = -2:2
    for cy = -2:2
        cshiftXY(:, :, end+1) = circshift(XY, [cx cy]);
    end
end

[lmx, lmy] = find(XY == max(cshiftXY, [], 3, 'omitnan') ...
    & XY >= thrI*max(XY, [],  'all'));