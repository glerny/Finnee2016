%% DESCRIPTION
%signFig a
%   PARAMETERS:
%
%   OUTPUT
%
%
%% Copyright 
% BSD 3-Clause License
% Copyright 2016-2017 G. Erny (guillaume@fe.up.pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function dataOut = trailRem(dataIn , n)

narginchk(2, 2)
%% FUNCTION CORE
 provMat      = [dataIn(2:end, n); 0];
 provMat(:,2) = dataIn(:, n);
 provMat(:,3) = [0; dataIn(1:end-1, n)];
 dataOut      = dataIn(sum(provMat, 2) > 0, :);