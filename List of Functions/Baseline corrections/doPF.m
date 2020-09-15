%% Description
% Polynomial fitting
% Unpublished results
%
%% Copyright 
% BSD 3-Clause License
% Copyright 2016-2017 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [z, w] = doPF(XY, n, w)
% input: XY 2xn array,
%        n order iof the polynomial
% output:
%        z basline model
%        bslPts false if peak pts, true in baseline pts

N = length(XY(:,1));
if nargin == 2
    w = true(N, 1);
end
L = 0;
iterMax = 20;

for ii = 1:iterMax
    
    p = polyfit(XY(w, 1), XY(w, 2), n);
    z = polyval(p, XY(:,1));
    d = XY(:,2) - z;
    
    dn = d(w);
    % dn is the part of the data that are baseline points
    
    m = mean(dn);
    s = std(dn);
    w = abs((m-d)/s) < 3;
    %%% ORIGINAL
    
    if sum(w) ==  L, break; end
    if sum(w) <= 0.5*N, break; end

    
    L = sum(w);
    % Check the baseline points and stop the loop if no more points are
    % detected as baseline
    
end

