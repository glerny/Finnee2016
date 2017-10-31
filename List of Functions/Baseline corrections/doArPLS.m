%% Description
% DOARPLS Estimate baseline with arPLS in MATLAB
% Function published in :
% Baek, S.-J., Park, A., Ahn, Y.-J., Choo, J. (2015) 
% "Baseline correction using asymmetrically reweighted penalized least 
% squares smoothing" Analyst, 140, 250-257. doi: 10.1039/c4an01061b
%
%% Copyright 
% Copyright 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [z, bslPts] = doArPLS(y, lambda, ratio)
% input: y y values, 
%       z basline model
% output:
%       bslPts false if peak pts, true in baseline pts

N = length(y);
D = diff(speye(N), 2);
H = lambda*(D'*D);
w = ones(N,1);
iterMax = 20;

for ii = 1:iterMax;
    W = spdiags(w, 0, N, N);
    % Cholesky decomposition
    C = chol(W + H);
    z = C \ (C'\(w.*y) );
    d = y - z;
    dn = d(d<0);
    m = mean(dn);
    s = std(dn);
    wt = 1./(1 + exp(2*(abs(d) - (2*s-m))/s));
   
    if norm(w - wt)/norm(w) < ratio, break; end
    w = wt;
    
end

bslPts = abs((m-d)/s) < 1.98;
