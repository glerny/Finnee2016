function [z, bslPts, parameters] = doArPLS(y, lambda)
% Estimate baseline with arPLS in MATLAB
% Function modified from the work of Baek and co-worker
% Analyst, 2015, 140, 250-257. doi: 10.1039/c4an01061b

N = length(y);
D = diff(speye(N), 2);
H = lambda*(D'*D);
w = ones(N,1);
L = 0;
iterMax = 20;

for ii = 1:iterMax;
    W = spdiags(w, 0, N, N);
    % Cholesky decomposition
    C = chol(W + H);
    z = C \ (C'\(w.*y) );
    d = y - z;
    
    %%% ORIGINAL 
    % dn = d(d<0);
    % replaced by:
    dn = d(w > 0.05);
    % dn is the part of the data that are baseline points
    %%% END MODIFICATION
    
    m = mean(dn);
    s = std(dn);
    w = 1./(1 + exp(2*(abs(d) - (2*s-m))/s));
    %%% ORIGINAL
    % if norm(w - wt)/norm(w) < ratio, break; end
    % w = wt;
    
    % replaced by
    if length(dn) ==  L, break; end
    L = length(dn);
    % Check the baseline points and stop the loop if no more points are
    % detected as baseline
    
end

bslPts = w > 0.05;
parameters(1) = lambda;