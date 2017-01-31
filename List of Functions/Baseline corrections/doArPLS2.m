%% Description
% DOARPLS2 modified arPLS algorithm. Function published in:
%
%% Copyright 
% Copyright 2016-2017 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [z, bslPts] = doArPLS2(y, lambda)
% input: 
%       y y values, 
% output:
%       z basline model
%       bslPts false if peak pts, true in baseline pts

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