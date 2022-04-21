%% Description
% DOARPLS2 modified arPLS algorithm. Function published in:
%
%% Copyright 
% Copyright 2016-2017 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [z, bslPts] = doArPLS2(y, lambda)

N = length(y);
D = diff(speye(N), 2);
H = lambda*(D'*D);
w = ones(N,1);
L = 0;
iterMax = 20;

for ii = 1:iterMax
    W = spdiags(w, 0, N, N);
    % Cholesky decomposition
    C = chol(W + H);
    z = C \ (C'\(w.*y) );
    d = y - z;
    dn = d(w > 0.05);
    
    m = mean(dn);
    s = std(dn);
    w = 1./(1 + exp(2*(abs(d) - (2*s-m))/s));
    if length(dn) ==  L, break; end
    L = length(dn);
end

bslPts = abs((m-d)/s) < 1.98;