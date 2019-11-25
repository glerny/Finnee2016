function [out, lmax, Nz, obj] = findPeaks(obj, Wdx, spkeSz)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

XY = [obj.AxisTm.Data, (sum(obj.Smoothed, 1))'];
XY(:,3) = XY(:,2);

if nnz(XY(:,2)) >= 0.9*length(XY(:,2))
    obj = obj.doBslCorrection;
    Nz = 0;
else
    Nz = 0;
end

% rem spikes
list = find(XY(:,2) == 0);
list(:,2) = [1; diff(find(XY(:,2) == 0))];
ix = find(list(:,2) > 1 & list(:,2) <= spkeSz);
for ii = 1:length(ix)
    XY(list(ix(ii)-1,1):list(ix(ii),1), 3) = 0;
end

[redMS, lmax]  = LocalMaxima(XY(:, [1 3]), Wdx, Nz);
 ctr = [];
 
for ii = 1:length(lmax)
    sumA = trapz(obj.AxisMZ.Data, obj.Smoothed(:, lmax(ii)));
    ctr(ii) = trapz(obj.AxisMZ.Data, ...
        obj.Smoothed(:, lmax(ii)).*obj.AxisMZ.Data)/sumA;
end
out = [ctr',redMS]; 


