function obj = doBslCorrection(obj)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% step 1

Area = obj.StoredData;

[p, q] = size(Area);
X = obj.AxisTm.Data;

for ii = 1:p
     [z(:,ii), w] = doPF([X, Area(ii, :)'], 1);
end

obj.StoredData = round(obj.StoredData-z');


end





