function Rs = GetResolution2(MeanX,StdX)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

for ii = 1:length(MeanX)
    Rs(:, ii) = abs((MeanX(ii)-MeanX)./(2*(StdX(ii)+StdX)));
end

end

