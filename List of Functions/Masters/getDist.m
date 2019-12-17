function [Dist, Std] = getDist(Clusters)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

for ii = 1:size(Clusters)
    Dist(ii, :) = abs(Clusters(:, 6) - Clusters(ii, 6));
    Std(ii, :) = sqrt(Clusters(:, 7).^2 + Clusters(ii, 7)^2);
end
end

