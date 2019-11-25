function Ds = GetDistance(X, val)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

switch val
    case 1
        for ii = 1:length(X)
            Ds(:, ii) = sqrt((X(ii)-X).^2);
        end
    case 2
        for ii = 1:length(X)
            Ds(:, ii) = sqrt((X(ii)-X).^2)/X(ii);
        end
        
        
end