function [mzAxis, r2, data4axis] = extrapolMZ(obj, n, Lim)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin == 1
    Lim = [0 inf];
    n   = 3;
elseif nargin == 2
    Lim = [0 inf];
end

MS            = obj.Data;
W2C           =  MS(:,2);
W2C(end+1)    = 0;
W2C(2:end, 2) = MS(:,2);
Id2c          = [0; find(sum(W2C, 2) == 0); size(MS, 1)+1];
if size(Id2c,1) > 2
    data4axis = [];
    for ii = 1:length(Id2c)-1
        cut = MS(Id2c(ii)+1:Id2c(ii+1)-1, 1);
        if size(cut, 1) > 1
            dt2add = [cut(1:end-1), cut(2:end) - cut(1:end-1)];
            data4axis = [data4axis; dt2add]; %#ok<AGROW>
        end
    end
end

[p, S, mu] = polyfit(data4axis(:,1), data4axis(:,2), n);
R = corrcoef(polyval(p, data4axis(:,1), S, mu), data4axis(:,2));
r2 = R(1,2);
if Lim(1) == 0
    mzAxis = MS(1,1);
else
    mzAxis = Lim(1);
end

if Lim(2) == inf
    Lim(2) = MS(end,1);
end

while mzAxis(end) < Lim(2)
    int = polyval(p, mzAxis(end), S, mu);
    mzAxis(end+1) = mzAxis(end) + int; %#ok<AGROW>
end
mzAxis = mzAxis';

end

