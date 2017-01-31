%% DESCRIPTION
%
%% Copyright 
% BSD 3-Clause License
% Copyright 2016-2017 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 



function A = spread( yy, wdz)

A = zeros(length(yy), 2*wdz + 1);
for ii = 1:2*wdz+1;
    ind1 = max(1, wdz+2-ii);
    ind2 = min(length(yy), length(yy) + wdz+1-ii);
    A(length(yy)-ind2+1:length(yy)-ind1+1,ii) = yy(ind1:ind2);
end


end

