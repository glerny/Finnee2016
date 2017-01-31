%% DESCRIPTION
%
%% Copyright 
% BSD 3-Clause License
% Copyright 2016-2017 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


function [bool, K, field] = fieldfind( structIn, label)
bool = 0; K = nan; field = {};
for ii = 1:length(structIn)
    if strcmp(structIn{ii}.label, label)
        bool = 1;
        break
    end
end
if bool
    K = ii;
    field = structIn{ii}.field;
end
end

