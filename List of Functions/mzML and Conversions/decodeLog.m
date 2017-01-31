%% DESCRIPTION
%
%% Copyright 
% BSD 3-Clause License
% Copyright 2016-2017 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


function strLog = decodeLog(Log)

A = strsplit(Log, ' ');
strLog.definition = A{1};

if length(A) > 1
    for ii = 2:length(A)
        B = strsplit(A{ii}, '=');
        strLog.(B{1}) = B{2};
    end
end
    

end

