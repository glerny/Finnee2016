%% DESCRIPTION
%
%% Copyright 
% BSD 3-Clause License
% Copyright 2016-2017 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function [Log, partial] = decipherLog(strLog)

Log = {};
partial = strsplit(strLog, '|');
for ii = 1:size(partial, 2)
    [token, remain] = strtok(strsplit(partial{ii}, ' '), '=');
    Log{ii}.dtsType = token{1};
    Log{ii}.dtsId   = str2double(remain{1}(2:end)); %#ok<*AGROW>
end
end

