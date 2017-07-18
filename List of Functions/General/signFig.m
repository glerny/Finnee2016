%% DESCRIPTION
%signFig a
%   PARAMETERS:
%
%   OUTPUT
%
%
%% Copyright 
% BSD 3-Clause License
% Copyright 2016-2017 G. Erny (guillaume@fe.up.pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function n = signFig(number)

narginchk(1, 1)
%% FUNCTION CORE
string =  sprintf('%e',number);
Id     = strfind(string, 'e');
n      = str2double(string(Id+1:end));
if n > 0
    n = 0;
else
    n = - n;
end