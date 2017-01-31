%% DESCRIPTION
%FINDCLOSER is used the closest valor to calorIn in the vector (Xx1 array)
%vectorIn
%   PARAMETERS:
%       valorIn
%           is a number
%
%       vectorIn
%           is a Xx1 array of numbers
%   
%   OUTPUT
%       indOut in the indice in vectorIn where the closest valor to valorIn
%       can be found
%
%% Copyright 
% BSD 3-Clause License
% Copyright 2016-2017 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function indOut = findCloser(valorIn, vectorIn)

%% FUNCTION CORE

if valorIn <= vectorIn(1)
    indOut = 1;
elseif valorIn >= vectorIn(end)
    indOut = length(vectorIn);
else
    [~, indOut] = min(abs(vectorIn - valorIn));
end
