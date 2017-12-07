%% DESCRIPTION 
% QUERYMZML is used to retrieve information from the mzML dump file. 
%  
%% INPUT PARAMETERS
% *Compulsory 
%   _obj_    : The Finnee object 
%   _string_ : The string to be looked for. It can either be 'all' to
%       print the full mzML part, or a key word. 
%  
%% EXAMPLES
% * myFinnee = myFinnee.doCentroid(3, 'LocalMax:2:100'); 
%       Will calculated, for each MS profile scan, their correcponsing MS
%       centroid scan using the LocalMax algorithm. 
%  
%% COPYRIGHT
% Copyright BSD 3-Clause License Copyright 2016-2017 G. Erny
% (guillaume@fe.up.pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function queryMZMLDump(obj, string)

fprintf('\n')
if strcmpi(string, 'all')
    for ii = 1:length(obj.MZMLDump)
        fprintf('\n%s', obj.MZMLDump{ii});
    end
else
    
    k  = strfind(obj.MZMLDump, string);
    Id = find(~cellfun(@isempty,k));
    if ~isempty(Id)
        for ii = 1:length(Id)
            fprintf('\n%s', obj.MZMLDump{Id(ii)});
        end
    else
        fprintf('\n %s was not found in MZMLDump', string)
    end
end
fprintf('\n')
