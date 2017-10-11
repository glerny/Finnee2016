%% DESCRIPTION
% Study will contains multiples peaks lists that described a full study
%
%% LIST OF THE CLASS'S PROPERTIES
%
%% LIST OF THE CLASS'S METHODS
%
%% Copyright
% BSD 3-Clause License
% Copyright 2016-2017 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Study
    
    properties
        Name
        DateOfCreation
        Description
        Replicates
        ListOfTags
        pkl4normal
        Path2Std    % where to save
    end
    
    methods
        function obj = Study(varargin)
            if nargin == 0
                obj.Name = '';
                obj.DateOfCreation = date;
                obj.Description = '';
                obj.Replicates = {};
                obj.ListOfTags = {};
                obj.pkl4normal = [];
                obj.Path2Std   = '';
            end
            
        end
        
        function save(obj)
            %% DESCRIPTION
            if isempty(obj.Path2Std)
                obj.Path2Std = uigetdir(pwd, 'Select the folder of destination');
                if ~ischar(obj.Path2Std)
                    error('myApp:argChk', 'Cancel by user');
                end
            end
            
            myStudy = obj; %#ok<*NASGU>
            save(fullfile(obj.Path2Std, 'myStudy.mat'), 'myStudy')
        end
        
        
    end
end