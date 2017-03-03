%% Description
% FINNEE is the main class of the Finnee2016 toolbox. It is used to record 
% information about the runs; to maintain the address to the .fin folder 
% and to the .mat object, and will contain all Dataset objects. Finned 
% only work with mzML files
%
%% List of Properties
% MZMLDump: is a structure that contains all the information that were 
% present in the mzML file.
%
% Path2Fin: is a string that is the path to the .fin directory where all 
% files are saved
%
% Options: is an Options4Finnee object were all parameters used by the 
% Finnee constructor method are saved.
%
% FileId: is a string that is a the name of the finnee object
%
% Datasets: is an array of Dataset object. Each dataset contains all data 
% related to the run. Multiple datasets are possible, each dataset would be 
% a different data transformation (smoothing, filtering, centroid data, 
% extracted ion profile).
%
%% List of Methods
% Finnee         : Constructor method. Will load an mzML file (centroid or 
% profile mode) and create from it the Finnee object.
%
% queryMZMl      : not implemented yet
%
% save           : record changes
%
% toMZML         : not implemented yet
%
% doBaselineInPrf: Use only with profile mode scan. Correct every channels 
% (i.e. m/z interval) from baseline drift and background noise. Input 
% parameters are obtained using the Options4bslCorPrf class.
%
% doCentroid     : not implemented yet
%
% doFilterInPrf  : not implemented yet
%
%% Copyright 
% BSD 3-Clause License
% Copyright 2016-2017 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef Finnee

    properties
        MZMLDump % general information from the mzML file
        Path2Fin % address of the .fin folder
        Options  % options used by the constructor method
        FileID   % name 
        Datasets % array of dataset
    end
    
    methods
        function obj = Finnee(varargin)
            % Constructor method, will create the Finnee object.
            
            % Check the options and create the Finnee object
            obj.Options = Options4Finnee(varargin);
            
            switch obj.Options.FileFormat
                case 'mzML'
                    obj = domzML2Finnee(obj);
            end
            
            myFinnee = obj;
            save(fullfile(obj.Path2Fin, 'myFinnee.mat'), 'myFinnee')
            
        end
        
        function queryMZML(obj, what)
            % TOBE DONE
        end
        
        function save(obj)
            myFinnee = obj;
            save(fullfile(obj.Path2Fin, 'myFinnee.mat'), 'myFinnee')
        end
        
        function toMZML(obj, dts)
            % TOBE DONE
        end
        
        function obj = doBaselinePrf(obj, dts, par4bas)
            % Allow to correct a full dataset from baseline drift.
            % Description and more information @ 
            % https://github.com/glerny/Finnee2016/wiki/Baseline-and_noise-correction
            %

            m = length(obj.Datasets) + 1;
            newDts = doBslCorPrf(obj.Datasets{dts}, par4bas, ...
                obj.Path2Fin, obj.Options.MaxFileSize, m);
            newDts.CreatedFrom = obj.Datasets{dts}.Log;
            newDts.Log = ['MSPROFDTS DTS=', num2str(m), ...
                ' BLC=1 ', ...
                par4bas.type, '=', num2str(par4bas.parameter),...
                ' NOI=', num2str(par4bas.noise)];
            
            obj.Datasets{end+1} = newDts;
            
            myFinnee = obj;
            save(fullfile(obj.Path2Fin, 'myFinnee.mat'), 'myFinnee')
        end
        
        function obj = doCentroid(obj, dts, method)
            % doCentroid is used to run the "centroid algorithms" defined 
            % in method. Method shopuld contain the algorithms and
            % associated parameters in the following format: 
            % 'thisMethos:par1:par2'. List of accepted methods:
            %   'LocalMax:par1': 
            %           LocalMax detect every local maxima as any
            %           points that is higher or equal than the par1 
            %           closest neighbours in each side, By default 
            %           par1 = 1
            %
            % examples: 
            %   myFinnee = myFinnee.doCentroid(1, 'LocalMax:2');
            %
            % Description and more information @ 
            % https://github.com/glerny/Finnee2016/wiki/Baseline-and_noise-correction
            
            m = length(obj.Datasets) + 1;
            narginchk(2, 3); 
            if nargin == 2
                method = 'LocalMax:1';
            end
            
            % 1. Checking the entry parameters
            
            if ~strcmp(obj.Datasets{dts}.Format, 'MS profile')
                error('Error. \n Target dataset should be in MS profile format. \n Dataset %i is in %s format',...
                    dts, obj.Datasets{dts}.Format)
            end
            
            mth = strsplit(method, ':');
            switch  lower(mth{1})
                case 'localmax'
                    algo.name = 'localmax';
                    if length(mth) == 1;
                        algo.prm{1} = 1;
                    else
                        algo.prm{1} = str2double(mth{2});
                    end
                otherwise
                    error('Error. \n Centroid algorithm not recognised');
            end
            
            newDts = doCtr(obj.Datasets{dts}, algo, ...
                obj.Path2Fin, obj.Options.MaxFileSize, m);
            newDts.CreatedFrom = obj.Datasets{dts}.Log;
            obj.Datasets{end+1} = newDts;
            
            myFinnee = obj;
            save(fullfile(obj.Path2Fin, 'myFinnee.mat'), 'myFinnee')
        end
         
    end
    
end
