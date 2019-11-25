%% DESCRIPTION
% AXIS is the class of the Finnee2016 toolbox that contain all 
% information associated with an axis.
%
%% LIST OF THE CLASS'S PROPERTIES
% *Label*          : Time, intensity, mass...
% *Unit*           : min, m/z...
% *DecimPlace*     : Number of decimal places to use when displaying data
%   related to this axis
%
%% LIST OF THE CLASS'S METHODS:
%
%% COPYRIGHT
% Copyright BSD 3-Clause License Copyright 2016-2017 G. Erny
% (guillaume@fe.up.pt), FEUP, Porto, Portugal
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef Axis
    
    properties
        Label   	 % Time, Intensity or other
        Unit         % min, s...
        DecimPlaces  % Only for display. Number of desimals
        DataStorage  % Define how and where the data are stored
                     %  'None'  : No Data
                     %  'inAxis' : data are within the class
                     %  'inFile': data are stored in an associated binary file
    end
    
    properties (Hidden = true, SetAccess = immutable)
        Precision    % 'inFile' only; Precision in the binary files
        Path2Dat     % 'InFile' only; Address to the file where
                     % the data are stored
        Index        % 'inFile' only; Index to the part of the file of interest
        StoredData   % 'inAxis' only; The data
        
    end
    
    properties (Dependent)
        Data         % The data  
        InfoAxis      % Get back the data of the axis
    end
    
    properties (Dependent, Hidden)
        bz           % Bit sizes, depend on the precision of the binry file    
        fo           % formatting operator from DecimPlaces
    end
    
    methods
        function obj = Axis(infoAxis, data2write)
            % Creator method.
            
            if nargin == 0
                obj.Label        = '';
                obj.Unit         = '';
                obj.DecimPlaces  = 0 ;
                obj.DataStorage  = 'None';
                obj.Precision    = '';
                obj.Path2Dat     = '';
                obj.Index        = [0 0];
                obj.StoredData   = [];
                
            else
                obj.Label       = infoAxis.Label;
                obj.Unit        = infoAxis.Unit;
                obj.DecimPlaces = infoAxis.dp;
                obj.DataStorage = infoAxis.Loc;
                
                % control of data2write
                if nargin == 1
                    obj.DataStorage = 'None';
                elseif isempty(data2write)
                    obj.DataStorage = 'None';
                end
                
                switch lower(obj.DataStorage)
                    case 'none'
                        obj.DataStorage = 'None';
                        obj.Precision   = '';
                        obj.Path2Dat    = '';
                        obj.Index       = [0 0];
                        obj.StoredData  = [];
                        
                    case 'infile'
                        obj.DataStorage = 'inFile';
                        obj.Precision   = infoAxis.Precision;
                        obj.Path2Dat    = infoAxis.Path2Dat;
                        obj.StoredData  = [];
                        
                        %write the data
                        fidWriteDat = fopen(obj.Path2Dat, 'ab');
                        obj.Index(1)= ftell(fidWriteDat);
                        fwrite(fidWriteDat, data2write, obj.Precision);
                        obj.Index(2)= ftell(fidWriteDat);
                        fclose(fidWriteDat);
                        
                    case 'inaxis'
                        obj.DataStorage  = 'inAxis';
                        obj.Precision    = '';
                        obj.Path2Dat     = '';
                        obj.Index        = [0 0];
                        
                        %write the data
                        obj.StoredData   = data2write;
                        
                    otherwise
                        error('Incorrect DataLocation type')
                end
            end
        end
        
        function bz = get.bz(obj)
            % determine number of bite sizes for fread
            
            switch obj.Precision
                case ''
                    bz = 0;
                case {'double', 'int64', 'uint64'}
                    bz = 8;
                case {'single', 'int32', 'uint32'}
                    bz = 4;
                case {'int16', 'uint16'}
                    bz = 2;
                case {'int8', 'uint8'}
                    bz = 1;
            end     
        end
        
        function fo = get.fo(obj)
            % make the adhoc formating opperator 
            if obj.DecimPlaces == 0
                fo = '%i';
            else
                fo = ['%.', num2str(obj.DecimPlaces), 'f'];
            end
        end
        
        function data = get.Data(obj)
            % recover the data depending on the way they are stored 
            
            switch obj.DataStorage
                case'None'
                    data = [];
                    
                case 'inFile'
                    if obj.Index(1) == obj.Index(2)
                        data = [];
                        
                    else
                        n    = obj.bz;
                        prcs = obj.Precision;
                        fidReadDat = fopen(obj.Path2Dat, 'rb');
                        fseek(fidReadDat,  obj.Index(1), 'bof');
                        data = fread(fidReadDat, ...
                            [( obj.Index(2)- obj.Index(1))/n 1], prcs);
                        fclose(fidReadDat);
                    end
                    
                case 'inAxis'
                    data = obj.StoredData;
            end
        end
        
        function infoAxis = get.InfoAxis(obj)
            % Pass backward InfoAxis
            
            infoAxis.Label     = obj.Label;
            infoAxis.Unit      = obj.Unit;
            infoAxis.dp        = obj.DecimPlaces;
            infoAxis.Loc       = obj.DataStorage;
            infoAxis.Precision = obj.Precision;
            infoAxis.Path2Dat  = obj.Path2Dat;
            
        end
        
        function  obj = set.Label(obj, value)
            % standardised the labels string (need to be improved)
            
            if isempty(value)
                obj.Label = '';
            else
                switch lower(value)
                    case {'time'}
                        obj.Label = 'Time';
                        
                    case {'mass'}
                        obj.Label = 'Mass';
                        
                    case {'intensity'}
                        obj.Label = 'Intensity';
                end
            end
        end
        
        function  obj = set.Unit(obj, value)
            % standardised the units string (need to be improved)

            switch lower(value)
                case {'minutes', 'min', 'minute'}
                    obj.Unit = 'min';
                    
                case {'secondes', 'sec', 's'}
                    obj.Unit = 's';
                    
                case{'number of detector counts', 'counts'}
                    obj.Unit = 'counts';
                    
                case{'m/z'}
                    obj.Unit = 'm/z';
                    
                otherwise
                    obj.Unit = '';
                    
            end
        end
    end
end
    
