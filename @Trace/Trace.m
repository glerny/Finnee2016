%% DESCRIPTION
% TRACE is the class that deals with two-dimensional representations. 
% Those can either be electropherograms, chromatograms, MS spectra or 
% others. 
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


classdef Trace
    
    properties
        Title       % Title of the Trace
        FigureTitle % Additional info (normally where the data come from'. 
        % Will be displayed in the figure Title 
        TraceType   % information about the typr of trace: 
        % SEP - Chromatogram, electropherogram or any separation profile 
        % CTR - MS centroided scan; 
        % PRF - MS profile scans
        % OTR - Other
        AxisX        % Information about AxisX (label, unit precision)
        AxisY        % Information about AxisY (label, unit precision)
        Path2Fin     % Link to the Finnee folder
        Precision    % 'inFile' only; Precision in the binary files
        Path2Dat     % 'InFile' only; Address to the file where
        % the data are stored
        Index        % 'inFile' only; Index to the part of the file of interest
        StoredData   % 'inTrace' only; The data
        DataStorage % Define how and where the data are stored
        %   'None'    : No Data
        %   'inTrace' : data are within the class
        %   'inFile'  : data are stored in an associated binary file
        AdiParam    % Place to put parameter for any normalisation/correction
        % should be a structure with two filed, Uses and Values
    end
    
    properties (Dependent)
        Data         % The data
        InfoTrc      % Get back the data of the Axis
        bz           % Bit sizes, depend on the precision of the binary file
    end
    
    methods
        
        function obj = Trace(infoTrc, data2write)
            % Creator method.
            
            if nargin == 0
                obj.Title       = '';
                obj.FigureTitle = '';
                obj.TraceType   = '';
                obj.AxisX       = Axis;
                obj.AxisY       = Axis;
                obj.Path2Fin    = '';
                obj.DataStorage = 'None';
                obj.AdiParam    = {};
                obj.Precision   = '';
                obj.Path2Dat    = '';
                obj.Index       = [0 0];
                obj.StoredData  = [];
                
            else
                obj.Title      	= infoTrc.Title;
                obj.FigureTitle = infoTrc.FT;
                obj.TraceType   = infoTrc.TT;
                obj.AxisX       = infoTrc.AxisX;
                obj.AxisY       = infoTrc.AxisY;
                obj.Path2Fin    = infoTrc.P2Fin;
                obj.DataStorage = infoTrc.Loc;
                obj.AdiParam    = infoTrc.AdiPrm;
                
                % control of data2write
                if nargin == 1, 
                    data2write = [];
               
                end
                if isempty(data2write), obj.DataStorage = 'None'; end
                
                switch lower(obj.DataStorage)
                    case 'none'
                        obj.DataStorage = 'None';
                        obj.Precision   = '';
                        obj.Path2Dat    = '';
                        obj.Index       = [0 0];
                        obj.StoredData  = [];
                        
                    case 'infile'
                        obj.DataStorage = 'inFile';
                        obj.Precision   = infoTrc.Precision;
                        obj.Path2Dat    = infoTrc.Path2Dat;
                        obj.StoredData  = [];
                        
                        % Writting the data
                        fidWriteDat  = fopen(obj.Path2Dat, 'ab');
                        obj.Index(1)= ftell(fidWriteDat);
                        fwrite(fidWriteDat, data2write, obj.Precision);
                        obj.Index(2)= ftell(fidWriteDat);
                        fclose(fidWriteDat);
                        
                    case 'intrace'
                        obj.DataStorage  = 'inTrace';
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
                            [( obj.Index(2)- obj.Index(1))/(n* 2) 2], prcs);
                        fclose(fidReadDat);
                    end
                    
                case 'inTrace'
                    data = obj.StoredData;
            end
        end
        
        % extra method in @trace folder
        % plot(obj)
        
        function InfoTrc = get.InfoTrc(obj)
            % Pass backward InfoAxis
            
            InfoTrc.Title     = obj.Title;
            InfoTrc.FT        = obj.FigureTitle;
            InfoTrc.TT        = obj.TraceType;
            InfoTrc.AxisX     = obj.AxisX;
            InfoTrc.AxisY     = obj.AxisY;
            InfoTrc.P2Fin     = obj.Path2Fin;
            InfoTrc.Loc       = obj.DataStorage;
            InfoTrc.AdiPrm    = obj.AdiParam;
            InfoTrc.Precision = obj.Precision;
        end
    end
    
end

