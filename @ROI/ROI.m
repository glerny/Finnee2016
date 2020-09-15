%% DESCRIPTION
% ROI is the class that deals with three-dimensional representations.
% Those can a cut from the dataset after mzAxis alignment to the Matser mz Axis
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


classdef ROI
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Title        % Title of the ROI
        TagOfDts     % Additional info (normally where the data come from'.
        TgtMz        % Target mz value
        MzWdw        % mz window in datapts (full interval 1+2*mzWdw
        TgtTm        % Target Time
        TgtVar       % Variance of the target Peak
        tmWdW        % time wdw
        AxisTm       % Time Axis
        AxisMZ       % MZ Axis
        Path2Fin     % Link to the Finnee folder
        StoredData   % The data
        %Filter  = {-3:3,-3:3,2,2};
        Noise = 0;
    end
    
    properties (Dependent)
        InfoROI     % Get back the data of the Axisthe binary file
        %Smoothed
        %Derived
    end
    
    
    methods
        function obj = ROI(InfoROI)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            if nargin == 0
                obj.Title      = '';
                obj.TagOfDts   = '';
                obj.TgtMz      = NaN;
                obj.MzWdw      = NaN;
                obj.TgtTm      = NaN;
                obj.tmWdW      = NaN;
                obj.TgtVar     = NaN;
                obj.AxisTm     = Axis;
                obj.AxisMZ     = Axis;
                obj.Path2Fin   = '';
                obj.StoredData = [];
                
            else
                obj.Title      = InfoROI.Title;
                obj.TagOfDts   = InfoROI.TagOfDts;
                obj.TgtMz      = InfoROI.TgtMz;
                obj.MzWdw      = InfoROI.MzWdw;
                obj.TgtTm      = InfoROI.TgtTm;
                obj.tmWdW      = InfoROI.tmWdW;
                obj.AxisTm     = InfoROI.AxisTm;
                obj.TgtVar     = InfoROI.TgtVarnc;
                obj.AxisMZ     = InfoROI.AxisMZ;
                obj.Path2Fin   = InfoROI.Path2Fin;
                obj.StoredData = InfoROI.StoredData;
                
            end
        end
        
        function InfoROI = get.InfoROI(obj)
            % Pass backward InfoAxis
            
            InfoROI.Title      = obj.Title;
            InfoROI.TagOfDts   = obj.TagOfDts;
            InfoROI.TgtMz      = obj.TgtMz;
            InfoROI.MzWdw      = obj.MzWdw;
            InfoROI.TgtTm      = obj.TgtTm;
            InfoROI.tmWdW      = obj.tmWdW;
            InfoROI.AxisTm     = obj.AxisTm;
            InfoROI.AxisMZ     = obj.AxisMZ;
            InfoROI.Path2Fin   = obj.Path2Fin;
            InfoROI.StoredData = obj.StoredData;
            InfoROI.SmoothWdw  = obj.SmoothWdw;
        end
        
%         function Smoothed = get.Smoothed(obj)
%             
%             h = sgsdf_2d(obj.Filter{1},obj.Filter{2},obj.Filter{3},obj.Filter{4},0,0);
%             Smoothed = imfilter(obj.StoredData,h);
%         end
%         
%         function Derived = get.Derived(obj)
%             
%             h = sgsdf_2d(obj.Filter{1},obj.Filter{2},obj.Filter{3},obj.Filter{4},1,0);
%             Derived = imfilter(obj.StoredData,h);
%         end
        
    end
end

