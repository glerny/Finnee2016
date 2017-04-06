%% DESCRIPTION 
% The DOCENTROID method allows converting profile scan to centroid scan.
% DOCENTROID can only be used with profile dataset and will create a new 
% dataset. The centroid algorithms and its associated parameters should be 
% defined in the method input parameter. For more information see 
% <https://github.com/glerny/Finnee2016/wiki/Basic-operation-with-Dataset>
%  
%% INPUT PARAMETERS
% *Compulsory 
%   _obj_   : The Finnee object 
%   _dts_   : The indices to the dataset to correct 
%       (i.e. myFinnee.Datasets{dts}
%   _method_: The centroid algorithm to be used. Method should contain the 
%       method's name followed by all needed parameter using the following
%       format: 'methodName:param1:param2:...:paramn'. At this date
%       (20/03/2017) the following methods have been implemented.
%       + 'LocalMax:CN:It'
%           LocalMax detects local maxima within a single profile. The
%           local maxima are recognised as any points whose intensities are
%           higher than their 2xCN closest neighbours (left and right) 
%           (CN is an integer between 1 and 10). If ‘It’ is higher than zero 
%           than only data points higher than It will be considered. The 
%           accurate masses and intensities are calculated at the maximum 
%           value of a polynomial of degree 2 that passes through the 
%           maximum values and its closest neighbours.
% * Optional 
%   None 
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

function obj = doCentroid(obj, dts, method, varargin)

%% CORE OF THE FUNCTION
% 1- Initialisation and options

dtsIn            = obj.Datasets{dts};
infoDts          = dtsIn.InfoDts;
infoDts.Path2Dat = {};
narginchk(2, inf)
if nargin == 2,
    method = 'LocalMax:2:0';
end
[options, MtU]      = checkVarargin(infoDts, method, varargin{:});

% Create new dat file
[~, rndStr]         = fileparts(tempname);
allProfiles(:,1)    = dtsIn.AxisX.Data;
allProfiles(:,3)    = 0;
infoDts.ListOfScans = {};
AxisMZ               = [];
fln                 = 1;
m                   = length(obj.Datasets)+1;
infoDts.Format      = 'centroid';

% 2- Load each scan and run the centroid algorithm
switch MtU{1}
    case 'LocalMax'
        
        infoDts.Title = 'Centroid dataset';
        infoDts.Log   = ['CTR=', num2str(m), ' CTA=', options.method,...
            '|', infoDts.Log];
        h = waitbar(0,'processing scans');
        for ii = 1:size(allProfiles, 1)
            waitbar(ii/size(allProfiles, 1))
            
            infoDts.Path2Dat{fln} = fullfile(obj.Path2Fin, rndStr);
            Scanii          = dtsIn.ListOfScans{ii};
            if ~isempty(Scanii.Data)
                LM          = LocalMaxima(Scanii.Data, MtU{2}, MtU{3});
                if ~isempty(LM)
                MSaccM      = LM(:,1);
                MSaccM(:,2) = round(LM(:,2));
                else
                MSaccM = [];
                end 
            else
                MSaccM = [];
            end
            infoScn           = Scanii.InfoTrc;
            Log               = decipherLog(infoDts.Log, 1);
            infoScn.FT        = Log{1};
            infoScn.Path2Dat  = infoDts.Path2Dat{fln};
            infoScn.Loc       = 'inFile';
            infoScn.Precision = 'single';
            infoScn.TT        = 'CTR';
            infoScn.Title     = ['Centroid scan #', num2str(ii)]; 
            if isempty(MSaccM)
                allProfiles(ii,2) = 0;
                allProfiles(ii,3) = 0;
            else
                allProfiles(ii,2) = sum(MSaccM(:,2));
                allProfiles(ii,3) = max(MSaccM(:,2));
            end
            
            % recorded each scans
            if isempty(MSaccM)
                infoDts.ListOfScans{ii} = Trace(infoScn);
            else
                infoDts.ListOfScans{ii} = Trace(infoScn, MSaccM);
            end
            
            % check the size of the dat file
            s = dir(infoDts.Path2Dat{fln});
            if ~isempty(s)
            if s.bytes > obj.Options.MaxFileSize;
                [~, rndStr]           = fileparts(tempname);
                fln                   = fln + 1;
                infoDts.Path2Dat{fln} = fullfile(obj.Path2Fin, rndStr);
                infoScn.Path2Dat      = infoDts.Path2Dat{fln};
            end
            end
        end
end
try close(h), catch, end

% 3- Create the new dataset and save
infoAxis           = dtsIn.AxisX.InfoAxis;
infoAxis.Loc       = 'inFile';
infoAxis.Precision = 'single';
infoAxis.Path2Dat  = infoDts.Path2Dat{fln};
if isempty(allProfiles)
    infoDts.AxisX  = Axis(infoAxis);
else
    infoDts.AxisX  = Axis(infoAxis, allProfiles(:,1));
end

infoAxis           = dtsIn.AxisY.InfoAxis;
infoAxis.Loc       = 'inFile';
infoAxis.Precision = 'single';
infoAxis.Path2Dat  = infoDts.Path2Dat{fln};
if isempty(AxisMZ)
    infoDts.AxisY  = Axis(infoAxis);
else
    infoDts.AxisY  = Axis(infoAxis, AxisMZ(:,1));
end

infoPrf.AxisX      = Axis(dtsIn.AxisX.InfoAxis);
infoPrf.AxisY      = Axis(dtsIn.AxisZ.InfoAxis);
infoPrf.Loc       = 'inFile';
infoPrf.Precision = 'single';
infoPrf.Path2Dat  = infoDts.Path2Dat{fln};
infoPrf.FT        = infoDts.Log;
infoPrf.TT        = 'SEP';
infoPrf.AdiPrm    = {};

infoPrf.Title     = 'Base Peak Profiles';
if isempty(allProfiles)
    infoDts.BPP   = Trace(infoPrf);
else
    infoDts.BPP   = Trace(infoPrf, [allProfiles(:,1), allProfiles(:,3)]);
end

infoPrf.Title     = 'Total Ion profiles';
if isempty(allProfiles)
    infoDts.TIP   = Trace(infoPrf);
else
    infoDts.TIP   = Trace(infoPrf, [allProfiles(:,1), allProfiles(:,2)]);
end

infoPrf.AxisX      = Axis(dtsIn.AxisY.InfoAxis);
infoPrf.Title     = 'Total Ion Spectrum';
if isempty(AxisMZ)
    infoDts.TIS   = Trace(infoPrf);
else
    infoDts.TIS   = Trace(infoPrf, [AxisMZ(:,1), AxisMZ(:,2)]);
end

infoPrf.Title     = 'Frequency Ion Spectrum';
if isempty(AxisMZ)
    infoDts.FIS   = Trace(infoPrf);
else
    infoDts.FIS   = Trace(infoPrf, [AxisMZ(:,1), AxisMZ(:,3)]);
end

infoPrf.Title     = 'Base Ion Spectrum';
if isempty(AxisMZ)
    infoDts.BIS   = Trace(infoPrf);
else
    infoDts.BIS   = Trace(infoPrf, [AxisMZ(:,1), AxisMZ(:,4)]);
end

infoDts.LAST      = Trace();
infoDts.Option4crt.function = '@Finnee/doCentroid';
infoDts.Option4crt.Options  = options;
obj.Datasets{end+1}         = Dataset(infoDts);
obj.save;


    %% SUB FUNCTIONS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% CHECKVARARGIN
    function [options, MtU] = checkVarargin(infoDts, method, varargin)
        % CHECKVARARGIN is used to check the input paramters
        % and create the options parameter.
                
        % 1- verify the input dataset
        options = {};
        if ~strcmp(infoDts.Format, 'profile')
            error('The target dataset should be in profile mode')
        end
        
        % 2- Verify the method and its parameters
        MtU =  regexp(method, ':', 'split');
        switch lower(MtU{1})
            case 'localmax'
                MtU{1} = 'LocalMax';
                if length(MtU) < 2, MtU{2} = '2'; end
                if length(MtU) < 3, MtU{3} = '0'; end
                MtU{2} = int16(str2double(MtU{2}));
                MtU{3} = int16(str2double(MtU{3}));
                if MtU{2} < 1
                    warning('The first parameter of LocalMax should be higher than 0. Parameters set to 1');
                    MtU{2} = 1;
                end
        end
        options.method = [MtU{1}, ':', num2str(MtU{2}), ':',  num2str(MtU{3})];
        
        % 3- Decipher varargin
    end
end