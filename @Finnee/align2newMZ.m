%% DESCRIPTION
% ALIGN2NEWMZ is used to normalise every MS scans in a datasets to a
% common mz axis. 
% 
%% INPUT PARAMETERS
% *Compulsory*
% *obj      : The Finnee object
% *dts      : The indice to the dataset that contains the original scans
%           (should be profile scans)
% *newAxis  : The mz axis that will be used as reference for every scans.
%           This parameter can be empty, in this case the optional parameter 
%           'masterAxis' should be used to defined the new mz Axis (see
%           *optional')
%
% *Optional*
% *tLim       : Followed by a 2x1 array of numbers (default [0 inf]). Only 
%             records scans between tLim(1) and tLim(2)
% *spikes     : Followed by an integer between 0 and 3 (default 2) (see the 
%             method @Finnee\FilterDataset) for additional information. 
%             Remove spikes in every MS scans If used, where spikes are any 
%             peaks in each MS of length equal or lower that the integer. 
%             'spikes' followed by 0 allows to to turn off spikes removal.
% *masterAxis': Followed by a string in the form 'mzStart:mzEnd:Id'.
%             mzStart and mzEnd will be used to defined the limit of the
%             axis, Id is the indice to the scen that will be used to
%             generate the mzAxis (used max if you want to use the most
%             aboundant scans). For more information see help at
%             @Trace\extrapolMZ
% *meth4int'  : Followed by a string to specify an alternative 
%             interpolation method: 'nearest', 'next', 'previous', 
%             'linear','spline','pchip', or 'cubic'. The default method is 
%             'linear'. help interp1 for more information on the
%             interpolation function. Changing this option is not advised.
%
%% OUTPUT PARAMETERS
% *obj      : The Finnee object. !Important, if not input parameter is used
%           the saved Finnnee object will still be modified.
%
%% EXAMPLES
% myFinnee = myFinnee(
%% COPYRIGHT
% Copyright BSD 3-Clause License Copyright 2016-2017 G. Erny
% (guillaume@fe.up.pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function obj = align2newMZ(obj, dts, newAxis, varargin)

%% CORE OF THE FUNCTION
% 1- Initialisation and options

dtsIn            = obj.Datasets{dts};
infoDts          = dtsIn.InfoDts;
infoDts.Path2Dat = {};
options          = checkVarargin(infoDts, newAxis, varargin{:});

% Create new dat file
[~, rndStr]         = fileparts(tempname);
allProfiles(:,1)    = dtsIn.AxisX.Data;
allProfiles(:,3)    = 0;

% Load or create master mz Axis
if options.MstAxis
    M4MA = strsplit(options.mth4mstAxis, ':');
    
    if strcmpi(M4MA{3}, 'max');
        [~,Id4MA] = max(dtsIn.TIP.Data(:,2));
    else
        Id4MA = str2double(M4MA{3});
    end
    
    MZlim(1) = min( str2double(M4MA{1}), str2double(M4MA{2}));
    MZlim(2) = max( str2double(M4MA{1}), str2double(M4MA{2}));
    
    newAxis = dtsIn.ListOfScans{Id4MA}.extrapolMZ(2, MZlim);
end

% Initiation of profiles
infoDts.ListOfScans = {};
AxisMZ(:,1)         = newAxis;
AxisMZ(:,4)         = 0;
fln                 = 1;
m                   = length(obj.Datasets)+1;
infoDts.Format      = 'profile';

% 2- Load each scan and run the centroid algorithm

infoDts.Title = 'profile dataset';
log2add   = ['PRF=', num2str(m), ' MMZ=true'];

if options.RemSpks
    log2add = [log2add, ' SPR=', num2str(options.SpkSz)];
end

infoDts.Log = [log2add, '|', infoDts.Log];

h = waitbar(0,'processing scans');
for ii = 1:length(dtsIn.AxisX.Data)
    waitbar(ii/size(allProfiles, 1))
    
    infoDts.Path2Dat{fln} = fullfile(obj.Path2Fin, rndStr);
    Scanii  = dtsIn.ListOfScans{ii};
    Id2Keep = Scanii.Data(:,1) >= newAxis(1) & Scanii.Data(:,1) <= newAxis(end);
    
    
    XMS = newAxis;
    XY  = Scanii.Data(Id2Keep, :);
    [~, ia, ~] = unique(XY(:,1));
    XY  = XY(ia, :);
    
    try
    XMS(:,2) =  interp1(XY(:,1), XY(:,2), XMS(:,1),...
        options.mth4interp1);
    catch
        [C, ia, ic] = unique(Scanii.Data(Id2Keep,1));
        disp('wtf')
    end
    XMS(isnan(XMS(:,2)), 2) = 0;
    XMS(:,2) = round(XMS(:,2));
    XMS(XMS(:,2) <0, 2) = 0;
    
    % Filter spikes if needed
    if obj.Options.RemSpks
        spkSz = obj.Options.SpkSz;
        XMS   = spikesRemoval(XMS, spkSz );
    end
    
    AxisMZ(:, 2) = AxisMZ(:,2) + XMS(:,2);
    AxisMZ(XMS(:,2) > 0, 3) = AxisMZ(XMS(:,2) > 0, 3) + 1;
    AxisMZ(:, 4) = max(AxisMZ(:,4), XMS(:,2));
    
    % reduced trailing zero in excess
    XMS = trailRem(XMS, 2);
    
    infoScn           = Scanii.InfoTrc;
    infoScn.FT        = log2add;
    infoScn.Path2Dat  = infoDts.Path2Dat{fln};
    infoScn.Loc       = 'inFile';
    infoScn.Precision = 'single';
    infoScn.TT        = 'PRF';
    infoScn.Title     = ['Profile scan #', num2str(ii)];
    infoScn.AdiPrm    = {};
    
    if isempty(XMS)
        allProfiles(ii,2) = 0;
        allProfiles(ii,3) = 0;
    else
        allProfiles(ii,2) = sum(XMS(:,2));
        allProfiles(ii,3) = max(XMS(:,2));
    end
    
    % recorded each scans
    if isempty(XMS)
        infoDts.ListOfScans{ii} = Trace(infoScn);
    else
        infoDts.ListOfScans{ii} = Trace(infoScn, XMS);
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
try close(h), catch, end


infoDts.Path2Dat{fln} = fullfile(obj.Path2Fin, rndStr);
% reduced trailing zero in AxisMZ to
AxisMZ        =  trailRem(AxisMZ, 2);


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

infoAxis           = dtsIn.AxisY.InfoAxis;
infoAxis.Loc       = 'inFile';
infoAxis.Precision = 'single';
infoAxis.Path2Dat  = infoDts.Path2Dat{fln};
masterMZAxis       = Axis(infoAxis, newAxis);

infoPrf.AxisX      = Axis(dtsIn.AxisX.InfoAxis);
infoPrf.AxisY      = Axis(dtsIn.AxisZ.InfoAxis);
infoPrf.Loc       = 'inFile';
infoPrf.Precision = 'single';
infoPrf.Path2Dat  = infoDts.Path2Dat{fln};
infoPrf.P2Fin     = obj.Path2Fin;
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
infoDts.Option4crt.function = '@Finnee/doAxeMZEqualized';
infoDts.Option4crt.Options  = {};
obj.Datasets{end+1}         = Dataset(infoDts);
obj.Datasets{end}.AddInfo   = {};
obj.Datasets{end}.AddInfo.masterMZAxis = masterMZAxis;
obj.save;


%% SUB FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CHECKVARARGIN
    function options = checkVarargin(infoDts, newAxis, varargin)
        % CHECKVARARGIN is used to check the input paramters
        % and create the options parameter.
        
        % 1- verify the input dataset
        options = {};
        if ~strcmp(infoDts.Format, 'profile')
            error('The target dataset should be in profile mode')
        end
        
        % 2- Defaults optional parameters
        options.RemSpks     = true;
        options.SpkSz       = 2;
        options.mth4interp1 = 'linear';
        options.MstAxis     = false;
        options.mth4mstAxis = '';
        options.XLim        = [0 inf];
        
        % 3- Decipher varargin
        input = @(x) find(strcmpi(varargin,x),1);
        
        tgtIx = input('tLim');
        if ~isempty(tgtIx)
            tLim         = varargin{tgtIx +1};
            options.XLim = [min(tLim) max(tLim)];
        end
        
        tgtIx = input('spikes');
        if ~isempty(tgtIx)
            spks = varargin{tgtIx +1};
            if spks == 0
                options.RemSpks = false;
            else
                options.RemSpks = true;
                options.SpkSz  =  spks;
            end
        end
        
        if isempty(newAxis)
            options.MstAxis      = true;
            options.mth4mstAxis  =  [num2str(infoDts.MZlim(1)), ':', ...
              num2str(infoDts.MZlim(2)), ':max'] ;
        end
        
        tgtIx = input('masterAxis');
        if ~isempty(tgtIx)
            options.MstAxis      = true;
            options.mth4mstAxis  =  varargin{tgtIx +1};
        end
        
    end
end