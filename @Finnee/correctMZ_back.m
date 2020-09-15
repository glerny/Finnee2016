%% DESCRIPTION

%% EXAMPLES
% myFinnee = myFinnee(
%% COPYRIGHT
% Copyright BSD 3-Clause License Copyright 2016-2017 G. Erny
% (guillaume@fe.up.pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function obj = correctMZ(obj, dts, backgIons, varargin)

%% CORE OF THE FUNCTION
% 1- Initialisation and options

dtsIn            = obj.Datasets{dts};
infoDts          = dtsIn.InfoDts;
infoDts.Path2Dat = {};
options          = checkVarargin(infoDts, varargin{:});

% Create new dat file
[~, rndStr]         = fileparts(tempname);
allProfiles(:,1)    = dtsIn.AxisX.Data;
allProfiles(:,3)    = 0;

% Initiation of profiles
infoDts.ListOfScans = {};
AxisMZ(:,1)         = dtsIn.AxisY.Data;
AxisMZ(:,4)         = 0;
fln                 = 1;
m                   = length(obj.Datasets)+1;
infoDts.Format      = 'profile';

% 2- Load each scan and run the centroid algorithm

infoDts.Title = 'profile dataset';
log2add   = ['PRF=', num2str(m), ' MZcor=true'];

if options.RemSpks
    log2add = [log2add, ' SPR=', num2str(options.SpkSz)];
end

infoDts.Log = [log2add, '|', infoDts.Log];

h = waitbar(0,'processing scans');
for ii = 1:length(dtsIn.AxisX.Data)
    waitbar(ii/size(allProfiles, 1))
    
    infoDts.Path2Dat{fln} = fullfile(obj.Path2Fin, rndStr);
    Scanii  = dtsIn.ListOfScans{ii};
    
    if ~isempty(Scanii.Data)
        
        XMS = dtsIn.AxisY.Data;
        XY  = Scanii.Data;
        XY(:,3) = polyval(backgIons.P(ii,:), XY(:,1));
        XY(:,4) = XY(:,1) - XY(:,3).*XY(:,1);
        XY(:,1) = XY(:,4);
        
        vq =  interp1(XY(:,1), XY(:,2), XMS(:,1),...
            options.mth4interp1);
        XMS(:,2) = vq;
        XMS(isnan(XMS(:,2)), 2) = 0;
        XMS(:,2) = round(XMS(:,2));
        XMS(XMS(:,2) <0, 2) = 0;
        
        % Filter spikes if needed
        if obj.Options.RemSpks
            spkSz = obj.Options.SpkSz;
            XMS   = spikesRemoval(XMS, spkSz );
        end
        
    else
        XMS = dtsIn.AxisY.Data;;
        XMS(:,2) = 0;
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
    infoScn.Precision = 'double';
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
        if s.bytes > obj.Options.MaxFileSize
            [~, rndStr]           = fileparts(tempname);
            fln                   = fln + 1;
            infoDts.Path2Dat{fln} = fullfile(obj.Path2Fin, rndStr);
            infoScn.Path2Dat      = infoDts.Path2Dat{fln};
        end
    end
end
try close(h), catch, end


infoDts.Path2Dat{fln} = fullfile(obj.Path2Fin, rndStr);

% 3- Create the new dataset and save
infoAxis           = dtsIn.AxisX.InfoAxis;
infoAxis.Loc       = 'inFile';
infoAxis.Precision = 'double';
infoAxis.Path2Dat  = infoDts.Path2Dat{fln};
if isempty(allProfiles)
    infoDts.AxisX  = Axis(infoAxis);
else
    infoDts.AxisX  = Axis(infoAxis, allProfiles(:,1));
end

infoAxis           = dtsIn.AxisY.InfoAxis;
infoAxis.Loc       = 'inFile';
infoAxis.Precision = 'double';
infoAxis.Path2Dat  = infoDts.Path2Dat{fln};
if isempty(AxisMZ)
    infoDts.AxisY  = Axis(infoAxis);
else
    infoDts.AxisY  = Axis(infoAxis, AxisMZ(:,1));
end

infoAxis           = dtsIn.AxisY.InfoAxis;
infoAxis.Loc       = 'inFile';
infoAxis.Precision = 'double';
infoAxis.Path2Dat  = infoDts.Path2Dat{fln};
masterMZAxis       = Axis(infoAxis, dtsIn.AxisY.Data);

infoPrf.AxisX      = Axis(dtsIn.AxisX.InfoAxis);
infoPrf.AxisY      = Axis(dtsIn.AxisZ.InfoAxis);
infoPrf.Loc       = 'inFile';
infoPrf.Precision = 'double';
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
    function options = checkVarargin(infoDts, varargin)
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
        
    end
end