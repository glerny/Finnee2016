%% DESCRIPTION 
% The BASELINECORRECTION method is used to correct a full profile dataset 
% from baseline drift. BASELINECORRECTION used a graphic user interface 
% (GUI)to select the m/z values that are most relevant for correcting and 
% to select and optimize the baseline correction methods.
% BASELINECORRECTION will create a new dataset. BASELINECORRECTION can be
% used with the original profile dataset. However better results will be
% obtained if the dataset has already been corrected from background noise.
% For more information see 
% <https://github.com/glerny/Finnee2016/wiki/Basic-operation-with-Dataset>
%
%% INPUT PARAMETERS
% *Compulsory 
%   _obj_, the Finnee object 
%   _dts_, the indice to the dts to correct (i.e. myFinnee.Datasets{dts}
% * Optional 
%   _ spikes_ (see the method @Finnee\FilterDataset) for additional 
%       information. Remove spikes in every MS scan. If used, the options 
%       should be followed by 0,1, 2 ot 3 to define the size of spikes to
%       be considered for removal. The default value is 2; 0 allows to turn
%       off spikes removal.
%
%% EXAMPLES
% * myFinnee = myFinnee.BaselineCorrection(2); Will correct dataset 2 for 
%   baseline drift. A new dataset will be created. 
%
%% COPYRIGHT
% Copyright BSD 3-Clause License Copyright 2016-2017 G. Erny
% (guillaume@fe.up.pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function obj = BaselineCorrection(obj, dts, varargin)
    
%% CORE OF THE FUNCTION
% 1- Initialisation and options
dtsIn            = obj.Datasets{dts};
infoDts          = dtsIn.InfoDts;
infoDts.Path2Dat = {};
options          = checkVarargin(infoDts, varargin{:});

% 2- Selection of m/z values 
FIS = dtsIn.FIS.Data(:,2) / size(dtsIn.AxisX.Data, 1)*100;
bool = false;

if isempty(options.FqThres)
    bool = true;
    f1 = figure;
    [N,edges] = histcounts(FIS ,1:1:100);
    bar(edges, [0, N]);
    axis([0 100 0 inf])
    hPlot = gcf;
    title('Frequency profile');
    h = msgbox('Select the frequency threshold for profile to be corrected', 'Correct','modal');
    uiwait(h)
    figure(hPlot);
    [fmax, ~] = ginput(1);
    if isempty(fmax)
        return
    end
    fmax = inputdlg('Confirm the value', 'Correct', 1, {num2str(fmax)});
    options.FqThres = str2double(fmax{1});
    try close(f1), catch, end
end

IndMax = FIS >= options.FqThres;
h = waitbar(0, 'Generating matrix of profiles');
profile(:,1) = dtsIn.AxisX.Data;

% 3- Creating the matrix of profiles to be corrected
for ii = 1:length(profile(:,1))
    waitbar(ii/length(profile(:,1)))
    MS  = dtsIn.xpend(dtsIn.ListOfScans{ii}).Data;
    matRes(:, ii) = MS(IndMax, 2); %#ok<*AGROW>
    profile(ii,2) = sum( MS(IndMax, 2));
    profile(ii,3) = max( MS(IndMax, 2));
    profile(ii,4) = sum( MS(~IndMax, 2));
    profile(ii,5) = max( MS(~IndMax, 2));
end
try close(h), catch, end

if bool
    f2 = figure('Name', 'Selected Profiles for Baseline Corrections');
    subplot(2, 1, 1)
    plot(profile(:,1), profile(:,3));
    title(['Selected profiles for correction (max plot): ', num2str(sum(IndMax))])
    ylabel([dtsIn.AxisZ.Label, ' / ',  dtsIn.AxisZ.Unit]);
    subplot(2, 1, 2)
    plot(profile(:,1), profile(:,5));
    title(['Non-selected profiles (max plot): ', num2str(sum(~IndMax))])
    xlabel([dtsIn.AxisX.Label, ' / ', dtsIn.AxisX.Unit]);
    ylabel([dtsIn.AxisZ.Label, ' / ', dtsIn.AxisZ.Unit]);
    uicontrol('Position',[20 20 200 40],'String','Continue',...
        'Callback','uiresume(gcbf)');
    uiwait(gcf);
    button = questdlg('Are you happy with those conditions?');
    if ~strcmp(button, 'Yes')
        return
    end

    try close(f2), catch, end
end

% 4- Trim the time axis for a more accurate baseline correction
if isempty(options.XLim)
    options.XLim(1) = 1;
    options.XLim(2) = length(dtsIn.AxisX.Data);
end

% 5- Select the baseline method
if isempty(options.baseline)
     warning('off') %#ok<*WNOFF>
     [options.baseline, ExitFlag] = gui4basCor( dtsIn.AxisX, dtsIn.AxisZ, matRes, options.XLim);
     if ExitFlag == -1
         return
     end
     warning('on') %#ok<*WNON>
end

% 6- corecting the profiles
corProf = zeros(size(matRes));
h = waitbar(0,'Correcting profiles');
ind2keep = dtsIn.AxisX.Data >= options.XLim(1) & dtsIn.AxisX.Data <= options.XLim(2);
infoTrc.Title  = '';
infoTrc.FT     = '';
infoTrc.TT     = 'PRF';
infoTrc.AxisX   = Axis(dtsIn.AxisX.InfoAxis);
infoTrc.AxisY   = Axis(dtsIn.AxisX.InfoAxis);
infoTrc.Loc    = 'inTrace';
infoTrc.AdiPrm = {};
infoTrc.P2Fin = '';

for ii = 1:size(corProf, 1)
    waitbar(ii/size(corProf, 1))
    
    prf      = dtsIn.AxisX.Data(ind2keep);
    prf(:,2) = matRes(ii, ind2keep);
    MtU = strsplit(options.baseline, ':');
    IdnZeros = prf(:,2) > 0;
    switch MtU{1}
        case 'PF'
            n = str2double(MtU{2});
            [z, bslPts] = doPF(prf(IdnZeros, :), n);
            
        case 'ArPLS'
            lambda = str2double(MtU{2});
            ratio = str2double(MtU{2}); 
            [z, bslPts] = doArPLS(prf(IdnZeros, 2), lambda, ratio);
            
        case 'ArPLS2'
            lambda = str2double(MtU{2});
            [z, bslPts] = doArPLS2(prf(IdnZeros, 2), lambda);
    end
    
    provPrf      = prf(IdnZeros, 2);
    noise(ii)    = 2*1.96*std(provPrf(bslPts) - z(bslPts));
    Yc           = prf(:,2);
    Yc(IdnZeros) = Yc(IdnZeros) - z ;
    
    Yc                    = round(Yc);
    Yc(Yc < 0)            = 0;
    corProf(ii, ind2keep) = Yc;
end
try close(h), catch, end


%% 7. Recording new dataset and saving
infoDts          = dtsIn.InfoDts;
infoDts.Path2Dat = {};

% Create new dat dile
[~, rndStr]         = fileparts(tempname);
allProfiles(:,1)    = dtsIn.AxisX.Data;
allProfiles(:,3)    = 0;
nnAllProfiles(:,1)  = dtsIn.AxisX.Data;
nnAllProfiles(:,3)  = 0;
infoDts.ListOfScans = {};
AxisMZ              = dtsIn.AxisY.Data;
AxisMZ(:,4)         = 0;
fln                 = 1;
m                   = length(obj.Datasets)+1;
P2PNoise            = dtsIn.AxisY.Data;
P2PNoise(IndMax,2)  = noise;
noise               = sort(noise);
minNoise            = mean(noise(1:round(0.1*length(noise))));
P2PNoise(P2PNoise(:,2) < minNoise ,2)  = minNoise;


infoDts.Title = 'Baseline corrected dataset';
log2add       = ['PRF=', num2str(m), ' BCR=', options.baseline];

if options.RemSpks
    log2add = [log2add, ' SPR=', num2str(options.SpkSz)];
end

h = waitbar(0,'processing scans');
for ii = 1:size(allProfiles, 1)
    waitbar(ii/size(allProfiles, 1))
    
    infoDts.Path2Dat{fln} = fullfile(obj.Path2Fin, rndStr);
    Scanii            = dtsIn.xpend(dtsIn.ListOfScans{ii});
    XMS               = Scanii.Data;
    XMS(IndMax,2)     = corProf(:, ii);
    
    % find and remove spikes
    if options.RemSpks
        spkSz  = options.SpkSz;
        XMS = spikesRemoval(XMS, spkSz );
    end
    
    infoScn           = Scanii.InfoTrc;
    infoScn.FT        = log2add;
    infoScn.Loc       = 'inFile';
    infoScn.Precision = 'single';
    infoScn.Title     = ['Profile scan #', num2str(ii)];
    infoScn.Path2Dat  = infoDts.Path2Dat{fln};
    AxisMZ(:,2)       = AxisMZ(:,2) + XMS(:,2);
    iNZ               = XMS(:,2) > 0;
    AxisMZ(iNZ, 3)    = AxisMZ(iNZ, 3) + 1;
    AxisMZ(:,4)       = max([AxisMZ(:,4), XMS(:,2)], [], 2);
    allProfiles(ii,2) = sum(XMS(:,2));
    allProfiles(ii,3) = max(XMS(:,2));
    
    nnXMS             = XMS(:,2)./P2PNoise(: ,2);
    nnAllProfiles(ii,2) = sum(nnXMS);
    nnAllProfiles(ii,3) = max(nnXMS);
    
    % reduced trailing zero in excess
    XMS = trailRem(XMS, 2);
    
    % recorded each scans
    if isempty(XMS)
        infoDts.ListOfScans{ii} = Trace(infoScn);
    else
        infoDts.ListOfScans{ii} = Trace(infoScn, XMS);
    end
    
    s = dir(infoDts.Path2Dat{fln});
    if isempty(s), continue, end
    if s.bytes > obj.Options.MaxFileSize;
        [~, rndStr]           = fileparts(tempname);
        fln                   = fln + 1;
    end
end

try close(h), catch, end

% creating Dataset
infoAxis           = dtsIn.AxisX.InfoAxis;
infoAxis.Loc       = 'inFile';
infoAxis.Precision = 'single';
infoAxis.Path2Dat  = infoDts.Path2Dat{fln};
AxisMZ             =  trailRem(AxisMZ, 2);

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

infoDts.Log = [log2add, '|', infoDts.Log];
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

infoPrf.Title     = 'Noise normalized base Peak Profiles';
if isempty(nnAllProfiles)
    nnBPP   = Trace(infoPrf);
else
    nnBPP   = Trace(infoPrf, [nnAllProfiles(:,1), nnAllProfiles(:,3)]);
end

infoPrf.Title     = 'Total Ion profiles';
if isempty(allProfiles)
    infoDts.TIP   = Trace(infoPrf);
else
    infoDts.TIP   = Trace(infoPrf, [allProfiles(:,1), allProfiles(:,2)]);
end

infoPrf.Title     = 'Noise normalized total Ion profiles';
if isempty(allProfiles)
    nnTIP   = Trace(infoPrf);
else
    nnTIP   = Trace(infoPrf, [nnAllProfiles(:,1), nnAllProfiles(:,2)]);
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

% record the noise profile
infoPrf.Title     = 'Peak to peak noise';
Peak2PeakNoise    = Trace(infoPrf, P2PNoise);

infoDts.LAST      = Trace();

infoDts.Option4crt.function = 'filterDataset';
infoDts.Option4crt.Options  = options;
obj.Datasets{end+1}         = Dataset(infoDts);
obj.Datasets{end}.AddInfo.Peak2PeakNoise = Peak2PeakNoise;
obj.Datasets{end}.AddInfo.nnBPP          = nnBPP;
obj.Datasets{end}.AddInfo.nnTIP          = nnTIP;
obj.save;

    %% SUB FUNCTIONS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% CHECKVARARGIN
    function options = checkVarargin(infoDts, varargin)
        % CHECKVARARGIN is used to check the input paramters
        % and create the options parameter.
                
        
        % 1- verify the input dataset
        if ~strcmp(infoDts.Format, 'profile')
            error('The target dataset should be in profile mode')
        end
        
        % 2- Default parameters
        options.RemSpks  = true;
        options.SpkSz    = 2;
        options.XLim     = [];
        options.FqThres  = []; 
        options.baseline = {};
        options.RemNoise = false;

        % 3- Decipher varargin
        input = @(x) find(strcmpi(varargin,x),1);
        tgtIx = input('spikes');
        if ~isempty(tgtIx)
            spks = varargin{tgtIx +1};
            if spks == 0
                options.RemSpks = false;
            else
                options.RemSpks = true;
                options.SpkSz   =  spks;
            end
        end
        
        tgtIx = input('XLim');
        if ~isempty(tgtIx)
            XmM = varargin{tgtIx +1};
            options.XLim(1) = min(XmM);
            options.XLim(2) = max(XmM);
            
        end
        
        tgtIx = input('Frequency');
        if ~isempty(tgtIx)
            options.FqThres = varargin{tgtIx +1};
        end
        
        tgtIx = input('Baseline');
        if ~isempty(tgtIx)
            options.baseline = varargin{tgtIx +1};
        end
    end
end
