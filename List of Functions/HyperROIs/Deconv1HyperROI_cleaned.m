function [FOM_HypROIs, FOM_MS, FOM_Prof]  = Deconv1HyperROI_cleaned(HyperROIs, ID, Saved)

%%% HYPER_ROIs Inside paramters (can be changed via options (TO BE DONE)
StoN                   = 5;
MS_Max_Resolution      = 0.5;
MS_Shift4Loop          = 0.95;
MSPtsPerPeaks          = 10;
ProfmaxShift           = 0.5; %min
Profile_Max_Resolution = 1.5;
ProfPtsPerPeaks        = 20;
Prof_Shift4Loop        = 0.95;
prct4noise             = 25;
minPtsPerPeaks         = 10;

% MODEL MS SPECTRA FOR DECONVOLUTION
RefMSSpectra = [HyperROIs.Deconvolution.RefMSSpectra(:, 1)  ...
    mean(HyperROIs.Deconvolution.RefMSSpectra(:, 2:end), 2, 'omitnan')];
% Initialisation

FOM_MS.Id_HyperROIs = [];
FOM_MS.MS_M1 = [];
FOM_MS.MS_M2 = [];
FOM_MS.MS_Crit1 = [];
FOM_MS.MS_Crit2 = {};
FOM_MS.MS_Model = {};
FOM_MS.ProfNoise = [];
FOM_MS.ProfSignal = [];
FOM_MS.ProfisBckg = [];
FOM_MS.ProfileMoments = [];
FOM_MS.shift = {};

FOM_Prof.Id_HyperROIs = [];
FOM_Prof.MS_M1 = [];
FOM_Prof.MS_Crit1 = [];
FOM_Prof.Prof_M1 = [];
FOM_Prof.Prof_M2 = [];
FOM_Prof.Prof_Crit1 = [];
FOM_Prof.Prof_Crit2 = {};
FOM_Prof.DeconvArea = {};
FOM_Prof.maxArea = [];
FOM_Prof.DeconvolutedProf = {};
FOM_Prof.Prof_Fitting = {};

% Load Data
FOM_HypROIs.ID = ID;
try
    [fidReadDat, errmsg]  = fopen(HyperROIs.Data{ID}.file, 'rb');
    aData = reshape(fread(fidReadDat, 'uint64'), HyperROIs.Data{ID}.size);
    fclose(fidReadDat);
    FOM_HypROIs.ErrorDataLoading = errmsg;
    FOM_HypROIs.MSnormalisation  = {};
    AxisTime = HyperROIs.Data{ID}.axisTM;
    AxisMz   = HyperROIs.Data{ID}.axisMZ;
    FOM_HypROIs.Moments = table();
    for ii = 1:size(aData, 3)
        FOM_HypROIs.Moments = [FOM_HypROIs.Moments; ...
            ChrMoment3D(AxisTime, AxisMz, squeeze(aData(:, :, ii)))];
    end
    FOM_HypROIs.Error = {};
    
catch ME
    FOM_HypROIs.ErrorDataLoading = {};
    FOM_HypROIs.MSnormalisation  = {};
    FOM_HypROIs.Moments = table();
    FOM_HypROIs.nbrMSPeaks = NaN;
    FOM_HypROIs.Error = ME;
    FOM_HypROIs = struct2table(FOM_HypROIs, 'AsArray', true);
    FOM_Prof = struct2table(FOM_Prof);
    FOM_MS = struct2table(FOM_MS);
    return
end
FOM_HypROIs.Parameters.StoN = StoN;
FOM_HypROIs.Parameters.MS_Max_Resolution = MS_Max_Resolution;
FOM_HypROIs.Parameters.MS_Shift4Loop = MS_Shift4Loop;
FOM_HypROIs.Parameters.MSPtsPerPeaks = MSPtsPerPeaks;
FOM_HypROIs.Parameters.ProfmaxShift = ProfmaxShift;
FOM_HypROIs.Parameters.Profile_Max_Resolution = Profile_Max_Resolution;
FOM_HypROIs.Parameters.ProfPtsPerPeaks = ProfPtsPerPeaks;
FOM_HypROIs.Parameters.Prof_Shift4Loop = Prof_Shift4Loop;
FOM_HypROIs.Parameters.prct4noise = prct4noise;
FOM_HypROIs.Parameters.minPtsPerPeaks = minPtsPerPeaks;
FOM_HypROIs.Parameters.Saved = Saved;

FOM_HypROIs = struct2table(FOM_HypROIs, 'AsArray', true);

% Normalise Model

% Step1: Deconv MS signals
[Model, FttParameters] = getMSModels_constrained_cleaned(RefMSSpectra,...
    HyperROIs.FOM.MS{ID}, aData, MS_Max_Resolution,...
    MS_Shift4Loop, MSPtsPerPeaks);
FOM_HypROIs.nbrMSPeaks = size(Model, 2);
if isempty(Model)
    FOM_Prof = struct2table(FOM_Prof);
    FOM_MS = struct2table(FOM_MS);
    return; 
end

[Crit1, Crit2] = getPeakPurity(Model, squeeze(sum(aData, 2, 'omitnan')));
for ii = 1:size(aData, 3)
    DeconvMS(:, :, ii) = squeeze(aData(:, :, ii))'/Model';
end

for ii = 1:size(DeconvMS, 2)
    [Aligned_Data, shift] = doAlignment_MinPearson_directed...
        ([HyperROIs.Data{ID}.axisTM squeeze(DeconvMS(:, ii, :))]...
        , ProfmaxShift, strcmp(HyperROIs.Files.Tags, 'QC'));
    DeconvMS(:, ii, :) = Aligned_Data(:, 2:end);
    FOM_MS.Id_HyperROIs(end+1, 1) = ID;
    FOM_MS.MS_M1(end+1, 1) = FttParameters(ii, 1);
    FOM_MS.MS_M2(end+1, 1) = FttParameters(ii, 2);
    FOM_MS.MS_Crit1(end+1, 1) = Crit1(ii);
    FOM_MS.MS_Crit2{end+1, 1} = squeeze(Crit2(:, :, ii));
    flatData = squeeze(DeconvMS(:, ii, :));
    for jj = 1:size(flatData, 2)
        M(jj, :) = ChrMoment([AxisTime, flatData(:, jj)]); 
    end
    FOM_MS.MS_Model{end+1, 1} = [AxisMz, Model(:, ii)];
    FOM_MS.ProfileMoments{end+1, 1} = M;
    FOM_MS.ProfisBckg(end+1, 1) = false;
    
    %%% FIRST CONDITION FOR BACKGROUND CENTER OF THE MAX MS PEAK NEAR EDGES
    Id = findCloser(FttParameters(ii, 1), AxisMz);
    if min(Id, length(AxisMz)-Id) <= 2
        FOM_MS.ProfisBckg(end, 1) = true;
    end
    
    XY = [AxisTime mean(flatData, 2, 'omitnan')];
    id = XY(:,2) ~= 0;
    if sum(id) < minPtsPerPeaks
        FOM_MS.ProfisBckg(end, 1) = true;
        FOM_MS.ProfNoise(end+1, 1) = NaN;
        FOM_MS.ProfSignal(end+1, 1) = NaN;
    else
        ix = XY(:,2) < prctile(XY(:,2), prct4noise);
        FOM_MS.ProfNoise(end+1, 1) = 3*std(XY(ix, 2));
        FOM_MS.ProfSignal(end+1, 1) = max(XY(:,2));
    end
    
    if FOM_MS.ProfSignal(end, 1)/FOM_MS.ProfNoise(end, 1) < StoN
        FOM_MS.ProfisBckg(end, 1) = true;
    end
    FOM_MS.shift{end+1, 1} = shift;
end

% Calculate missing FOM
FOM_MS = struct2table(FOM_MS);

id2deconv = find(~FOM_MS.ProfisBckg);
if ~isempty(id2deconv)
    for jj = 1:length(id2deconv)
        cData = [HyperROIs.Data{ID}.axisTM squeeze(DeconvMS(:, id2deconv(jj), :))];
        [Model, PeakModels] =  fittPMG_cleaned(cData...
            , Profile_Max_Resolution...
            , ProfPtsPerPeaks...
            , Prof_Shift4Loop);
        
        Areas = (cData(:, 2:end))'/Model';
        [Crit1, Crit2, DeconvSignal] = getPeakPurity(Model, cData(:, 2:end));
        
        for kk = 1:size(Model, 2)
            FOM_Prof.Id_HyperROIs(end+1, 1) = ID;
            FOM_Prof.MS_M1(end+1, 1) = FOM_MS.MS_M1(id2deconv(jj));
            FOM_Prof.MS_Crit1(end+1, 1) = FOM_MS.MS_Crit1(id2deconv(jj));
            
            id = Model(:, kk) > 0.0001*max(Model(:, kk));
            DS = squeeze(DeconvSignal(id, kk, :));
            XY = [AxisTime(id) mean(DS, 2, 'omitnan')];
            M = ChrMoment(XY);
            FOM_Prof.Prof_M1(end+1, 1) = M(2);
            FOM_Prof.Prof_M2(end+1, 1) = M(3);
            FOM_Prof.DeconvArea{end+1, 1} = Areas(:, kk);
            FOM_Prof.maxArea(end+1, 1) = max(Areas(:, kk));
            FOM_Prof.Prof_Crit1(end+1, 1) = Crit1(kk);
            FOM_Prof.Prof_Crit2{end+1, 1} = squeeze(Crit2(:, :, kk));
            
            switch Saved
                case 'small'
                    FOM_Prof.DeconvolutedProf{end+1, 1} = [];
                    
                case 'all'
                    
                    FOM_Prof.DeconvolutedProf{end+1, 1} = [AxisTime(id) DS];
            end
            FOM_Prof.Prof_Fitting{end+1, 1} = PeakModels{kk};
        end
    end
end
FOM_Prof = struct2table(FOM_Prof);

