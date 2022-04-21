function [FOM_MS, FOM_Prof]  = Deconv1HyperROI(HyperROIs, ID)

StoN = 50;
method = 'unconstrained'; %'constrained' or 'unconstrained'

% Load Data

[fidReadDat, errmsg]  = fopen(HyperROIs.Data{ID}.file, 'rb');
cData = reshape(fread(fidReadDat, 'uint64'), HyperROIs.Data{ID}.size);
fclose(fidReadDat);

%%% Align Data with reference profile (time direction)
[aData, Shifts] = doAlignment3D_MinPearson_WithRef(cData, HyperROIs.FOM.Prof{ID}, 10);
HyperROIs.FOM.TimeShift{ID} = Shifts;

% Step1: Deconv MS signals
RefMSSpectra = [HyperROIs.Deconvolution.RefMSSpectra(:, 1)  ...
    mean(HyperROIs.Deconvolution.RefMSSpectra(:, 2:end), 2, 'omitnan')];
[Model, f] = getMSModels(RefMSSpectra, HyperROIs.FOM.MS{ID}, ...
    aData, HyperROIs.Deconvolution.p);
iD2Rem = sum(Model, 1) == 0;
Model(:, iD2Rem) = [];
Model(isnan(Model)) = 0;

Id = sum(Model ~= 0) <= 5;
Model(:, Id) = [];

% Normalise Model
for ii = 1:size(aData, 3)
    switch method
        case 'unconstrained'
            DeconvMS(:, :, ii) = squeeze(aData(:, :, ii))'/Model';
        case 'constrained'
            DeconvMS(:, :, ii) = constrainedDeconv(squeeze(aData(:, :, ii))', Model');
    end
end

FOM_MS.Id_HyperROIs = [];
FOM_MS.MS_M1 = [];
FOM_MS.MS_M2 = [];
FOM_MS.MS_Res_Closest = [];
FOM_MS.MS_max_PearsonCorrelation = [];
FOM_MS.MS_Model_lm = [];
FOM_MS.MS_Model = {};
FOM_MS.ProfNoise = [];
FOM_MS.ProfSignal = [];
FOM_MS.ProfisBckg = [];
for ii = 1:size(DeconvMS, 2)
    XY = [ HyperROIs.FOM.MS{ID}(:, 1), Model(:, ii)];
    M = ChrMoment(XY);
    lm = LocalMaxima(XY, 5, max(XY(:,2))*0.05);
    FOM_MS.Id_HyperROIs(end+1, 1) = ID;
    FOM_MS.MS_M1(end+1, 1) = M(2);
    FOM_MS.MS_M2(end+1, 1) = M(3);
    FOM_MS.MS_Res_Closest(end+1, 1) = NaN;
    FOM_MS.MS_max_PearsonCorrelation(end+1, 1) = NaN;
    FOM_MS.MS_Model_lm(end+1, 1) = size(lm, 1);
    FOM_MS.MS_Model{end+1, 1} = XY;
    FOM_MS.ProfisBckg(end+1, 1) = false;
    Id = findCloser(M(2), XY(:,1));
    if min(Id, length(XY(:,1))-Id) <= 2
        FOM_MS.ProfisBckg(end, 1) = true;
    end
    
    XY = [HyperROIs.Data{ID}.axisTM sum(squeeze(DeconvMS(:, ii, :)), 2)];
    BCKG = diff(XY(:, 2));
    BCKG = MADPeaks([(1:size(BCKG, 1))' BCKG]);
    BCKG = MADPeaks([XY(BCKG(:,1), 1) XY(BCKG(:,1), 2)]);
    FOM_MS.ProfNoise(end+1, 1) = 3*std(BCKG(:,2));
    FOM_MS.ProfSignal(end+1, 1) = max(XY(:,2));
    
    if FOM_MS.ProfSignal(end, 1)/FOM_MS.ProfNoise(end, 1) < StoN
        FOM_MS.ProfisBckg(end, 1) = true;
    end
end

% Calculate missing FOM
FOM_MS = struct2table(FOM_MS);
[FOM_MS, Id4Sort] = sortrows(FOM_MS, 'MS_M1');
DeconvMS = DeconvMS(:, Id4Sort, :);
sum_DeconvMS = squeeze(sum(DeconvMS, 3));

if size(FOM_MS, 1) > 1
    for ii = 1:size(FOM_MS, 1)
        CC = corrcoef(sum_DeconvMS);
        FOM_MS.MS_max_PearsonCorrelation(ii) = max(CC(1:end ~= ii, ii));
        
        if ii > 1
            FOM_MS.MS_Res_Closest(ii) = (FOM_MS.MS_M1(ii) - FOM_MS.MS_M1(ii-1))/(4*sqrt(FOM_MS.MS_M2(ii)));
        end
    end
end

id2deconv = find(~FOM_MS.ProfisBckg);
FOM_Prof.Id = [];
FOM_Prof.Id_HyperROIs = [];
FOM_Prof.Id_MS = [];
FOM_Prof.MS_M1 = [];
FOM_Prof.MS_M2 = [];
FOM_Prof.MS_Model = {};
FOM_Prof.Prof_M1 = [];
FOM_Prof.Prof_M2 = [];
FOM_Prof.Prof_Res_Closest = {};
FOM_Prof.Prof_max_PearsonCorrelation = [];
FOM_Prof.Area_max_PearsonCorrelation = [];
FOM_Prof.DeconvArea = {};
FOM_Prof.DeconvRSSR = {};
FOM_Prof.Prof_Model = {};
FOM_Prof.Prof_Fitting = {};

if ~isempty(id2deconv)
    for jj = 1:length(id2deconv)
        cData = [HyperROIs.Data{ID}.axisTM squeeze(DeconvMS(:, id2deconv(jj), :))];
        % alData = doAlignment_MinPearson(cData, 20);
        
        % Model = fittPMG2(alData);
        % Model = fittPLMG(alData);
        [Model, PeakModels] =  fittPMG(cData, true);
        Areas = (cData(:, 2:end))'/Model';
        Predict = Areas*Model';
        SSR = sum((cData(:,2:end) - Predict').^2);
        
        for kk = 1:size(Model, 2)
            FOM_Prof.Id(end+1, 1) = length(FOM_Prof.Id)+1;
            FOM_Prof.Id_HyperROIs(end+1, 1) = ID;
            FOM_Prof.Id_MS(end+1, 1) = id2deconv(jj);
            FOM_Prof.MS_M1(end+1, 1) = FOM_MS.MS_M1(id2deconv(jj));
            FOM_Prof.MS_M2(end+1, 1) = FOM_MS.MS_M2(id2deconv(jj));
            FOM_Prof.MS_Model{end+1, 1} = FOM_MS.MS_Model{id2deconv(jj)};
            
            XY = [HyperROIs.Data{ID}.axisTM Model(:, kk)];
            M = ChrMoment(XY);
            FOM_Prof.Prof_M1(end+1, 1) = M(2);
            FOM_Prof.Prof_M2(end+1, 1) = M(3);
            FOM_Prof.Prof_Res_Closest{end+1, 1} = [];
            FOM_Prof.Prof_max_PearsonCorrelation(end+1, 1) = NaN;
            FOM_Prof.Area_max_PearsonCorrelation(end+1, 1) = NaN;
            FOM_Prof.DeconvArea{end+1, 1} = Areas(:, kk);
            FOM_Prof.DeconvRSSR{end+1, 1} = SSR;
            FOM_Prof.Prof_Model{end+1, 1} = XY;
            FOM_Prof.Prof_Fitting{end+1, 1} = PeakModels{kk};
        end
    end
end

% Calculate missing FOM
FOM_Prof = struct2table(FOM_Prof);
%%%% MISSING FOM TO BE DONE %%%%

