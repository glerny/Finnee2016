function [FOM_MS, FOM_Prof, HRShift]  = Deconv1HyperROI(HyperROIs, ID)

StoN = 100;
maxShift = 0.5; %min
method = 'constrained'; %'fullconstrained, 'constrained' or 'unconstrained'
% Initialisation

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
FOM_MS.isMerged = [];
FOM_MS.CenterTime = [];
FOM_MS.shift = {};

FOM_Prof.Id = [];
FOM_Prof.Id_HyperROIs = [];
FOM_Prof.Id_MS = [];
FOM_Prof.MS_M1 = [];
FOM_Prof.MS_M2 = [];
FOM_Prof.isMergedMS = [];
FOM_Prof.MS_Model = {};
FOM_Prof.Prof_M1 = [];
FOM_Prof.Prof_M2 = [];
FOM_Prof.isMergedProf = [];
FOM_Prof.DeconvArea = {};
FOM_Prof.maxArea = [];
FOM_Prof.DeconvRSSR = {};
FOM_Prof.Prof_Model = {};
FOM_Prof.Prof_Fitting = {};

% Load Data

[fidReadDat, errmsg]  = fopen(HyperROIs.Data{ID}.file, 'rb');
cData = reshape(fread(fidReadDat, 'uint64'), HyperROIs.Data{ID}.size);
fclose(fidReadDat);

%%% Align Data with reference profile (time direction)
aData = cData;
%[aData, Shifts] = doAlignment3D_MinPearson_WithRef(cData, HyperROIs.FOM.Prof{ID}, 150, 1);
%HyperROIs.FOM.TimeShift{ID} = Shifts;

%%% Align Data with reference profile (mz direction)
%[aData, Shifts] = doAlignment3D_MinPearson_WithRef(aData, HyperROIs.FOM.MS{ID}, 1, 2);
%HyperROIs.FOM.mzShift{ID} = Shifts;

% Step1: Deconv MS signals
RefMSSpectra = [HyperROIs.Deconvolution.RefMSSpectra(:, 1)  ...
    mean(HyperROIs.Deconvolution.RefMSSpectra(:, 2:end), 2, 'omitnan')];

% Step2: Aligned MS
sData = squeeze(sum(cData, 2));
[~, shift] = doAlignment_MinPearson_directed...
    ([HyperROIs.Data{ID}.axisMZ sData]...
    , 0.001, strcmp(HyperROIs.Files.Tags, 'QC'), false);

X = HyperROIs.Data{ID}.axisTM;
Y = HyperROIs.Data{ID}.axisMZ;
[Xq, Yq] = meshgrid(X, Y);
for ii = 1:size(cData, 3)
    [Xn, Yn] = meshgrid(X, Y+shift(ii));
    aData(:, :, ii) = interp2(Xn, Yn, cData(:, :, ii), Xq, Yq);
end
aData(isnan(aData)) = 0;
HRShift = shift;
    
switch method
    case 'unconstrained'
        [Model, FttData] = getMSModels(RefMSSpectra, HyperROIs.FOM.MS{ID}, ...
            aData, HyperROIs.Deconvolution.p);
        
    case 'constrained'
        [Model, FttData] = getMSModels_constrained(RefMSSpectra, HyperROIs.FOM.MS{ID}, ...
            aData, HyperROIs.Deconvolution.p);
        
    case 'fullconstrained'
        [Model, FttData] = getMSModels_fullconstrained(RefMSSpectra, HyperROIs.FOM.MS{ID}, ...
            aData, HyperROIs.Deconvolution.p);
        
end

if isempty(Model), return; end

iD2Rem = sum(Model, 1) == 0;
Model(:, iD2Rem) = [];
FttData(:, iD2Rem) = [];
Model(isnan(Model)) = 0;

Id = sum(Model ~= 0) <= 5;
Model(:, Id) = [];
FttData(:, Id) = [];

% Normalise Model
for ii = 1:size(aData, 3)
    switch method
        case {'unconstrained', 'constrained'}
            DeconvMS(:, :, ii) = squeeze(aData(:, :, ii))'/Model';
        case 'fullconstrained'
            DeconvMS(:, :, ii) = constrainedDeconv(squeeze(aData(:, :, ii))', Model');
    end
end

for ii = 1:size(DeconvMS, 2)
    [Aligned_Data, shift] = doAlignment_MinPearson_directed...
        ([HyperROIs.Data{ID}.axisTM squeeze(DeconvMS(:, ii, :))]...
        , maxShift, strcmp(HyperROIs.Files.Tags, 'QC'));
    DeconvMS(:, ii, :) = Aligned_Data(:, 2:end);
    XY = [ HyperROIs.FOM.MS{ID}(:, 1), Model(:, ii)];
    M = ChrMoment(XY);
    lm = LocalMaxima(XY, 5, max(XY(:,2))*0.05);
    FOM_MS.Id_HyperROIs(end+1, 1) = ID;
    FOM_MS.MS_M1(end+1, 1) = FttData(1, ii);
    FOM_MS.MS_M2(end+1, 1) = FttData(2, ii);
    FOM_MS.isMerged(end+1, 1) = FttData(3, ii);
    FOM_MS.MS_Res_Closest(end+1, 1) = NaN;
    FOM_MS.MS_max_PearsonCorrelation(end+1, 1) = NaN;
    FOM_MS.MS_Model_lm(end+1, 1) = size(lm, 1);
    FOM_MS.MS_Model{end+1, 1} = XY;
    FOM_MS.ProfisBckg(end+1, 1) = false;
    Id = findCloser(M(2), XY(:,1));
    if min(Id, length(XY(:,1))-Id) <= 2
        FOM_MS.ProfisBckg(end, 1) = true;
    end
    
    XY = [Aligned_Data(:,1) sum(Aligned_Data(:,2:end), 2)];
    BCKG = diff(XY(:, 2));
    BCKG = MADPeaks([(1:size(BCKG, 1))' BCKG]);
    BCKG = MADPeaks([XY(BCKG(:,1), 1) XY(BCKG(:,1), 2)]);
    FOM_MS.ProfNoise(end+1, 1) = 3*std(BCKG(:,2));
    FOM_MS.ProfSignal(end+1, 1) = max(XY(:,2));
    
    if FOM_MS.ProfSignal(end, 1)/FOM_MS.ProfNoise(end, 1) < StoN
        FOM_MS.ProfisBckg(end, 1) = true;
    end
    M = ChrMoment([XY(:,1), XY(:, 2)]);
    FOM_MS.CenterTime(end+1, 1) = M(2);
    FOM_MS.shift{end+1, 1} = shift;
end

% Calculate missing FOM
FOM_MS = struct2table(FOM_MS);
if ~isempty(FOM_MS)
    
end

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

if ~isempty(id2deconv)
    for jj = 1:length(id2deconv)
        cData = [HyperROIs.Data{ID}.axisTM squeeze(DeconvMS(:, id2deconv(jj), :))];
        [Model, PeakModels, isMerged] =  fittPMG(cData, false);
        
        Areas = (cData(:, 2:end))'/Model';
        Predict = Areas*Model';
        SSR = sum((cData(:,2:end) - Predict').^2);
        
        for kk = 1:size(Model, 2)
            FOM_Prof.Id(end+1, 1) = length(FOM_Prof.Id)+1;
            FOM_Prof.Id_HyperROIs(end+1, 1) = ID;
            FOM_Prof.Id_MS(end+1, 1) = id2deconv(jj);
            FOM_Prof.MS_M1(end+1, 1) = FOM_MS.MS_M1(id2deconv(jj));
            FOM_Prof.MS_M2(end+1, 1) = FOM_MS.MS_M2(id2deconv(jj));
            FOM_Prof.isMergedMS(end+1, 1) = FOM_MS.isMerged(id2deconv(jj));
            FOM_Prof.MS_Model{end+1, 1} = FOM_MS.MS_Model{id2deconv(jj)};
            
            XY = [HyperROIs.Data{ID}.axisTM Model(:, kk)];
            M = ChrMoment(XY);
            FOM_Prof.Prof_M1(end+1, 1) = M(2);
            FOM_Prof.Prof_M2(end+1, 1) = M(3);
            try
                FOM_Prof.isMergedProf(end+1, 1) = isMerged(kk);
            catch
                disp('pp')
                error(['warning at ' num2str(ID)])
            end
            FOM_Prof.DeconvArea{end+1, 1} = Areas(:, kk);
            FOM_Prof.maxArea(end+1, 1) = max(Areas(:, kk));
            FOM_Prof.DeconvRSSR{end+1, 1} = SSR;
            FOM_Prof.Prof_Model{end+1, 1} = XY;
            FOM_Prof.Prof_Fitting{end+1, 1} = PeakModels{kk};
        end
    end
end

% Calculate missing FOM
FOM_Prof = struct2table(FOM_Prof);
%%%% MISSING FOM TO BE DONE %%%%

