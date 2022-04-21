function [tgtHyperROIs, sumFOM] = optimHyperROIs(obj, dts, timeWdw, MzWdW, QCAnalysis, filter)

% TODO: Change mkHyperROIs to "quick" approach to find limits and outliers
if nargin == 5
    filter = 'none';
    txt4filter = 'None';
    
elseif strcmpi(filter, 'dofilter')
    wx_filter = 3;
    wy_filter = 1;
    nx_filter = 2;
    ny_filter = 2;
    h_filter = sgsdf_2d(-wx_filter:wx_filter, -wy_filter:wy_filter, nx_filter, ny_filter);
    txt4filter = sprintf('sgsdf_2d(%i:%i, %i:%i, %i, %i)', ...
        -wx_filter, wx_filter, -wy_filter, wy_filter, nx_filter, ny_filter);
else
    error('');
end

% TODO: make an object of HyperROIs
% 1. Introduction
tgtHyperROIs.parameters.dts = dts;
tgtHyperROIs.parameters.MzWdW = MzWdW;
tgtHyperROIs.parameters.timeWdw = timeWdw;
tgtHyperROIs.parameters.master = fullfile(obj.Path, obj.Name);
QCFiles = obj.QC.Files;
tgtHyperROIs.Files = table();

switch QCAnalysis
    case 1
        old_FOMs = obj.QC.Method1.FOM;
        
    case 2
        old_FOMs = obj.QC.Method2.FOM;
        
    case 3
        old_FOMs = obj.QC.Method3.FOM;
end

fprintf('\n CREATING HYPERROIS \n\n')

% 2. Setting Limits for ROIs
% 2.1. Open first QC
MF = load(fullfile(QCFiles{1}, 'myFinnee.mat'));
myFinnee = MF.myFinnee;
tgtHyperROIs.Files.ID(1) = 1;
tgtHyperROIs.Files.Name{1} = myFinnee.FileID;
tgtHyperROIs.Files.Tags{1} = 'QC';
tgtHyperROIs.Files.Path{1} = fullfile(QCFiles{1}, 'myFinnee.mat');

fprintf('\t Initialisation with %s\n', myFinnee.FileID)
fprintf('\t Finding ROIs Limits\n')

InfoAxis = myFinnee.Datasets{dts}.AxisX.InfoAxis;
InfoAxis.Loc = 'inAxis';
TimeAxis = Axis(InfoAxis, myFinnee.Datasets{dts}.AxisX.Data);

InfoAxis = myFinnee.Datasets{dts}.AxisY.InfoAxis;
InfoAxis.Loc = 'inAxis';
MzAxis = Axis(InfoAxis, myFinnee.Datasets{dts}.AxisY.Data);

FOMs = table();
FOMs.ID = (1:size(old_FOMs, 1))';
FOMs.centroid_Time = old_FOMs.CtrTime;
FOMs.Time_min = max(old_FOMs.Time_Start - timeWdw, ...
    min(TimeAxis.Data));
FOMs.Time_max = min(old_FOMs.Time_End + timeWdw, ...
    max(TimeAxis.Data));
FOMs.Accurate_Mass = old_FOMs.AccurateMass;

FOMs.mz_min = zeros(size(FOMs, 1), 1);
FOMs.mz_max = zeros(size(FOMs, 1), 1);
for ii = 1:size(FOMs, 1)
    id = findCloser(FOMs.Accurate_Mass(ii), MzAxis.Data);
    IdS = max(1, id - MzWdW);
    FOMs.mz_min(ii) = MzAxis.Data(IdS);
    IdE = min(length(MzAxis.Data), id + MzWdW);
    FOMs.mz_max(ii) = MzAxis.Data(IdE);
end

tgtHyperROIs.parameters4ROIs = FOMs;

% TODO: Modify and merge mkMnROI and mkMnROI_2 => mkMnROI_Finnee &
% mkMnROI_mzmL
fprintf('\t Adding the first ROIs\n')
[ROI, X, Y] = myFinnee.Datasets{dts}.mkMnROI_2...
    (FOMs.mz_min, FOMs.mz_max, FOMs.Time_min, FOMs.Time_max);

for ii = 1:size(FOMs, 1)
    tgtHyperROIs.Data{ii}.axisMZ = X{ii};
    tgtHyperROIs.Data{ii}.axisTM = Y{ii}';
    
    cROI = ROI{ii};
    cROI(isnan(cROI)) = 0;
    if strcmpi(filter, 'dofilter')
        cROI = filter2(h_filter, cROI, 'same');
        cROI(isnan(cROI)) = 0;
    end
    
    M = ChrMoment3D(Y{ii}', X{ii}, cROI);
    tgtHyperROIs.Data{ii}.axisMZ(:, end+1) = mean(cROI, 2);
    lm = LocalMaxima([tgtHyperROIs.Data{ii}.axisMZ(:,1), tgtHyperROIs.Data{ii}.axisMZ(:,end)]...
        , 3, 0.1*max(tgtHyperROIs.Data{ii}.axisMZ(:,end)));
    M.mzLocMax = size(lm, 1);
    M.mzIntPerPeak = size(tgtHyperROIs.Data{ii}.axisMZ, 1)/size(lm, 1);
    
    tgtHyperROIs.Data{ii}.axisTM(:, end+1) = mean(cROI, 1);
    lm = LocalMaxima([tgtHyperROIs.Data{ii}.axisTM(:,1), tgtHyperROIs.Data{ii}.axisTM(:,end)]...
        , 3, 0.1*max(tgtHyperROIs.Data{ii}.axisTM(:,end)));
    M.tmLocMax = size(lm, 1);
    M.tmIntPerPeak = size(tgtHyperROIs.Data{ii}.axisTM, 1)/size(lm, 1);
    
    tgtHyperROIs.Data{ii}.FOM = M;
end


% 2.3. Add remaining QCs

for cF = 2:length(QCFiles)
    MF = load(fullfile(QCFiles{cF}, 'myFinnee.mat'));
    myFinnee = MF.myFinnee;
    fprintf('\n\t Adding %s\n', myFinnee.FileID)
    
    warning off
    tgtHyperROIs.Files.ID(cF) = cF;
    tgtHyperROIs.Files.Name{cF} = myFinnee.FileID;
    tgtHyperROIs.Files.Tags{cF} = 'QC';
    tgtHyperROIs.Files.Path{cF} = fullfile(QCFiles{cF}, 'myFinnee.mat');
    warning on
    
    [ROI, X, Y] = myFinnee.Datasets{dts}.mkMnROI_2...
        (FOMs.mz_min, FOMs.mz_max, FOMs.Time_min, FOMs.Time_max);
    
    fprintf('\t Aligning ROIs\n')
    for ii = 1:size(FOMs, 1)
        [Xq, Yq] = meshgrid(tgtHyperROIs.Data{ii}.axisTM(:,1), tgtHyperROIs.Data{ii}.axisMZ(:,1));
        [Xn, Yn] = meshgrid(Y{ii}', X{ii});
        cROI = interp2(Xn, Yn, ROI{ii}, Xq, Yq);
        cROI(isnan(cROI)) = 0;
        if strcmpi(filter, 'dofilter')
            cROI = filter2(h_filter, cROI, 'same');
            cROI(isnan(cROI)) = 0;
        end
        
        M = ChrMoment3D(tgtHyperROIs.Data{ii}.axisTM(:,1), tgtHyperROIs.Data{ii}.axisMZ(:,1), cROI);
        tgtHyperROIs.Data{ii}.axisMZ(:, end+1) = mean(cROI, 2);
        lm = LocalMaxima([tgtHyperROIs.Data{ii}.axisMZ(:,1), tgtHyperROIs.Data{ii}.axisMZ(:,end)]...
            , 3, 0.1*max(tgtHyperROIs.Data{ii}.axisMZ(:,end)));
        M.mzLocMax = size(lm, 1);
        M.mzIntPerPeak = size(tgtHyperROIs.Data{ii}.axisMZ, 1)/size(lm, 1);
        
        tgtHyperROIs.Data{ii}.axisTM(:, end+1) = mean(cROI, 1);
        lm = LocalMaxima([tgtHyperROIs.Data{ii}.axisTM(:,1), tgtHyperROIs.Data{ii}.axisTM(:,end)]...
            , 3, 0.1*max(tgtHyperROIs.Data{ii}.axisTM(:,end)));
        M.tmLocMax = size(lm, 1);
        M.tmIntPerPeak = size(tgtHyperROIs.Data{ii}.axisTM, 1)/size(lm, 1);
        
        
        tgtHyperROIs.Data{ii}.FOM = [tgtHyperROIs.Data{ii}.FOM; M];
        
    end
end

for ii = 1:size(tgtHyperROIs.Data, 2)
    cFOM = tgtHyperROIs.Data{ii}.FOM;
    sumFOM(ii, 1) = mean(cFOM.M0_X, 'omitnan');
    sumFOM(ii, 2) = std(cFOM.M0_X, [], 'omitnan');
    sumFOM(ii, 3) = mean(cFOM.M1_X, 'omitnan');
    sumFOM(ii, 4) = std(cFOM.M1_X, [], 'omitnan');
    sumFOM(ii, 5) = mean(cFOM.M2_X, 'omitnan');
    sumFOM(ii, 6) = std(cFOM.M2_X, [], 'omitnan');
    sumFOM(ii, 7) = mean(cFOM.M1_Y, 'omitnan');
    sumFOM(ii, 8) = std(cFOM.M1_Y, [], 'omitnan');
    sumFOM(ii, 9) = mean(cFOM.M2_Y, 'omitnan');
    sumFOM(ii, 10) = std(cFOM.M2_Y, [], 'omitnan');
    sumFOM(ii, 11) = mean(cFOM.mzLocMax, 'omitnan');
    sumFOM(ii, 12) = std(cFOM.mzLocMax, [], 'omitnan');
    sumFOM(ii, 13) = mean(cFOM.mzIntPerPeak, 'omitnan');
    sumFOM(ii, 14) = mean(cFOM.tmLocMax, 'omitnan');
    sumFOM(ii, 15) = std(cFOM.tmLocMax, [], 'omitnan');
    sumFOM(ii, 16) = mean(cFOM.tmIntPerPeak, 'omitnan');
end


