function [HyperROIs, DeconvolvedPeaks, BackgroundSignals, myToc]  = SampleAnalysis_1(HyperROIs, FilterToKeep)

tic
path = fileparts(HyperROIs.myPath);
HyperROIs.Deconvolution.link2PeakList = fullfile(path, 'links2PeakList.mat');
HyperROIs.Deconvolution.link2PeakListDat = fullfile(path, 'peakdata.dat');

[fidWriteDat, errmsg]  = fopen(HyperROIs.Deconvolution.link2PeakListDat, 'wb');

% Check if peak with QCs
Options.maxPeaks = 20;
Options.Function = 'PMG1';
Options.LoopMe = 25;
Options.RecursiveLoop = 0.95;
Options.InitialFactor = [0.8 0.4];
Options.MinResolution = 0.8;
Options.Penalisation = true;
Options.PenalisationWeight = 1.2;
Options.Constrained.SharedParameters = "None";
Options.PointsPerPeaks = [10 50];
Options.MinMax = 0.05;
Options.Robust = false;
Options.MaxLeeway  = 0.2;
thresI = 10;


% Check in decon profiles are likely to
TBL = load(HyperROIs.Deconvolution.link2DeconvTbl);
links2deconv = TBL.links2deconv; clear TBL;
[fidRead, errmsg]  = fopen(HyperROIs.Deconvolution.link2DeconvDat, 'r');
links2deconv = sortrows(links2deconv,'ID','ascend');

if nargin == 1
    FilterToKeep = (1:size(links2deconv, 1))';
end
count = 0;
count_BS = 0;

for II = 1:size(FilterToKeep, 1)
    ii = FilterToKeep(II);
    
    Pos = links2deconv.Position(ii, :);
    Size = links2deconv.SizeDeconv(ii, :);
    Format = links2deconv.Format{ii};
    fseek(fidRead,  Pos(1), 'bof');
    data = fread(fidRead,...
        (Pos(2)- Pos(1))/8, Format);
    data = reshape(data, Size);
    Xdata = data(:, 2:end);
    Xdata(Xdata < thresI) = 0;
    data(:, 2:end) = Xdata; clear Xdata;
    
    [Aligned_Data, shift, leeway] = doAlignment_MinPearson...
        (data,  strcmp(HyperROIs.Files.Tags, 'QC'), Options.MaxLeeway, true);
    
    Aligned_Data(isnan(Aligned_Data)) = 0;
    
    [Model, FittedChannels, Stats, FittedProfile] = IterativeMethod1_fminsearch(Aligned_Data(:,1), Aligned_Data(:,2:end), Options);
    
    
    if Stats.NbrPeaks == 0
        count_BS = count_BS + 1;
        BackgroundSignals.ID(count_BS, 1) = count_BS;
        BackgroundSignals.ID2Deconv(count_BS, 1) = links2deconv.ID(ii);
        BackgroundSignals.ID2HROI(count_BS, 1) = links2deconv.ID2HROIs(ii);
        BackgroundSignals.Amplitude(count_BS, 1) = links2deconv.Amplitude(ii);
        BackgroundSignals.AlignmentLeeway(count_BS, 1) = leeway;
        
        continue
    end
    
    FittedProfile(isnan(FittedProfile)) = 0;
    
    % Merge one and 2 for baseline corr.
    Baseline = squeeze(FittedChannels(:, 1, :)) + squeeze(FittedChannels(:, 2, :));
    FittedChannels(:, 1:2, :) = [];
    
    FittedChannels(FittedChannels < thresI) = 0; FittedChannels = round(FittedChannels);
    % Measured goodness of fit
    clear Deconvoluted_Profiles
    for jj = 1:size(Aligned_Data, 2)-1
        filter = squeeze(FittedChannels(:, :, jj));
        filter(filter(:,1) < thresI, 1) = thresI;
        filter = filter./sum(filter, 2);
        Deconvoluted_Profiles(:, :, jj) = Aligned_Data(:, jj+1).*filter - Baseline(:, jj);
    end
    
    Deconvoluted_Profiles(isnan(Deconvoluted_Profiles)) = 0;
    Deconvoluted_Profiles(Deconvoluted_Profiles < thresI) = 0;
    
    clear Cor2Profile Cor2Model
    for jj = 1:size(Deconvoluted_Profiles, 2)
        Cor2Model(jj) = ...
            1-pdist2(FittedProfile(:, jj+2)', mean(Deconvoluted_Profiles(:, jj, :), 3)', 'correlation');
    end
    
    
    for jj = 1:Stats.NbrPeaks
        count = count+1;
        DeconvolvedPeaks.ID(count, 1) = count;
        DeconvolvedPeaks.ID2Deconv(count, 1) = links2deconv.ID(ii);
        DeconvolvedPeaks.ID2HROI(count, 1) = links2deconv.ID2HROIs(ii);
        DeconvolvedPeaks.relRMSE_mz(count, 1) = links2deconv.relRMSE(ii);
        DeconvolvedPeaks.Amplitude(count, 1) = links2deconv.Amplitude(ii);
        DeconvolvedPeaks.AlignmentMethod{count, 1} = 'doAlignment_MinPearson';
        DeconvolvedPeaks.AlignmentShift{count, 1} = shift;
        DeconvolvedPeaks.AlignmentLeeway(count, 1) = leeway;
        DeconvolvedPeaks.Method{count, 1} = 'MultiVariateFittingMethod_2';
        DeconvolvedPeaks.FittingFunction{count, 1} = Model.Peaks.Function;
        DeconvolvedPeaks.Model{count, 1} = Model;
        DeconvolvedPeaks.nbrPeaks(count, 1) = Stats.NbrPeaks;
        DeconvolvedPeaks.RMSE(count, 1) = Stats.RMSE;
        DeconvolvedPeaks.relRMSE(count, 1) = Stats.RMSE/links2deconv.Amplitude(ii);
        DeconvolvedPeaks.meanBasLvl(count, 1:2) = ...
            [mean(Baseline, 'all'), std(Baseline, [], 'all')];
        DeconvolvedPeaks.Options{count, 1} = Options;
        DeconvolvedPeaks.PeakID(count, 1) = jj;
        DeconvolvedPeaks.MZ(count, :) = links2deconv.MZModel(ii);
        DeconvolvedPeaks.FittingParameters(count, :) = Model.Peaks.FittingParameters(:, jj)';
        DeconvolvedPeaks.Intensities(count, :) = Model.Peaks.FittingIntensities(:, jj)';
        lx = Model.Peaks.FittingIntensities(strcmp(HyperROIs.Files.Tags, 'QC'), jj);
        DeconvolvedPeaks.QC_RSD(count, 1) = std(lx)/mean(lx)*100;
        DeconvolvedPeaks.MeanQCIntensities(count, 1) = mean(lx);
        
        % Trim Peaks;
        XY = max([squeeze(FittedChannels(:,jj, :)),  squeeze(Deconvoluted_Profiles(:,jj, :))],[], 2);
        XY(:, 2) = circshift(XY(:,1), 1); XY(1, 2) = 0;
        XY(:, 3) = circshift(XY(:,1), -1); XY(end, 3) = 0;
        Id2rec = sum(XY < 10, 2) < 3;
        
        data2write = [Aligned_Data(Id2rec,1), FittedProfile(Id2rec,jj+2)];
        DeconvolvedPeaks.ModelProfile(count, 1:2) = size(data2write);
        DeconvolvedPeaks.ModelProfile(count, 3) = ftell(fidWriteDat);
        fwrite(fidWriteDat, data2write(:), 'double');
        DeconvolvedPeaks.ModelProfile(count, 4) = ftell(fidWriteDat);
        
        data2write = [Aligned_Data(Id2rec,1), squeeze(Deconvoluted_Profiles(Id2rec,jj, :))];
        DeconvolvedPeaks.DeconvolvedProfile(count, 1:2) = size(data2write);
        DeconvolvedPeaks.DeconvolvedProfile(count, 3) = ftell(fidWriteDat);
        fwrite(fidWriteDat, data2write(:), 'double');
        DeconvolvedPeaks.DeconvolvedProfile(count, 4) = ftell(fidWriteDat);
        
        DeconvolvedPeaks.Deconvolution_Pearson2Model(count, 1) = Cor2Model(jj);
    end
end

fclose(fidRead);
fclose(fidWriteDat);
DeconvolvedPeaks = struct2table(DeconvolvedPeaks);
DeconvolvedPeaks = sortrows(DeconvolvedPeaks,'ID','ascend');
save(HyperROIs.Deconvolution.link2PeakList, 'DeconvolvedPeaks')
save(HyperROIs.myPath, 'HyperROIs');
myToc = toc;

end
