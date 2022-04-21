function HyperROIs = mkHyperROIs(obj, dts, timeWdw, MzWdW, QCAnalysis)

%%%%!!!!!!!
%%%% TO BE CHANGED BY MEASURINH mzstart AND mzend DURING CENTROID
%%%% TRANSFORMATIOM
%%%% !!!!!

% 1. Introduction
HyperROIs.parameters.dts = dts;
HyperROIs.parameters.MzWdW = MzWdW;
HyperROIs.parameters.timeWdw = timeWdw;
HyperROIs.parameters.master = fullfile(obj.Path, obj.Name);
QCFiles = obj.QC.Files;
HyperROIs.Files = table();

switch QCAnalysis
    case 1
        HyperROIs.parameters.Initial_FOM = obj.QC.Method1.FOM;
        
    case 2
        HyperROIs.parameters.Initial_FOM = obj.QC.Method2.FOM;
        
    case 3
        HyperROIs.parameters.Initial_FOM = obj.QC.Method3.FOM;
end

[~, name] = fileparts(obj.Name);
nameDir = sprintf('HyperROIs4%s', name);
status = mkdir(fullfile(obj.Path, nameDir));
HyperROIs.parameters.Path = obj.Path; %fullfile(obj.Path, nameDir);
HyperROIs.parameters.name = [nameDir '.mat']; % obj.Path; %fullfile(obj.Path, nameDir);
HyperROIs.ListOfTags = {'QC'};

fprintf('\n CREATING HYPERROIS \n\n')

% 2. Setting Limits for ROIs
% 2.1. Open first QC
MF = load(fullfile(QCFiles{1}, 'myFinnee.mat'));
myFinnee = MF.myFinnee;
HyperROIs.Files.ID(1) = 1;
HyperROIs.Files.Name{1} = myFinnee.FileID;
HyperROIs.Files.Tags{1} = 'QC';
HyperROIs.Files.Path{1} = fullfile(QCFiles{1}, 'myFinnee.mat');

fprintf('\t Initialisation with %s\n', myFinnee.FileID)
fprintf('\t Finding ROIs Limits\n')

InfoAxis = myFinnee.Datasets{dts}.AxisX.InfoAxis;
InfoAxis.Loc = 'inAxis';
HyperROIs.MasterAxis.TimeAxis = Axis(InfoAxis, myFinnee.Datasets{dts}.AxisX.Data);

InfoAxis = myFinnee.Datasets{dts}.AxisY.InfoAxis;
InfoAxis.Loc = 'inAxis';
HyperROIs.MasterAxis.MzAxis = Axis(InfoAxis, myFinnee.Datasets{dts}.AxisY.Data);

FOMs = table();
old_FOMs = HyperROIs.parameters.Initial_FOM;
FOMs.ID = (1:size(old_FOMs, 1))';
FOMs.centroid_Time = old_FOMs.CtrTime;
FOMs.Time_min = max(old_FOMs.Time_Start - timeWdw, ...
    min(HyperROIs.MasterAxis.TimeAxis.Data));
FOMs.Time_max = min(old_FOMs.Time_End + timeWdw, ...
    max(HyperROIs.MasterAxis.TimeAxis.Data));
FOMs.Accurate_Mass = old_FOMs.AccurateMass;

FOMs.mz_min = zeros(size(FOMs, 1), 1);
FOMs.mz_max = zeros(size(FOMs, 1), 1);
for ii = 1:size(FOMs, 1)
    id = findCloser(FOMs.Accurate_Mass(ii), ...
        HyperROIs.MasterAxis.MzAxis.Data);
    IdS = max(1, id - MzWdW);
    FOMs.mz_min(ii) = HyperROIs.MasterAxis.MzAxis.Data(IdS);
    IdE = min(length(HyperROIs.MasterAxis.MzAxis.Data), id + MzWdW);
    FOMs.mz_max(ii) = HyperROIs.MasterAxis.MzAxis.Data(IdE);
    
end

% 2.2. Check for supperposition
fprintf('\t Dealing with overlaping ROIs\n')
myCtr = size(FOMs, 1);
while 1
    Id = find(((FOMs.Time_min <= FOMs.Time_min(myCtr) & FOMs.Time_max >= FOMs.Time_min(myCtr)) | ...
        (FOMs.Time_min <= FOMs.Time_max(myCtr) & FOMs.Time_max >= FOMs.Time_max(myCtr))) & ...
        ((FOMs.mz_min <= FOMs.mz_min(myCtr) & FOMs.mz_max >= FOMs.mz_min(myCtr)) | ...
        (FOMs.mz_min <= FOMs.mz_max(myCtr) & FOMs.mz_max >=  FOMs.mz_max(myCtr))));
    if size(Id, 1) > 1
        FOMs.Time_min(Id(end-1)) = min(FOMs.Time_min(Id(end-1)), ...
            FOMs.Time_min(Id(end)));
        FOMs.mz_min(Id(end-1)) = min(FOMs.mz_min(Id(end-1)), ...
            FOMs.mz_min(Id(end)));
        FOMs.Time_max(Id(end-1)) = max(FOMs.Time_max(Id(end-1)), ...
            FOMs.Time_max(Id(end)));
        FOMs.mz_max(Id(end-1)) = max(FOMs.mz_max(Id(end-1)), ...
            FOMs.mz_max(Id(end)));
        FOMs(Id(end), :) = [];
        myCtr = min(myCtr, size(FOMs, 1));
    else
        myCtr = myCtr - 1
        if myCtr <= 1
            break
        end
    end
end
FOMs.ID = (1:size(FOMs, 1))';
HyperROIs.FOM = FOMs;


fprintf('\t Adding the first ROIs\n')
[ROI, X, Y] = myFinnee.Datasets{dts}.mkMnROI_2( FOMs.mz_min,  ...
    FOMs.mz_max, FOMs.Time_min, FOMs.Time_max);

for ii = 1:size(FOMs)
    HyperROIs.Data{ii}.axisMZ = X{ii};
    HyperROIs.Data{ii}.axisTM = Y{ii}';
    HyperROIs.Data{ii}.size = [size(ROI{ii}) 1];
    HyperROIs.Data{ii}.file = fullfile(obj.Path, nameDir, ['HyperROI#' num2str(ii) '.dat']);
    [fidWriteDat, errmsg]  = fopen(HyperROIs.Data{ii}.file, 'wb');
    cROI = ROI{ii}(:);
    cROI(isnan(cROI)) = 0;
    fwrite(fidWriteDat, ROI{ii}(:), 'uint64');
    fclose(fidWriteDat);
end


% 2.3. Add remaining QCs

for cF = 2:length(QCFiles)
    MF = load(fullfile(QCFiles{cF}, 'myFinnee.mat'));
    myFinnee = MF.myFinnee;
    fprintf('\n\t Adding %s\n', myFinnee.FileID)
    
    warning off
    HyperROIs.Files.ID(cF) = cF;
    HyperROIs.Files.Name{cF} = myFinnee.FileID;
    HyperROIs.Files.Tags{cF} = 'QC';
    HyperROIs.Files.Path{cF} = fullfile(QCFiles{cF}, 'myFinnee.mat');
    warning on
    
    [ROI, X, Y] = myFinnee.Datasets{dts}.mkMnROI_2( FOMs.mz_min,  ...
        FOMs.mz_max, FOMs.Time_min, FOMs.Time_max);
    
    fprintf('\t Aligning ROIs\n')
    for ii = 1:size(FOMs)
        [Xq, Yq] = meshgrid(HyperROIs.Data{ii}.axisTM, HyperROIs.Data{ii}.axisMZ);
        [Xn, Yn] = meshgrid(Y{ii}', X{ii});
        newROI = interp2(Xn, Yn, ROI{ii}, Xq, Yq);
        
        HyperROIs.Data{ii}.size(1, 3) = HyperROIs.Data{ii}.size(1, 3) + 1;
        [fidWriteDat, errmsg]  = fopen(HyperROIs.Data{ii}.file, 'ab');
        cROI = newROI(:);
        cROI(isnan(cROI)) = 0;
        fwrite(fidWriteDat, cROI, 'uint64');
        fclose(fidWriteDat);
    end
end

fprintf('1n\n\t GENERATING MODEL SPECTRA AND PROFILE\n')
for ii = 1:size(HyperROIs.Data, 2)
    [fidReadDat, errmsg]  = fopen(HyperROIs.Data{ii}.file, 'rb');
    cData = reshape(fread(fidReadDat, 'uint64'), HyperROIs.Data{ii}.size);
    fclose(fidReadDat);
    
    allProf = [HyperROIs.Data{ii}.axisTM squeeze(mean(cData, 1))];
    allProf = doAlignment_MinPearson(allProf, 100);
    Prof = [allProf(:,1), sum(allProf(:, 2:end), 2)];
    HyperROIs.FOM.Prof{ii} = Prof;
    lm = LocalMaxima(Prof, 1, max(Prof(:,2))*0.05);
    HyperROIs.FOM.Nlm_Prof(ii) = size(lm, 1);
    HyperROIs.FOM.lm_Prof{ii} = lm;    
    
    allMS = [HyperROIs.Data{ii}.axisMZ, squeeze(mean(cData, 2))];
    allMS = doAlignment_MinPearson(allMS, 2);
    MS = [allMS(:,1), sum(allMS(:, 2:end), 2)];
    HyperROIs.FOM.MS{ii} = MS;
    lm = LocalMaxima(MS, 1, max(MS(:,2))*0.05);
    HyperROIs.FOM.Nlm_MS(ii) = size(lm, 1);
    HyperROIs.FOM.lm_MS{ii} = lm;    
end

% 3. Make model Ion peaks
% 3.1 Normalise
Id = find(HyperROIs.FOM.Nlm_MS == 1);
RefMSSpectra = (-6:0.1:6)';
XYMSSpectra  = [];
for ii = 1:length(Id)
   MSSpectra =  HyperROIs.FOM.MS{Id(ii)};
   M = ChrMoment(MSSpectra);
   XYMSSpectra(end+1, :) = M;
   MSSpectra(:,1) = (MSSpectra(:,1) - M(2))/sqrt(M(3));
   MSSpectra(:,2) = MSSpectra(:,2)/max(MSSpectra(:,2));
   vq = interp1(MSSpectra(:,1), MSSpectra(:,2), RefMSSpectra(:,1));
   RefMSSpectra(:, end+1) = vq;
end

%%3.2 Find and remove outliers
while 1
    aveSpectra = mean( RefMSSpectra(:,2:end), 2, 'omitnan');
    R = corrcoef([aveSpectra,  RefMSSpectra(:,2:end)], 'rows', 'complete');
    io = isoutlier(R(2:end, 1));
    
    if any(io)
        XYMSSpectra(io, :) = [];
        RefMSSpectra(:, [false; io]) = [];
    else
        break
    end
end

HyperROIs.Deconvolution.RefMSSpectra = RefMSSpectra;
HyperROIs.Deconvolution.DataForMS = XYMSSpectra;

file = fullfile(obj.Path, [nameDir '.mat']);
save(file, 'HyperROIs');








    
% 
% %% 1. Initialisation and checks
% p_tm     = {};
% p_mz     = {};
% options  = checkVarargin(varargin{:});
% STOP = false;
% 
% if options.doSmoothing
%     switch options.SmoothingMethod
%         case '2-D Savitzky-Golay'
%             Filter = options.SmoothingFilter;
%             h = sgsdf_2d(Filter{1}, Filter{2}, Filter{3}, Filter{4},0,0);
%     end
% else
%     h = [];
% end
% perc = options.prctile;
% 
% parameters.PathMaster = obj.Path;
% parameters.NameMaster = obj.Name;
% parameters.ListSamples = obj.Samples.Files;
% parameters.IDSamples = obj.Samples.ID;
% parameters.TagsSamples = obj.Samples.Tags;
% parameters.ParametersSamples.dts = obj.Samples.Parameters.dts;
% parameters.ParametersSamples.timeWdw = obj.Samples.Parameters.timeWdw;
% parameters.ParametersSamples.mzWdw = obj.Samples.Parameters.mzWdw;
% parameters.ParametersSamples.Filter = obj.Samples.Parameters.Filter;
% parameters.AnalysisMethod.Files = obj.QC.Files;
% parameters.AnalysisMethod.perameters = obj.QC.parameters;
% parameters.AnalysisMethod.Axis = obj.QC.Axis;
% 
% 
% %% 0- Alignment of all features
% FOM_QC  = obj.Samples.Parameters.FOM;
% FOM_STm = table2array(obj.Samples.FOM.Tm);
% FOM_SMZ = table2array(obj.Samples.FOM.MZ);
% FOM_SAr = table2array(obj.Samples.FOM.Area);
% perc = options.prctile;
% ord_time = options.Time_ord;
% ord_mz = options.mz_ord;
% 
% % 0.1 time dimension
% for ii = 1:length(obj.Samples.Files)
%     IdTft = FOM_SAr(:, 1) >= prctile(FOM_SAr(:, 1), perc);
%     
%     if options.TimeAlign
%         
%         tXY = [FOM_STm(IdTft, ii), FOM_QC.mean_M1(IdTft) - FOM_STm(IdTft, ii)];
%         tXY(tXY(:,1) == 0, :) = [];
%         cp = polyfit(tXY(:,1), tXY(:,2), ord_time);
%         io = isoutlier(tXY(:,2) - polyval(cp, tXY(:,1)));
%         XY = tXY(~io, :);
%         p_tm{ii} = polyfit(XY(:,1), XY(:,2), ord_time);
%         
%         if options.check_TimeAlign
%             [~, cName] = fileparts(fileparts(obj.Samples.Files{ii}));
%             f1 = figure('name', 'Time Alignment', 'units', 'normalized','outerposition',[0 0 1 1]);
%             hPlot = gcf;
%             ax1 = axes('units', 'normalized', 'Position',[0.05 0.25 0.9 0.7]);
%             title(sprintf('Time alignment of %s to averaged QC samples', cName));
%             hold on
%             scatter(ax1, tXY(:,1), tXY(:,2), 'r');
%             Leg{1} = 'outliers';
%             scatter(XY(:,1), XY(:,2), 'k')
%             Leg{2} = 'Experimental data';
%             XYp = (min(tXY(:,1)):(max(tXY(:,1))-min(tXY(:,1)))/100:max(tXY(:,1)))';
%             XYp(:,2) = polyval(p_tm{ii}, XYp(:,1));
%             plot(XYp(:,1), XYp(:,2), 'r');
%             Leg{3} = sprintf('Fitted polynomial (order:%i)', ord_time);
%             legend(Leg)
%             xlabel('Time')
%             ylabel('Delta(Time) (Sample - averaged QC time)')
%             
%             ctrl1 = uicontrol('units', 'normalized', 'Position',[0.05 0.1 0.1 0.05],'Style','radiobutton',...
%                 'String','Use to all');
%             ctrl1.Callback = @selection_ctrl1;
%             ctrl1.Tag = 'Ctr1';
%             
%             ctrl2 = uicontrol('units', 'normalized', 'Position',[0.15 0.1 0.1 0.05],'String','Continue',...
%                 'Callback','uiresume(gcbf)');
%             ctrl2.Callback = @selection_ctrl1;
%             ctrl2.Tag = 'Ctr2';
%             
%             ctrl3 = uicontrol('units', 'normalized', 'Position',[0.25 0.1 0.1 0.05],'String','Stop',...
%                 'Callback','uiresume(gcbf)');
%             ctrl3.Callback = @selection_ctrl1;
%             ctrl3.Tag = 'Ctr3';
%             
%             uiwait(f1);
%             if STOP
%                 warning('Process terminated by user')
%                 return
%             end
%             
%             
%         end
%     else
%         p_tm{ii} = 0;
%     end
% end
% 
% % 0.1 mz dimension
% for ii = 1:length(obj.Samples.Files)
%     IdTft = FOM_SAr(:, 1) >= prctile(FOM_SAr(:, 1), perc);
%     
%     if options.mzAlign
%         
%         tXY = [FOM_SMZ(IdTft, ii), ...
%             (FOM_QC.mean_AccMass(IdTft) - FOM_SMZ(IdTft, ii))./FOM_SMZ(IdTft, ii)];
%         tXY(tXY(:,1) == 0, :) = [];
%         cp = polyfit(tXY(:,1), tXY(:,2), ord_mz);
%         io = isoutlier(tXY(:,2) - polyval(cp, tXY(:,1)));
%         XY = tXY(~io, :);
%         p_mz{ii} = polyfit(XY(:,1), XY(:,2), ord_mz);
%         
%         if options.check_mzAlign
%             [~, cName] = fileparts(fileparts(obj.Samples.Files{ii}));
%             f1 = figure('name', 'mass Alignment', 'units', 'normalized','outerposition',[0 0 1 1]);
%             hPlot = gcf;
%             ax1 = axes('units', 'normalized', 'Position',[0.05 0.25 0.9 0.7]);
%             title(sprintf('mz alignment of %s to averaged QC samples', cName));
%             hold on
%             scatter(ax1, tXY(:,1), tXY(:,2), 'r');
%             Leg{1} = 'outliers';
%             scatter(XY(:,1), XY(:,2), 'k')
%             Leg{2} = 'Experimental data';
%             XYp = (min(XY(:,1)):(max(XY(:,1))-min(XY(:,1)))/100:max(XY(:,1)))';
%             XYp(:,2) = polyval(p_mz{ii}, XYp(:,1));
%             plot(XYp(:,1), XYp(:,2), 'r');
%             Leg{3} = sprintf('Fitted polynomial (order:%i)', ord_mz);
%             legend(Leg)
%             xlabel('mz (sample)')
%             ylabel('error(mz) ((Sample - averaged QC mz)/sample mz)')
%             
%             ctrl1 = uicontrol('units', 'normalized', 'Position',[0.05 0.1 0.1 0.05],'Style','radiobutton',...
%                 'String','Use to all');
%             ctrl1.Callback = @selection_ctrl2;
%             ctrl1.Tag = 'Ctr1';
%             
%             ctrl2 = uicontrol('units', 'normalized', 'Position',[0.15 0.1 0.1 0.05],'String','Continue',...
%                 'Callback','uiresume(gcbf)');
%             ctrl2.Callback = @selection_ctrl2;
%             ctrl2.Tag = 'Ctr2';
%             
%             ctrl3 = uicontrol('units', 'normalized', 'Position',[0.25 0.1 0.1 0.05],'String','Stop',...
%                 'Callback','uiresume(gcbf)');
%             ctrl3.Callback = @selection_ctrl2;
%             ctrl3.Tag = 'Ctr3';
%             
%             uiwait(f1);
%             if STOP
%                 warning('Process terminated by user')
%                 return
%             end
%             
%             
%         end
%     else
%         p_mz{ii} = 0;
%     end
% end
% 
% 
% 
% 
% 
% 
% %% 1- Load ROIs and correct tm & MZ
% ROIs.Data  = {};
% ROIs.axeMz = {};
% ROIs.axeTm = {};
% 
% for ii = 1:length(obj.Samples.Files)
%     path = fileparts(obj.Samples.Files{ii});
%     ROI = load(fullfile(path, options.ROIsname));
%     
%     for jj = 1:length(ROI.tgtROIs.ROI)
%         cROI = ROI.tgtROIs.ROI{jj};
%         
%         if ii == 1
%             ROIs.axeMz{jj, 1} = cROI.AxisMZ.Data;
%             ROIs.axeTm{jj, 1} = cROI.AxisTm.Data;
%             ROIs.Data{jj, 1}  = NaN(length(ROIs.axeMz{jj}), length(ROIs.axeTm{jj}), length(obj.Samples.Files));
%         end
%         [Xq, Yq] = meshgrid(ROIs.axeTm{jj}, ROIs.axeMz{jj});
%         [X, Y] = meshgrid(cROI.AxisTm.Data + polyval(p_tm{ii}, cROI.AxisTm.Data), ...
%             cROI.AxisMZ.Data + polyval(p_mz{ii}, cROI.AxisMZ.Data).* cROI.AxisMZ.Data);
%         cData = imfilter(cROI.StoredData, h);
%         ROIs.Data{jj}(:, :, ii) = interp2(X, Y, cData, Xq, Yq);
%         
%     end
% end
% 
% ROIs.M0_tm = [];
% ROIs.M1_tm = [];
% ROIs.M2_tm = [];
% ROIs.M3_tm = [];
% ROIs.M0_mz = [];
% ROIs.M1_mz = [];
% ROIs.M2_mz = [];
% ROIs.M3_mz = [];
% 
% for ii = 1:size(ROIs.Data, 1)
%     Data_3D = ROIs.Data{ii};
%     axeMz      =  ROIs.axeMz{ii};
%     axeMz(:,2) = squeeze(sum(sum(Data_3D, 2, 'omitnan'), 3, 'omitnan'));
%     axeMz(:,3) = squeeze(std(sum(Data_3D, 2, 'omitnan'), [], 3, 'omitnan'));
%     ROIs.axeMz{ii} = axeMz;
%     M = ChrMoment(axeMz(:, 1:2));
%     ROIs.M0_mz(end+1, 1) = M(1);
%     ROIs.M1_mz(end+1, 1) = M(2);
%     ROIs.M2_mz(end+1, 1) = M(3);
%     ROIs.M3_mz(end+1, 1) = M(4);
%     
%     axeTm      =  ROIs.axeTm{ii};
%     axeTm(:,2) = squeeze(sum(sum(Data_3D, 1, 'omitnan'), 3, 'omitnan'));
%     axeTm(:,3) = squeeze(std(sum(Data_3D, 1, 'omitnan'), [], 3, 'omitnan'));
%     ROIs.axeTm{ii} = axeTm;
%     M = ChrMoment(axeTm(:, 1:2));
%     ROIs.M0_tm(end+1, 1) = M(1);
%     ROIs.M1_tm(end+1, 1) = M(2);
%     ROIs.M2_tm(end+1, 1) = M(3);
%     ROIs.M3_tm(end+1, 1) = M(4);
%     
% end
% HyperROIs.HyperROIs = ROIs;
% HyperROIs.Parameters = parameters;
% HyperROIs.Options = options;
% HyperROIs.Alignment.p_mz = p_mz;
% HyperROIs.Alignment.p_tm = p_tm;
% 
%     function options = checkVarargin(varargin)
%         % CHECKVARARGIN is used to check the input paramters
%         % and create the options parameter.
%         
%         % 1- Default parameters
%         
%         options.ROIsname  = 'ROIs4Quant.mat';
%         options.TimeAlign = true;
%         options.Time_ord  = 1;
%         options.mzAlign = true;
%         options.mz_ord  = 1;
%         options.check_TimeAlign = true;
%         options.check_mzAlign = true;
%         options.doSmoothing   = true;
%         options.SmoothingMethod = '2-D Savitzky-Golay';
%         options.SmoothingFilter = {-3:3,-3:3,2,2};
%         options.prctile = 25;
%         
%         % 3- Decipher varargin
%         input = @(x) find(strcmpi(varargin,x),1);
%         tgtIx = input('ROIsname');
%         if ~isempty(tgtIx)
%             options.ROIsname = varargin{tgtIx +1};
%         end
%         
%         tgtIx = input('Time_ord');
%         if ~isempty(tgtIx)
%             n = varargin{tgtIx +1};
%             if n == 0
%                 options.TimeAlign = false;
%             else
%                 options.TimeAlign = n;
%             end
%         end
%         
%         tgtIx = input('prctile');
%         if ~isempty(tgtIx)
%             options.prctile = varargin{tgtIx +1};
%         end
%         
%         tgtIx = input('mz_ord');
%         if ~isempty(tgtIx)
%             n = varargin{tgtIx +1};
%             if n == 0
%                 options.mzAlign = false;
%             else
%                 options.mz_ord = n;
%             end
%         end
%     end
% 
%     function selection_ctrl1(src, event)
%         switch src.Tag
%             case 'Ctr1'
%                 if src.Value == 1
%                     options.check_TimeAlign = false;
%                     
%                 else
%                     options.check_TimeAlign = true;
%                 end
%                 
%             case 'Ctr2'
%                 close(f1)
%                 
%             case'Ctr3'
%                 STOP = true;
%                 close(f1);
%         end
%     end
% 
%     function selection_ctrl2(src, event)
%         switch src.Tag
%             case 'Ctr1'
%                 if src.Value == 1
%                     options.check_mzAlign = false;
%                     
%                 else
%                     options.check_mzAlign = true;
%                 end
%                 
%             case 'Ctr2'
%                 close(f1)
%                 
%             case'Ctr3'
%                 STOP = true;
%                 close(f1);
%         end
%     end
% end
