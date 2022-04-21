function [HyperROIs, links2dat] = mkHyperROIs(obj, dts, Limits)
% NOTE: filter should be extended

%% 2D SG filter
wtm_filter = 3;
wmz_filter = 2;
ntm_filter = 1;
nmz_filter = 1;
h_filter = sgsdf_2d...
    (-wtm_filter:wtm_filter, -wmz_filter:wmz_filter, ntm_filter, nmz_filter);

% TODO: make an object of HyperROIs
% 1. Introduction
HyperROIs.parameters.dts = dts;
HyperROIs.parameters.master = fullfile(obj.Path, obj.Name);
HyperROIs.parameters.limits = Limits;
QCFiles = obj.QC.Files;
HyperROIs.Files = table();
[~, name] = fileparts(obj.Name);

%TODO: allow choosing the name of the folder of HROIs
nameDir = sprintf('HyperROIs4%s', name);

% TODO: condition in error of mkdir
try
    status = mkdir(fullfile(obj.Path, nameDir));
catch EM
    status
    retrow(EM)
end

HyperROIs.parameters.Path = obj.Path; %fullfile(obj.Path, nameDir);
HyperROIs.parameters.name = [nameDir '.mat']; % obj.Path; %fullfile(obj.Path, nameDir);
HyperROIs.ListOfTags = {'QC'};

fprintf('\n CREATING HYPERROIS \n\n')

% 2. Setting Limits for ROIs
% 2.1. Open first QC
MF = load(fullfile(QCFiles{1}, 'myFinnee.mat'));
myFinnee = MF.myFinnee; clear MF;
HyperROIs.Files.ID(1) = 1;
HyperROIs.Files.Name{1} = myFinnee.FileID;
HyperROIs.Files.Tags{1} = 'QC';
HyperROIs.Files.Path{1} = fullfile(QCFiles{1}, 'myFinnee.mat');
HyperROIs.myPath     = fullfile(fullfile(obj.Path, nameDir), 'HyperROIs.mat');
HyperROIs.links2dat  = fullfile(fullfile(obj.Path, nameDir), 'links2dat.mat');
HyperROIs.Path2Dat   = fullfile(fullfile(obj.Path, nameDir), 'Data4HyperROIs');

fprintf('\t Initialisation with %s\n', myFinnee.FileID)
fprintf('\t Finding ROIs Limits\n')

InfoAxis = myFinnee.Datasets{dts}.AxisX.InfoAxis;
InfoAxis.Loc = 'inAxis';
HyperROIs.MasterAxis.TimeAxis = Axis(InfoAxis, myFinnee.Datasets{dts}.AxisX.Data);

InfoAxis = myFinnee.Datasets{dts}.AxisY.InfoAxis;
InfoAxis.Loc = 'inAxis';
HyperROIs.MasterAxis.MzAxis = Axis(InfoAxis, myFinnee.Datasets{dts}.AxisY.Data);

% NOTE: The limits are not check anymore at this point. In the first pass
% HROIs are constructed using the whole data from the QC analysis. Limits
% and normalisation will be obtain using the HROIs


% TODO: Modify and merge mkMnROI and mkMnROI_2 => mkMnROI_Finnee &
% mkMnROI_mzmL

% TODO: condition in error of mkdir
try
    status = mkdir(HyperROIs.Path2Dat);
catch EM
    status
    retrow(EM)
end


fprintf('\t Adding the first ROIs\n')
% TODO: Change .MkMnROI_2 to .MkMnROI_det for improved speed.
[ROI, X, Y] = myFinnee.Datasets{dts}.mkMnROI_2...
    (Limits.mz_min, Limits.mz_max, Limits.Time_min, Limits.Time_max);

path = obj.Path;
Xlabel = myFinnee.Datasets{1}.AxisX.Label;
Xunit = myFinnee.Datasets{1}.AxisX.Unit;
Ylabel = myFinnee.Datasets{1}.AxisX.Label;
Yunit = myFinnee.Datasets{1}.AxisX.Unit;
FileID = myFinnee.FileID;


for ii = 1:size(Limits, 1)
    file2write = fullfile(HyperROIs.Path2Dat, ['HROI#' num2str(Limits.ID(ii))]);
    Datafiles{ii, 1} = file2write;
    [fidWriteDat, errmsg]  = fopen(file2write, 'wb');
    cROI = ROI{ii};
    cROI(isnan(cROI)) = 0;
    fROI = filter2(h_filter, cROI, 'same');
    fROI(isnan(fROI)) = 0;
    
    Moment3D{ii, 1} = ChrMoment3D( Y{ii}', X{ii}, fROI);
    ID(ii, 1)            = Limits.ID(ii);
    mz_Interval(ii, :)   = [Limits.mz_min(ii), Limits.mz_max(ii)];
    time_Interval(ii, :) = [Limits.Time_min(ii), Limits.Time_max(ii)];
    Size(ii, :)          = size(fROI);
    format{ii, 1}        = 'double';
    fwrite(fidWriteDat, Y{ii}', 'double');
    fwrite(fidWriteDat, X{ii}, 'double');
    fwrite(fidWriteDat, fROI(:), 'double');
    fclose(fidWriteDat);
end

% 2.3. Add remaining QCs
for cF = 2:length(QCFiles)
    MF = load(fullfile(QCFiles{cF}, 'myFinnee.mat'));
    myFinnee = MF.myFinnee; clear MF;
    fprintf('\n\t Adding %s\n', myFinnee.FileID)
    
    HyperROIs.Files.ID(cF) = cF;
    HyperROIs.Files.Name{cF} = myFinnee.FileID;
    HyperROIs.Files.Tags{cF} = 'QC';
    HyperROIs.Files.Path{cF} = fullfile(QCFiles{cF}, 'myFinnee.mat');
    clear ROI X Y
    [ROI, X, Y] = myFinnee.Datasets{dts}.mkMnROI_2...
        (Limits.mz_min, Limits.mz_max, Limits.Time_min, Limits.Time_max);
    
    fprintf('\t Aligning ROIs\n')
    for ii = 1:size(Limits, 1)
        file2add = Datafiles{ii};
        [fidWriteDat, errmsg]  = fopen(file2add, 'a+');
        while fidWriteDat == -1
            [fidWriteDat, errmsg]  = fopen(file2add, 'a+');
        end
        
        
        fseek(fidWriteDat,  0, 'bof');
        data = fread(fidWriteDat,inf, format{ii});
        X_Tm = data(1:Size(ii,2));
        Y_mz = data(Size(ii,2)+1:Size(ii,2)+Size(ii,1));
        
        [Xq, Yq] = meshgrid(X_Tm, Y_mz);
        [Xn, Yn] = meshgrid(Y{ii}', X{ii});
        newROI = interp2(Xn, Yn, ROI{ii}, Xq, Yq);
        newROI(isnan(newROI)) = 0;
        fROI = filter2(h_filter, newROI, 'same');
        fROI(isnan(fROI)) = 0;
        
        Moment3D{ii, 1} = [Moment3D{ii, 1}; ChrMoment3D(X_Tm, Y_mz, fROI)];
        % Write data set
        fwrite(fidWriteDat, fROI, 'double');
        fclose(fidWriteDat);
    end
end

links2dat = table(ID, mz_Interval, time_Interval, Moment3D, Size, format, Datafiles);
save(HyperROIs.links2dat, 'links2dat')
save(HyperROIs.myPath, 'HyperROIs')