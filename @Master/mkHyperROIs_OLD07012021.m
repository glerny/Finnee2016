function HyperROIs = mkHyperROIs(obj, dts, timeWdw, MzWdW, QCAnalysis, filter)
% NOTE: filter should be extended
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
% NOTE: The limits are not check anymore at this point. In the first pass
% HROIs are constructed using the whole data from the QC analysis. Limits
% and normalisation will be obtain using the HROIs

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

HyperROIs.FOM = FOMs;

% TODO: Modify and merge mkMnROI and mkMnROI_2 => mkMnROI_Finnee &
% mkMnROI_mzmL

fprintf('\t Adding the first ROIs\n')
[ROI, X, Y] = myFinnee.Datasets{dts}.mkMnROI_2( FOMs.mz_min,  ...
    FOMs.mz_max, FOMs.Time_min, FOMs.Time_max);

for ii = 1:size(FOMs, 1)
    
    cROI = ROI{ii};
    cROI(isnan(cROI)) = 0;
    if strcmpi(filter, 'dofilter')
        cROI = filter2(h_filter, cROI, 'same');
        cROI(isnan(cROI)) = 0;
    end
    
    %% Create HD5 files
    file = fullfile(obj.Path, nameDir, ['HyperROI#' num2str(ii) '.h5']);
    
    % Check if file exist, if yes delete
    % DOTO: give choice
    if isfile(file)
        delete(file)
    end
    
    % Axis Time
    loc = '/ROIs/AxisTime';
    h5create(file, loc, size(Y{ii}'), 'Datatype', 'double')
    h5write(file, loc, Y{ii}')
    h5writeatt(file, loc, 'Label', myFinnee.Datasets{1}.AxisX.Label)
    h5writeatt(file, loc, 'Unit', myFinnee.Datasets{1}.AxisX.Unit)
    h5writeatt(file, loc, 'DateOfCreation', datestr(datetime))
    
    % Axis mz
    loc = '/ROIs/AxisMz';
    h5create(file, loc, size(X{ii}), 'Datatype', 'double')
    h5write(file, loc,X{ii})
    h5writeatt(file, loc, 'Label', myFinnee.Datasets{1}.AxisY.Label)
    h5writeatt(file ,loc, 'Unit', myFinnee.Datasets{1}.AxisY.Unit)
    h5writeatt(file, loc, 'DateOfCreation', datestr(datetime))
    
    % First data set
    loc = sprintf('/ROIs/ROI4%s', myFinnee.FileID);
    h5create(file, loc, size(cROI), 'Datatype', 'double')
    h5write(file, loc, cROI)
    h5writeatt(file, loc, 'Normalisation@mz', 1)
    h5writeatt(file, loc, 'Normalisation@time', 1)
    h5writeatt(file, loc, 'TAG', 'QC')
    h5writeatt(file, loc, 'Filter', txt4filter)
    h5writeatt(file, loc, 'DateOfCreation', datestr(datetime))
    M = ChrMoment3D(Y{ii}', X{ii}, cROI);
    h5writeatt(file, loc, 'Label4FOM', ['Time_M0, Time_M1, Time_M2, Time_M3, ', ...
        'MZ_M0, MZ_M1, MZ_M2, MZ_M3'])
    h5writeatt(file, loc, 'FOM', table2array(M))
    
    % Remaining data
    loc = '/';
    h5writeatt(file, loc, 'DateOfCreation', datestr(datetime))
    h5writeatt(file, loc, 'MzInterval', [FOMs.mz_min(ii) FOMs.mz_max(ii)])
    h5writeatt(file, loc, 'TimeInterval', [FOMs.Time_min(ii) FOMs.Time_max(ii)])
    h5writeatt(file, loc, 'MasterObj', fullfile(obj.Path, obj.Name))
    
    HyperROIs.Data{ii}.file = file;
    HyperROIs.Data{ii}.FOM = M;
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
    for ii = 1:size(FOMs, 1)
        file = fullfile(obj.Path, nameDir, ['HyperROI#' num2str(ii) '.h5']);
        [Xq, Yq] = meshgrid(HyperROIs.Data{ii}.axisTM, HyperROIs.Data{ii}.axisMZ);
        [Xn, Yn] = meshgrid(Y{ii}', X{ii});
        newROI = interp2(Xn, Yn, ROI{ii}, Xq, Yq);
        newROI(isnan(newROI)) = 0;
        if strcmpi(filter, 'dofilter')
            newROI = filter2(h_filter, newROI, 'same');
            newROI(isnan(newROI)) = 0;
        end
        
         % Write data set
         loc = sprintf('/ROIs/ROI4%s', myFinnee.FileID);
         h5create(file, loc, size(cROI), 'Datatype', 'double')
         h5write(file, loc, cROI)
         h5writeatt(file, loc, 'Normalisation@mz', 1)
         h5writeatt(file, loc, 'Normalisation@time', 1)
         h5writeatt(file, loc, 'TAG', 'QC')
         h5writeatt(file, loc, 'Filter', txt4filter)
         h5writeatt(file, loc, 'DateOfCreation', datestr(datetime))
         M = ChrMoment3D(Y{ii}', X{ii}, cROI);
         h5writeatt(file, loc, 'Label4FOM', ['Time_M0, Time_M1, Time_M2, Time_M3, ', ...
             'MZ_M0, MZ_M1, MZ_M2, MZ_M3'])
         h5writeatt(file, loc, 'FOM', table2array(M))
         HyperROIs.Data{ii}.FOM = [HyperROIs.Data{ii}.FOM; M];
    end
end

file = fullfile(obj.Path, [nameDir '.mat']);
save(file, 'HyperROIs');