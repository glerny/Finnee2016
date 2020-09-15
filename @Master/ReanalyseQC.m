function obj = ReanalyseQC(obj, dts1, timeWdw, mzWdw, varargin)

%% 1- Parameters and calculating the ROIs(to be changed and add to varadin)
name      = 'ROIs4QC';
name2     = 'ROIs4QC_opti';
overwrite = false;
Mw        = 100000;

obj.mkROIs(dts1, name, mzWdw, timeWdw, overwrite)


%% 2- Getting the FOM based on results of method 1 -> Method2
% 2.1. Measuring FOM
% REDO IT ALL DTS1/DTS2 = 4 (bsl corrected, noiser removed)

fprintf('\n\nSTEP 1: Getting limits\n\n')
allFOM               = table();
allFOM.IDFile        = zeros(0);
allFOM.IDFeature     = zeros(0);
allFOM.M0            = zeros(0);
allFOM.M1            = zeros(0);
allFOM.M2            = zeros(0);
allFOM.M3            = zeros(0);
allFOM.TimeAtPeakMax = zeros(0);
allFOM.IntAtPeakMax  = zeros(0);
allFOM.nbrPts        = zeros(0);
allFOM.Lm1           = zeros(0);
allFOM.Lm2           = zeros(0);
allFOM.mean_mz       = zeros(0);
allFOM.std_mz        = zeros(0);
allFOM.accuMass      = zeros(0);


for ii = 1:length(obj.QC.Files)
    
    fprintf('\t processing %s\n', obj.QC.Files{ii})
    TR{ii} = load(fullfile(obj.QC.Files{ii}, [name '.mat']));
    
    h = waitbar(0,'Processing ROIs');
    for jj = 1:length(TR{ii}.tgtROIs.ROI)
        waitbar(jj/length(TR{ii}.tgtROIs.ROI))
        cROI = TR{ii}.tgtROIs.ROI{jj};
        ChromMom =  cROI.getFOM( 3, 11, 5, []);
        
        if isempty(ChromMom)
            warning off
            allFOM.IDFile(end+1)  = ii;
            allFOM.IDFeature(end) = jj;
            warning on
            
        else
            [~, it] = min(((ChromMom.TimeAtPeakMax - cROI.TgtTm)/100).^2 +...
                (ChromMom.accuMass - cROI.TgtMz).^2);
            head = table(ii, jj, 'VariableNames', {'IDFile', 'IDFeature'});
            allFOM = [allFOM; [head ChromMom(it, :)]];
        end
        
    end
    try close(h); catch, end %#ok<CTCH>
    
end

%2.2. Getting profiles

mrgFOM.IDFeature    = [];
mrgFOM.nbrDetec     = [];
mrgFOM.mean_tMax    = [];
mrgFOM.std_tMax     = [];
mrgFOM.mean_M1      = [];
mrgFOM.std_M1       = [];
mrgFOM.mean_AccMass = [];
mrgFOM.str_AccMass  = [];
mrgFOM.mean_Area    = [];
mrgFOM.std_Area     = [];
mrgFOM.RSD_Area     = [];
mrgFOM.Lm1          = [];
mrgFOM.Lm2          = [];
allPrf = {};


for ii = 1:length(TR{1}.tgtROIs.ROI)
    cFOM = allFOM(allFOM.IDFeature == ii, :);
    
    Profiles = obj.QC.Axis.AxisX.Data;
    for jj = 1:length(TR)
        if cFOM.M0(jj) == 0
            Profiles(:, end+1) = 0;
            
        else
            cROI = TR{jj}.tgtROIs.ROI{ii};
            
            try
            XYZ = cROI.AxisTm.Data;
            XYZ(:,2) = trapz(cROI.AxisMZ.Data, cROI.StoredData);
            XYZ(isnan(XYZ)) = 0;
            
            is = findCloser(cFOM.Lm1(jj), XYZ(:,1));
            ie = findCloser(cFOM.Lm2(jj), XYZ(:,1));
            
            vq = interp1(XYZ(is:ie, 1), XYZ(is:ie, 2), Profiles(:,1));
            vq(isnan(vq)) = 0;
            Profiles(:, end+1) = vq;
            catch
                disp('wth')
            end
        end
    end
    
    Profiles(1,2:end) = zeros(size(Profiles(1,2:end)));
    Profiles(end,2:end) = zeros(size(Profiles(1,2:end)));
    PrfMinus = circshift(Profiles(:, 2:end), -1);
    PrfPlus  = circshift(Profiles(:, 2:end), 1);
    Profiles(sum([PrfMinus, Profiles(:, 2:end), PrfPlus] == 0, 2) == 3*jj, :) = [];
    allPrf{end+1} = [Profiles(:,1), mean(Profiles(:, 2:end), 2)];
    
    mrgFOM.IDFeature(ii, 1) = ii;
    mrgFOM.nbrDetec(ii, 1)  = size(cFOM, 1) - sum(cFOM.M0 == 0);
    mrgFOM.mean_tMax(ii, 1) = mean(nonzeros(cFOM.TimeAtPeakMax));
    mrgFOM.std_tMax(ii, 1) = std(nonzeros(cFOM.TimeAtPeakMax));
    mrgFOM.mean_M1(ii, 1) = mean(nonzeros(cFOM.M1));
    mrgFOM.std_M1(ii, 1) = std(nonzeros(cFOM.M1));
    mrgFOM.mean_AccMass(ii, 1) = mean(nonzeros(cFOM.accuMass));
    mrgFOM.str_AccMass(ii, 1) = std(nonzeros(cFOM.accuMass));
    mrgFOM.mean_Area(ii, 1) = mean(nonzeros(cFOM.M0));
    mrgFOM.std_Area(ii, 1) = std(nonzeros(cFOM.M0));
    mrgFOM.RSD_Area(ii, 1) = std(nonzeros(cFOM.M0))/ mean(nonzeros(cFOM.M0))*100;
    
    VLim = nonzeros(cFOM.Lm1);
    io = isoutlier(VLim);
    mrgFOM.Lm1(ii, 1) = mean(VLim(~io)) - 3*std(VLim(~io));
    
    VLim = nonzeros(cFOM.Lm2);
    io = isoutlier(VLim);
    mrgFOM.Lm2(ii, 1) = mean(VLim(~io)) + 3*std(VLim(~io));
end

obj.QC.Method2.allFOM = allFOM;
obj.QC.Method2.FOM    = struct2table(mrgFOM);
obj.QC.Method2.Profiles = allPrf;

myMaster = obj;
save(fullfile(obj.Path, obj.Name), 'myMaster')


%% 3- Getting the FOM based on results of method 2 -> Method3

clear FOM
fprintf('\n\nSTEP 2: Quantitative analysis\n\n')

Lim = [obj.QC.Method2.FOM.Lm1 obj.QC.Method2.FOM.Lm2];
LimTgt = zeros(size(Lim, 1), 4);

% 3.1 FInd superimposing ROI
for ii = 1:size(obj.QC.Method2.FOM, 1)
    
    %3.3.1 FInd superposing ROI with 1*MW;
    IdX = abs(obj.QC.Method2.FOM.mean_AccMass - obj.QC.Method2.FOM.mean_AccMass(ii))...
        <= obj.QC.Method2.FOM.mean_AccMass(ii)/Mw;
    
    if sum(IdX) > 1
        LimX = [find(IdX), Lim(IdX,:), zeros(sum(IdX), 1)];
        LimX = sortrows(LimX, 2);
        for jj = 2:size(LimX, 1)
            if LimX(jj-1, 3) >= LimX(jj, 2)
                LimX(jj, 2) = min(LimX(jj-1:jj, 2));
                LimX(jj, 3) = max(LimX(jj-1:jj, 3));
                LimX(jj, 4) = 1;
                LimX(jj-1, 2:4) = [NaN NaN NaN];
            end
        end
        
        LimTgt(LimX(:,1),:) = LimX(:, :);
    else
        LimTgt(ii, :) = [ii Lim(ii, :) 0];
    end
end

LimTgt(isnan(LimTgt(:,4)), :) = [];
TimeList    = (LimTgt(:,2) + LimTgt(:,3))/2;
timeWdw_opt = (LimTgt(:,3) - LimTgt(:,2) + timeWdw/5)/2;

for ii = 1:length(obj.QC.Files)
    fprintf('\nChecking %s\n', obj.QC.Files{ii})
    fprintf('\tCreating %s\n', name2)
    MF = load(fullfile(obj.QC.Files{ii}, 'myFinnee.mat'));
    myFinnee = MF.myFinnee;   
    
    if isfile(fullfile(obj.QC.Files{ii}, [name2, '.mat'])) && ~overwrite
        fprintf('\t %s already exist\n', name)
        continue
    else
        myFinnee.Datasets{dts1}.mkMnROI(obj.QC.Method2.FOM.mean_AccMass(LimTgt(:,1)), mzWdw, ...
            TimeList, timeWdw_opt, name2, ...
            obj.QC.Method1.FOM.std_CtrTime);
        
    end
end

clear TR

allFOM               = table();
allFOM.IDFile        = zeros(0);
allFOM.IDFeature     = zeros(0);
allFOM.M0            = zeros(0);
allFOM.M1            = zeros(0);
allFOM.M2            = zeros(0);
allFOM.M3            = zeros(0);
allFOM.TimeAtPeakMax = zeros(0);
allFOM.IntAtPeakMax  = zeros(0);
allFOM.nbrPts        = zeros(0);
allFOM.Lm1           = zeros(0);
allFOM.Lm2           = zeros(0);
allFOM.mean_mz       = zeros(0);
allFOM.std_mz        = zeros(0);
allFOM.accuMass      = zeros(0);

for ii = 1:length(obj.QC.Files)
    
    TR{ii} = load(fullfile(obj.QC.Files{ii}, [name2 '.mat']));
    
    h = waitbar(0,'Processing ROIs');
    for jj = 1:length(TR{ii}.tgtROIs.ROI)
        waitbar(jj/length(TR{ii}.tgtROIs.ROI))
        cROI = TR{ii}.tgtROIs.ROI{jj};
        
        ChromMom =  cROI.getFOM( 3, 11, 5, LimTgt(jj, 2:3));
        
        if isempty(ChromMom)
            warning off
            allFOM.IDFile(end+1)  = ii;
            allFOM.IDFeature(end) = jj;
            warning on
            
        else
            head = table(ii, jj, 'VariableNames', {'IDFile', 'IDFeature'});
            allFOM = [allFOM; [head ChromMom]];
        end
    end
    try close(h); catch, end %#ok<CTCH>
end


mrgFOM.IDFeature    = [];
mrgFOM.nbrDetec     = [];
mrgFOM.mean_tMax    = [];
mrgFOM.std_tMax     = [];
mrgFOM.mean_M1      = [];
mrgFOM.std_M1       = [];
mrgFOM.mean_AccMass = [];
mrgFOM.str_AccMass  = [];
mrgFOM.mean_Area    = [];
mrgFOM.std_Area     = [];
mrgFOM.RSD_Area     = [];
mrgFOM.Lm1          = [];
mrgFOM.Lm2          = [];
allPrf = {};


for ii = 1:length(TR{1}.tgtROIs.ROI)
    cFOM = allFOM(allFOM.IDFeature == ii, :);
    
    Profiles = obj.QC.Axis.AxisX.Data;
    for jj = 1:length(TR)
        if cFOM.M0(jj) == 0
            Profiles(:, end+1) = 0;
            
        else
            cROI = TR{jj}.tgtROIs.ROI{ii};
            
            try
            XYZ = cROI.AxisTm.Data;
            XYZ(:,2) = trapz(cROI.AxisMZ.Data, cROI.StoredData);
            XYZ(isnan(XYZ)) = 0;
            
            is = findCloser(cFOM.Lm1(jj), XYZ(:,1));
            ie = findCloser(cFOM.Lm2(jj), XYZ(:,1));
            
            vq = interp1(XYZ(is:ie, 1), XYZ(is:ie, 2), Profiles(:,1));
            vq(isnan(vq)) = 0;
            Profiles(:, end+1) = vq;
            catch
                disp('wth')
            end
        end
    end
    
    Profiles(1,2:end) = zeros(size(Profiles(1,2:end)));
    Profiles(end,2:end) = zeros(size(Profiles(1,2:end)));
    PrfMinus = circshift(Profiles(:, 2:end), -1);
    PrfPlus  = circshift(Profiles(:, 2:end), 1);
    Profiles(sum([PrfMinus, Profiles(:, 2:end), PrfPlus] == 0, 2) == 3*jj, :) = [];
    allPrf{end+1} = [Profiles(:,1), mean(Profiles(:, 2:end), 2)];
    
    mrgFOM.IDFeature(ii, 1) = ii;
    mrgFOM.nbrDetec(ii, 1)  = size(cFOM, 1) - sum(cFOM.M0 == 0);
    mrgFOM.mean_tMax(ii, 1) = mean(nonzeros(cFOM.TimeAtPeakMax));
    mrgFOM.std_tMax(ii, 1) = std(nonzeros(cFOM.TimeAtPeakMax));
    mrgFOM.mean_M1(ii, 1) = mean(nonzeros(cFOM.M1));
    mrgFOM.std_M1(ii, 1) = std(nonzeros(cFOM.M1));
    mrgFOM.mean_AccMass(ii, 1) = mean(nonzeros(cFOM.accuMass));
    mrgFOM.str_AccMass(ii, 1) = std(nonzeros(cFOM.accuMass));
    mrgFOM.mean_Area(ii, 1) = mean(nonzeros(cFOM.M0));
    mrgFOM.std_Area(ii, 1) = std(nonzeros(cFOM.M0));
    mrgFOM.RSD_Area(ii, 1) = std(nonzeros(cFOM.M0))/ mean(nonzeros(cFOM.M0))*100;
    
    VLim = nonzeros(cFOM.Lm1);
    io = isoutlier(VLim);
    mrgFOM.Lm1(ii, 1) = mean(VLim(~io)) - 3*std(VLim(~io));
    
    VLim = nonzeros(cFOM.Lm2);
    io = isoutlier(VLim);
    mrgFOM.Lm2(ii, 1) = mean(VLim(~io)) + 3*std(VLim(~io));
end

obj.QC.Method3.allFOM = allFOM;
obj.QC.Method3.FOM    = struct2table(mrgFOM);
obj.QC.Method3.Profiles = allPrf;

myMaster = obj;
save(fullfile(obj.Path, obj.Name), 'myMaster');

end
