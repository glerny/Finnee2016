function obj = AnalyseSamples(obj, name, dts)


%% 1- Initialisation
% 1.1- Load Samples
fprintf('\n\n STARTING \n')
if isempty(obj.Samples)
    obj.Samples.Files = {};
    obj.Samples.Id    = [];
    obj.Samples.FOM.Tm   = table;
    obj.Samples.FOM.MZ   = table;
    obj.Samples.FOM.Area = table;
    obj.Samples.Parameters.dts   = dts;
end

%1 . Get tgt mz/tm
FOM1     = obj.QC.Method1.FOM; FOM1 = sortrows(FOM1, 'ID');
FOM2     = obj.QC.Method2.FOM; FOM2 = sortrows(FOM2, 'ID');
FOM3     = obj.QC.Method3.FOM; FOM3 = sortrows(FOM3, 'ID');


Lim      = [FOM3.TmStart FOM3.TmEnd];
tgTm     = FOM1.CtrTime;
tgMz     = FOM1.AccurateMass;
tgVar    = FOM1.std_CtrTime;
obj.Samples.Parameters.eTM = inf;
obj.Samples.Parameters.eMz = inf;
TimeList    = (FOM3.TmEnd + FOM3.TmStart)/2;
timeWdw_opt = FOM3.TmEnd - FOM3.TmStart;
timeWdw_opt(isnan(TimeList)) = NaN;
TimeList(isnan(TimeList)) = (max(TimeList) + min(TimeList))/2;
Ic = isnan(timeWdw_opt);
timeWdw_opt = timeWdw_opt + 0.5;
timeWdw_opt(Ic) = 0;

dirs = uigetdirs(pwd, 'select Finnee Samples folders');
for ii = 1:length(dirs)
    disp(datetime), disp(dirs{ii});
    [~, N, ~] = fileparts(dirs{ii});
    disp(N);
    %1. Check if the files has been already done
    FileName = fullfile(dirs{ii}, 'myFinnee.mat');
    if isempty(obj.Samples.Files)
        obj.Samples.Files{end+1} = FileName;
        [~, fName, ~] = fileparts(dirs{1});
    else
        if contains(obj.Samples.Files, FileName)
            error('file exist')
        else
            obj.Samples.Files{end+1} = FileName;
            [~, fName, ~] = fileparts(dirs{ii});
        end
    end
    
    %2. CHeck if name (ROIs) exist otherwise create it
    
    ROIname = fullfile(dirs{ii}, name);
    fprintf('\nChecking %s\n', dirs{ii})
%     if exist(ROIname, 'file') ~= 2
        fprintf('\tCreating %s\n', ROIname)
        MF = load( FileName);
        myFinnee = MF.myFinnee;
        TR = myFinnee.Datasets{dts}.mkMnROI(tgMz, 5, TimeList, ...
            timeWdw_opt, name, tgVar);
%     else
%         fprintf('\t%s already exists\n', name);
%         TR = load(ROIname);
%         TR = TR.tgtROIs;
%     end
    
    eTm = obj.Samples.Parameters.eTM;
    eMz = obj.Samples.Parameters.eMz;
    Res = [];
    h = waitbar(0,'Processing ROIs');
    for jj = 1:length(TR.ROI)
        waitbar(jj/length(TR.ROI))
        cROI = TR.ROI{jj};
        Res( jj, :) = cROI.getFOM(tgTm(jj), eTm, cROI.TgtMz, eMz, Lim(jj, :));
    end
    if ishandle(h), close(h); end
    
    files = ['s', fName];
    files = regexprep(files,'\W','');
    obj.Samples.FOM.Tm   = addvars(obj.Samples.FOM.Tm, Res(:,2), 'NewVariableNames',files);
    obj.Samples.FOM.MZ   = addvars(obj.Samples.FOM.MZ, Res(:,11), 'NewVariableNames',files);
    obj.Samples.FOM.Area = addvars(obj.Samples.FOM.Area, Res(:,1), 'NewVariableNames',files);
    myMaster = obj;
    save(fullfile(obj.Path, obj.Name), 'myMaster')
end
end
