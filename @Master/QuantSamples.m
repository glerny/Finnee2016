function obj = QuantSamples(obj, name, dts , mzWdw, tmWdw, eTm, eMz, CV)

if isempty(obj.SamplesFiles)
    obj.SamplesFiles = {};
end

if ~isfield(obj.QuantitAnalysis, 'Samples')
    obj.QuantitAnalysis.Samples.Tm   = table;
    obj.QuantitAnalysis.Samples.MZ   = table;
    obj.QuantitAnalysis.Samples.Area = table;
end
%1 . Get tgt mz/tm
FOM = obj.QuantitAnalysis.QC.SummaryStep2;
id = FOM.NonNull == max(FOM.NonNull) & FOM.RSD <= CV;
tgTm = FOM.tm(id);
tgMz = FOM.mz(id);
Lim = obj.QuantitAnalysis.QC.Lim(id, :);

dirs = uigetdirs(pwd, 'select Finnee Samples folders');
for ii = 1:length(dirs)
    disp(datetime), disp(dirs{ii});
    [~, N, ~] = fileparts(dirs{ii});
    disp(N);
    %1. Check if the files has been already done
    FileName = fullfile(dirs{ii}, 'myFinnee.mat');
    if isempty(obj.SamplesFiles)
        obj.SamplesFiles{end+1} = FileName;
        [~, fName, ~] = fileparts(dirs{1});
    else
        if contains(obj.SamplesFiles, FileName)
            error('file exist')
        else
            obj.SamplesFiles{end+1} = FileName;
            [~, fName, ~] = fileparts(dirs{ii});
        end
    end
    
    %2. CHeck if name (ROIs) exist otherwise create it
    
    ROIname = fullfile(dirs{ii}, name);
    fprintf('\nChecking %s\n', dirs{ii})
    if exist(ROIname, 'file') ~= 2
        fprintf('\tCreating %s\n', ROIname)
        MF = load( FileName);
        myFinnee = MF.myFinnee;
        TR = myFinnee.Datasets{dts}.mkMnROI(tgMz, mzWdw, tgTm, tmWdw, name);
    else
        fprintf('\t%s already exists\n', name);
        TR = load(ROIname);
        TR = TR.tgtROIs;
    end
    
    Res = [];
    h = waitbar(0,'Processing ROIs');
    for jj = 1:length(TR.ROI)
        waitbar(jj/length(TR.ROI))
        cROI = TR.ROI{jj};
        Res( jj, :) = cROI.getFOM(cROI.TgtTm, eTm, cROI.TgtMz, eMz, Lim(jj, :));
    end
    if ishandle(h), close(h); end
    
    obj.QuantitAnalysis.Samples.Tm   = addvars(obj.QuantitAnalysis.Samples.Tm, Res(:,2), 'NewVariableNames',['s', fName]);
    obj.QuantitAnalysis.Samples.MZ   = addvars(obj.QuantitAnalysis.Samples.MZ, Res(:,7), 'NewVariableNames',['s', fName]);
    obj.QuantitAnalysis.Samples.Area = addvars(obj.QuantitAnalysis.Samples.Area, Res(:,1), 'NewVariableNames',['s', fName]);
    myMaster = obj;
    save(fullfile(obj.Path, obj.Name), 'myMaster')
end
end
