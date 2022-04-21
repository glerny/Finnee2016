function obj = QuantSamples(obj, dts, timeWdw, mzWdw, TAG, Filter, QCAnalysis)


name      = 'ROIs4Quant';
Filter   = dcdFilter(Filter);

switch QCAnalysis
    case 1
        disp('TOBEDONE')
        
    case 2
        disp('TOBEDONE')
        
    case 3
        FOM      = obj.QC.Method3.FOM;
        id = FOM.nbrDetec < Filter.freq*max(FOM.nbrDetec)/100 |...
            FOM.RSD_Area > Filter.RSD;
        FOM(id, :) = [];
        FOM = sortrows(FOM, 'IDFeature');
        mzList = FOM.mean_AccMass;
        tmList = (FOM.Lm1 + FOM.Lm2)/2;
        tmWdW  = (FOM.Lm2 - FOM.Lm1 + timeWdw)/2;
        Lim = [FOM.Lm1, FOM.Lm2];
        
    otherwise
        error
end


%% 1- Initialisation
% 1.1- Load Samples
fprintf('\n\n STARTING \n')
if isempty(obj.Samples)
    obj.Samples.Files = {};
    obj.Samples.ID    = [];
    obj.Samples.Tags  = {};
    obj.Samples.allFOM   = {};
    obj.Samples.FOM.IDFeat = table;
    obj.Samples.FOM.Tm   = table;
    obj.Samples.FOM.MZ   = table;
    obj.Samples.FOM.Area = table;
    obj.Samples.Parameters.dts = dts;
    obj.Samples.Parameters.timeWdw = timeWdw;
    obj.Samples.Parameters.mzWdw = mzWdw;
    obj.Samples.Parameters.Filter = Filter;
    obj.Samples.Parameters.FOM = FOM;
    obj.Samples.Parameters.ListTags = {};
    obj.Samples.Parameters.ListFiles = {};
    
end

obj.Samples.Parameters.ListTags{end+1} = TAG;

%1 . Get tgt mz/tm

dirs = uigetdirs(pwd, 'select Finnee Samples folders');
obj.Samples.Parameters.ListFiles = {obj.Samples.Parameters.ListFiles, dirs};

for ii = 1:length(dirs)
    [~, N, ~] = fileparts(dirs{ii});
    fprintf('\t processing %s\n', N)
    
    %1. Check if the files has been already done
    FileName = fullfile(dirs{ii}, 'myFinnee.mat');
    obj.Samples.Files{end+1} = FileName;
    [~, fName, ~] = fileparts(dirs{ii});
    obj.Samples.ID{end+1} = fName;
    obj.Samples.Tags{end+1} = TAG;
    
    ROIname = fullfile(dirs{ii}, name);
    fprintf('\tCreating %s\n', ROIname)
    MF = load( FileName);
    myFinnee = MF.myFinnee;
    TR = myFinnee.Datasets{dts}.mkMnROI(mzList, mzWdw, tmList, ...
        tmWdW, name, []);
    
    Res = table;
    h = waitbar(0,'Processing ROIs');
    for jj = 1:length(TR.ROI)
        waitbar(jj/length(TR.ROI))
        cROI = TR.ROI{jj};
        cFOM = cROI.getFOM(3, 11, 5, Lim(jj, :));
        
        if isempty(cFOM)
            warning off
            if isempty(Res)
                Res.M0(1, 1) = 0;
                Res.M1(1, 1) = 0;
                Res.M2(1, 1) = 0;
                Res.M3(1, 1) = 0;
                Res.TimeAtPeakMax(1, 1) = 0;
                Res.IntAtPeakMax(1, 1) = 0;
                Res.nbrPts(1, 1) = 0;
                Res.Lm1(1, 1) = 0;
                Res.Lm2(1, 1) = 0;
                Res.mean_mz(1, 1) = 0;
                Res.std_mz(1, 1) = 0;
                Res.accuMass(1, 1) = 0;
            else
                Res.M0(end+1, 1) = 0;
            end
            warning on
            
        else
            Res = [Res; cFOM];
        end
    end
    if ishandle(h), close(h); end
    
    Res = [FOM(:, 1), Res];
    files = ['s_', fName];
    files = regexprep(files,'\W','');
    obj.Samples.FOM.IDFeat = addvars(obj.Samples.FOM.IDFeat, Res.IDFeature, 'NewVariableNames',files);
    obj.Samples.FOM.Tm   = addvars(obj.Samples.FOM.Tm, Res.M1, 'NewVariableNames',files);
    obj.Samples.FOM.MZ   = addvars(obj.Samples.FOM.MZ, Res.accuMass, 'NewVariableNames',files);
    obj.Samples.FOM.Area = addvars(obj.Samples.FOM.Area, Res.M0, 'NewVariableNames',files);
    obj.Samples.allFOM{end+1}   = Res;
end

myMaster = obj;
save(fullfile(obj.Path, obj.Name), 'myMaster')


    function FltOut =  dcdFilter(FltIn)
        Stoppers = strfind(FltIn, ';');
        if length(Stoppers) ~= 2
            error('Filter incorrect')
        end
        
        IxF = strfind(FltIn, 'RSD');
        FltOut.RSD  = str2double(FltIn(IxF + 3:Stoppers(1)-1));
        IxF = strfind(FltIn, 'freq');
        FltOut.freq = str2double(FltIn(IxF + 4:Stoppers(2)-1));
    end

end
