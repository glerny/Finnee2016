function obj = mkROIs(obj, dts, name, mzWdw, tmWdw, overwrite, varargin)

% 1- Initialisation and options
narginchk(2, inf)
maxRoi = 20000;

%options = checkVarargin(CluLvl, varargin);


for ii = 1:length(obj.QC.Files)
    fprintf('\nChecking %s\n', obj.QC.Files{ii})
    
    if isfile(fullfile(obj.QC.Files{ii}, [name, '.mat'])) && ~overwrite
        
        fprintf('\t %s already exist\n', name)
        continue
    else
        
        fprintf('\tDeleting %s \n', name)
    end

    fprintf('\tCreating %s\n', name)
    MF = load(fullfile(obj.QC.Files{ii}, 'myFinnee.mat'));
    myFinnee = MF.myFinnee;
    
    if length(obj.QC.Method1.FOM.AccurateMass) <= maxRoi
        myFinnee.Datasets{dts}.mkMnROI(obj.QC.Method1.FOM.AccurateMass, mzWdw, ...
            obj.QC.Method1.FOM.CtrTime, tmWdw, name, ...
            obj.QC.Method1.FOM.std_CtrTime);
        
    else
        is = 1;
        ie = maxRoi;
           tgtROIs =  myFinnee.Datasets{dts}.mkMnROI(...
               obj.QC.Method1.FOM.AccurateMass(is:ie), mzWdw, ...
               obj.QC.Method1.FOM.CtrTime(is:ie), tmWdw, '', ...
               obj.QC.Method1.FOM.std_CtrTime(is:ie));
           
        for jj = maxRoi+1:maxRoi:length(obj.QC.Method1.FOM.AccurateMass)
           is = jj;
           ie = min(jj + maxRoi -1, length(obj.QC.Method1.FOM.AccurateMass));
           cRoi =  myFinnee.Datasets{dts}.mkMnROI(...
               obj.QC.Method1.FOM.AccurateMass(is:ie), mzWdw, ...
               obj.QC.Method1.FOM.CtrTime(is:ie), tmWdw, '', ...
               obj.QC.Method1.FOM.std_CtrTime(is:ie));
           tgtROIs.mzList = [tgtROIs.mzList; cRoi.mzList];
           tgtROIs.tmList = [tgtROIs.tmList; cRoi.tmList];
           tgtROIs.ROI    = [tgtROIs.ROI, cRoi.ROI];
        end
    save(fullfile(myFinnee.Path2Fin, name), 'tgtROIs')
    end
    
end




%% SUB FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

