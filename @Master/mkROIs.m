function obj = mkROIs(obj, dts, name, mzWdw, tmWdw, varargin)

% 1- Initialisation and options
narginchk(2, inf)

%options = checkVarargin(CluLvl, varargin);


for ii = 1:length(obj.QCFiles)
    fprintf('\nChecking %s\n', obj.QCFiles{ii})
    if exist(fullfile(obj.QCFiles{ii}, name), 'file') ~= 2
        fprintf('\tCreating %s\n', name)
        MF = load(fullfile(obj.QCFiles{ii}, 'myFinnee.mat'));
        myFinnee = MF.myFinnee;
        myFinnee.Datasets{dts}.mkMnROI(obj.ClusteredQC.FOM(:,6), mzWdw, obj.ClusteredQC.FOM(:,4), tmWdw, name);
    else
        fprintf('\t%s already exists\n', name)
        
    end
end




%% SUB FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

