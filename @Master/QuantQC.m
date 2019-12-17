function obj = QuantQC(obj, name, eTm, eMz, varargin)

% 1- Initialisation and options

fprintf('\n\nSTEP 1: Getting limits\n\n')
for ii = 1:length(obj.QCFiles)

    fprintf('\t processing %s\n', obj.QCFiles{ii})
    TR = load(fullfile(obj.QCFiles{ii}, name));

    h = waitbar(0,'Processing ROIs');
    for jj = 1:length(TR.tgtROIs.ROI)
        waitbar(jj/length(TR.tgtROIs.ROI))
        cROI = TR.tgtROIs.ROI{jj};
        FOM(ii, jj, :) = cROI.getFOM(cROI.TgtTm, eTm, cROI.TgtMz, eMz);
    end
    try close(h); catch, end %#ok<CTCH>
end
obj.QuantitAnalysis.QC.Step1 = FOM;

nthre = obj.ClusteredQC.paramters.nthre;
FOM = obj.QuantitAnalysis.QC.Step1;
Lim = zeros(size(FOM, 2), 2);
for ii = 1:size(FOM, 2)
    cFOM = squeeze(FOM(:,ii,:));

    if sum(isnan(cFOM(:,1))) < (1-nthre/100)*size(cFOM, 1)

        TF = isoutlier(cFOM(:,9));
        Lim(ii, 1) = mean(cFOM(~TF,9), 'omitnan') -...
            2* std(cFOM(~TF,9), [], 'omitnan');
        TF = isoutlier(cFOM(:,10));
        Lim(ii, 2) = mean(cFOM(~TF,10), 'omitnan') +...
            2* std(cFOM(~TF,10), [], 'omitnan');
    else
        Lim(ii, :) = [NaN, NaN];
    end
end


clear FOM
fprintf('\n\nSTEP 2: Quantitative analysis\n\n')
for ii = 1:length(obj.QCFiles)

    fprintf('\t processing %s\n', obj.QCFiles{ii})
    TR = load(fullfile(obj.QCFiles{ii}, name));

    h = waitbar(0,'Processing ROIs');
    for jj = 1:length(TR.tgtROIs.ROI)
        waitbar(jj/length(TR.tgtROIs.ROI))
        cROI = TR.tgtROIs.ROI{jj};
        FOM(ii, jj, :) = cROI.getFOM(cROI.TgtTm, eTm, cROI.TgtMz, eMz, Lim(jj, :));
    end
    try close(h); catch, end %#ok<CTCH>
end

obj.QuantitAnalysis.QC.Lim = Lim;
obj.QuantitAnalysis.QC.Step2 = FOM;
myMaster = obj;
save(fullfile(obj.Path, obj.Name), 'myMaster')

FOM = obj.QuantitAnalysis.QC.Step2;
sz = [size(FOM, 2), 8];
varTypes = {'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'};
varNames = {'Id', 'NonNull', 'tm', 'std_tm', 'mz', 'std_mz', 'Area', 'RSD'};
SumStep1 = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
SumStep2 = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);

for ii = 1:size(FOM, 2)
    cFOM = squeeze(obj.QuantitAnalysis.QC.Step1(:,ii,:));
    if sum(isnan(cFOM(:,1))) < (1-nthre/100)*size(cFOM, 1)
        SumStep1(ii, :) = {ii, size(cFOM(:,1), 1) - sum(isnan(cFOM(:,1))), ...
            mean(cFOM(:,2), 'omitnan'), std(cFOM(:,2), [], 'omitnan'), ...
            mean(cFOM(:,11), 'omitnan'), std(cFOM(:,11), [], 'omitnan'), ...
            mean(cFOM(:,1), 'omitnan'), ...
            std(cFOM(:,1), [], 'omitnan')./mean(cFOM(:,1), 'omitnan')*100};
    else
        SumStep1(ii, :) = {NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN}; 
    end
    
    cFOM = squeeze(obj.QuantitAnalysis.QC.Step2(:,ii,:));
    if sum(isnan(cFOM(:,1))) < (1-nthre/100)*size(cFOM, 1)
        SumStep2(ii, :) = {ii, size(cFOM(:,1), 1) - sum(isnan(cFOM(:,1))), ...
            mean(cFOM(:,2), 'omitnan'), std(cFOM(:,2), [], 'omitnan'), ...
            mean(cFOM(:,11), 'omitnan'), std(cFOM(:,11), [], 'omitnan'), ...
            mean(cFOM(:,1), 'omitnan'), ...
            std(cFOM(:,1), [], 'omitnan')./mean(cFOM(:,1), 'omitnan')*100};
    else
        SumStep2(ii, :) = {NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN}; 
    end
end

obj.QuantitAnalysis.QC.SummaryStep1 = SumStep1;
obj.QuantitAnalysis.QC.SummaryStep2 = SumStep2;
myMaster = obj;
save(fullfile(obj.Path, obj.Name), 'myMaster')

end
