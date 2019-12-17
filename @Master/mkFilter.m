function obj = mkFilter(obj, N, RSD)

% 1- Initialisation and options


obj.ClusteredQC.Filter = obj.QuantitAnalysis.QC.SummaryStep2.NonNull >=...
    max(obj.QuantitAnalysis.QC.SummaryStep2.NonNull)*N/100 &...
    obj.QuantitAnalysis.QC.SummaryStep2.RSD < RSD;

myMaster = obj;
try
    save(fullfile(obj.Path, obj.Name), 'myMaster')
catch
    uisave('myMaster')
end







%% SUB FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

