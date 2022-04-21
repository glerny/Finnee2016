function [FOM_HROIs, FOM_MS, FOM_Prof]  = DeconvAllHyperROI_cleaned(HyperROIs)

FOM_HROIs = table();
FOM_MS    = table();
FOM_Prof  = table();

fprintf('\n ANALYSING HYPERROIS \n\n')
datetime
parfor ii = 1:size(HyperROIs.Data, 2)
    try
        %fprintf('\t Processing HyperROI %i out of %i (@%s)\n', ii, size(HyperROIs.Data, 2), datetime);
        [c_FOM_HROIs, c_FOM_MS, c_FOM_Prof]  = Deconv1HyperROI_cleaned(HyperROIs, ii, 'small');
        FOM_HROIs = [FOM_HROIs; c_FOM_HROIs];
        
        if ~isempty(c_FOM_MS)
            FOM_MS = [FOM_MS; c_FOM_MS];
        end
        
        if ~isempty(c_FOM_Prof)
            FOM_Prof = [FOM_Prof; c_FOM_Prof];
        end
    catch ME
        fprintf('\t Processing HyperROI %i out of %i (@%s)\n', ii, size(HyperROIs.Data, 2), datetime);
        ii
        rethrow(ME)
    end
end

save FOMsOrbiE_end.mat FOM_MS FOM_Prof FOM_HROIs
datetime