function [FOM_MS, FOM_Prof, HRShift]  = DeconvAllHyperROI(HyperROIs)

FOM_MS = table();
FOM_Prof = table();

fprintf('\n ANALYSING HYPERROIS \n\n')
datetime
parfor ii = 1:size(HyperROIs.Data, 2)
    fprintf('\t Processing HyperROI %i out of %i (@%s)\n', ii, size(HyperROIs.Data, 2), datetime);
    try
        [c_FOM_MS, c_FOM_Prof, HRShift(ii, :)]  = Deconv1HyperROI(HyperROIs, ii);
        try
            if ~isempty(c_FOM_Prof.Id)
                FOM_MS = [FOM_MS; c_FOM_MS];
                FOM_Prof = [FOM_Prof; c_FOM_Prof];
            end
        catch  ME
            disp(ii)
            rethrow(ME)
        end
    catch ME
        disp(ii)
        rethrow(ME)
    end
end

save FOMsOrbiE_end.mat FOM_MS FOM_Prof HRShift
datetime