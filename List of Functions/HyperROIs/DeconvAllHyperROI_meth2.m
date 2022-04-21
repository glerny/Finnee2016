function DeconvAll  = DeconvAllHyperROI_meth2(HyperROIs)

DeconvAll = table();

% fprintf('\n ANALYSING HYPERROIS \n\n')
datetime
parfor ii = 1:size(HyperROIs.Data, 2)
   % fprintf('\t Processing HyperROI %i out of %i (@%s)\n', ii, size(HyperROIs.Data, 2), datetime);
    DeconvMe  = Deconv1HyperROI_meth2(HyperROIs, ii);
    
    if ~isempty(DeconvMe.Id_HyperROIs)
        DeconvAll = [DeconvAll; DeconvMe];
        
    end
    
end
datetime