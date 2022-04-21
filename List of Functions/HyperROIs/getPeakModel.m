function [HyperROIs, myPeakModel] = getPeakModel(HyperROIs, TgtId)

peakModel = (-6:0.1:6)';
MassDep = [];
if ~isfield(HyperROIs, 'Deconvolution')
    HyperROIs.Deconvolution.Model = [];
    path = fileparts(HyperROIs.myPath);
    HyperROIs.Deconvolution.link2DeconvDat = fullfile(path, 'deconv.dat');
    HyperROIs.Deconvolution.link2DeconvTbl = fullfile(path, 'links2deconv.mat');
end


LD = load(HyperROIs.links2dat);
links2dat = LD.links2dat; clear LD;

if isempty(TgtId)
    TgtId = (1:size(links2dat, 1))';
end

for tgt = 1:size(TgtId, 1)
    ii = find(links2dat.ID == TgtId(tgt));
    file2add = links2dat.Datafiles{ii};
    [fidWriteDat, errmsg]  = fopen(file2add, 'r');
    fseek(fidWriteDat,  0, 'bof');
    data = fread(fidWriteDat,inf, links2dat.format{ii});
    X_Tm = data(1:links2dat.Size(ii,2));
    Y_mz = data(links2dat.Size(ii,2)+1:links2dat.Size(ii,2)+links2dat.Size(ii,1));
    id1 = links2dat.Size(ii,2)+links2dat.Size(ii,1)+1;
    cHROI = reshape(data(id1:end), ...
        links2dat.Size(ii,1), links2dat.Size(ii,2), []);
    proj2MZ = mean(mean(cHROI, 3), 2);
    if size(LocalMaxima([Y_mz, proj2MZ], 3, 0.05*max(proj2MZ)), 1) == 1
        M = ChrMoment([Y_mz, proj2MZ]);
        MassDep(end+1, :) = [ii, M];
        Y_mz = (Y_mz - M(2))/sqrt(M(3));
        proj2MZ = proj2MZ/max(proj2MZ);
        peakModel(:, end+1) = interp1(Y_mz, proj2MZ, peakModel(:,1));
    end
    fclose(fidWriteDat);
end

while 1
    myPeakModel = mean(peakModel(:, 2:end), 2, 'omitnan');
    peakModel(isnan(peakModel)) = 0;
    test = pdist2(myPeakModel', peakModel(:, 2:end)', 'correlation');
    if any(isoutlier(test))
        MassDep(isoutlier(test), :) = [];
        peakModel(:, [false, isoutlier(test)]) = [];
        peakModel(peakModel == 0) = nan;
    else
        break
    end
end

myPeakModel = peakModel(:,1);
myPeakModel(:, 2) = mean(peakModel(:,2:end), 2, 'omitnan');
myPeakModel(:, 3) = std(peakModel(:,2:end), [], 2, 'omitnan');
HyperROIs.Deconvolution.PeakModel = myPeakModel;
HyperROIs.Deconvolution.PeakModel_FOM = MassDep;

save(HyperROIs.myPath, 'HyperROIs')
end

