function obj = doSpectralDeconv(obj, option)

name     = 'ROIs4Quant.mat';
perc     = 25;
ord_time = 1;
ord_mz   = 1;
p_tm     = {};
p_mz     = {};

%% 0- Alignment of all features
FOM_QC  = obj.QC.Method3.FOM;
FOM_STm = table2array(obj.Samples.FOM.Tm);
FOM_SMZ = table2array(obj.Samples.FOM.MZ);
FOM_SAr = table2array(obj.Samples.FOM.Area);

for ii = 1:length(obj.Samples.Files)
    IdTft = FOM_SAr(:, 1) >= prctile(FOM_SAr(:, 1), perc);
    
    % 0.1 time dimension
    XY = [FOM_STm(IdTft,1), FOM_QC.mean_M1(IdTft) - FOM_STm(IdTft,1)];
    cp = polyfit(XY(:,1), XY(:,2), ord_time);
    io = isoutlier(XY(:,2) - polyval(cp, XY(:,1)));
    XY(io, :) = [];
    p_tm{ii} = polyfit(XY(:,1), XY(:,2), ord_time);
    
    % 0.1 mz dimension
    XY = [FOM_SMZ(IdTft,1), ...
        (FOM_QC.mean_AccMass(IdTft) - FOM_SMZ(IdTft,1))./FOM_STm(IdTft,1)];
    cp = polyfit(XY(:,1), XY(:,2), ord_mz);
    io = isoutlier(XY(:,2) - polyval(cp, XY(:,1)));
    XY(io, :) = [];
    p_mz{ii} = polyfit(XY(:,1), XY(:,2), ord_mz);
    
end


%% 1- Load ROIs and correct tm & MZ
allROIs = {};
ModelSpectra(:,1) = -6:0.1:6;
XY = [];
HyperMat.Data  = {};
HyperMat.axeMz = {};
HyperMat.axeTm = {};

for ii = 1:length(obj.Samples.Files)
    path = fileparts(obj.Samples.Files{ii});
    ROI = load(fullfile(path, name));
    
    for jj = 1:length(ROI.tgtROIs.ROI)
        cROI = ROI.tgtROIs.ROI{jj};
        
        if ii == 1
            HyperMat.axeMz{jj, 1} = cROI.AxisMZ.Data;
            HyperMat.axeTm{jj, 1} = cROI.AxisTm.Data;
            HyperMat.Data{jj, 1}  = NaN(length(HyperMat.axeMz{jj}), length(HyperMat.axeTm{jj}), length(obj.Samples.Files));
        end
        [Xq, Yq] = meshgrid(HyperMat.axeTm{jj}, HyperMat.axeMz{jj});
        [X, Y] = meshgrid(cROI.AxisTm.Data + polyval(p_tm{ii}, cROI.AxisTm.Data), ...
            cROI.AxisMZ.Data + polyval(p_mz{ii}, cROI.AxisMZ.Data).* cROI.AxisMZ.Data);
        cData = cROI.StoredData;
        HyperMat.Data{jj}(:, :, ii) = interp2(X, Y, cData, Xq, Yq);
        
    end
end

%% Find model spectra and prediction for MS spectra
HyperMat = struct2table(HyperMat);

ModelSpectra(:, 1) = -6:0.1:6;
XY = [];
for ii = 1:size(HyperMat, 1)
    cSpectra = [HyperMat.axeMz{ii} sum(sum(HyperMat.Data{ii}, 2, 'omitnan'), 3, 'omitnan')];
    
    if size(LocalMaxima(cSpectra, 1, max(cSpectra(:,2))*0.01), 1) == 1
        M = ChrMoment(cSpectra);
        XY(end+1, :) = M;
        cSpectra(:,1) = (cSpectra(:,1) - M(2))/sqrt(M(3));
        cSpectra(:,2) = cSpectra(:,2)/max(cSpectra(:,2));
        vq = interp1(cSpectra(:,1), cSpectra(:,2), ModelSpectra(:,1));
        ModelSpectra(:, end+1) = vq;
    end
end

%%2- Find and remove outliers
while 1
    aveSpectra = mean( ModelSpectra(:,2:end), 2, 'omitnan');
    R = corrcoef([aveSpectra,  ModelSpectra(:,2:end)], 'rows', 'complete');
    io = isoutlier(R(2:end, 1));
    
    if any(io)
        XY(io, :) = [];
        ModelSpectra(:, [false; io]) = [];
    else
        break
    end
end

aveSpectra = [ModelSpectra(:, 1) mean( ModelSpectra(:,2:end), 2, 'omitnan')];


%% 3- Spectral deconvolution
%%!!!! DO XY PROJ FOR NOW
M2 = mean(XY(:,3));
for ii = 2:size(HyperMat, 1)
    tgtSpectra = [HyperMat.axeMz{ii} ...
        sum(sum(HyperMat.Data{ii}, 2, 'omitnan'), 3, 'omitnan')];
    lm = LocalMaxima(tgtSpectra, 1, max(cSpectra(:,2))*0.01);
    
    Model = [];
    clear x0;
    for jj = 1:size(lm, 1)
        cModel = aveSpectra;
        cModel(:,1) = cModel(:,1)*sqrt(M2) + lm(jj,1);
        cModel(:,2) = cModel(:,2)* lm(jj,2);
        Model(:, end+1) = interp1(cModel(:,1), cModel(:,2), tgtSpectra(:,1));
        x0(2*jj-1) = 0.00005;
        x0(2*jj) = 1;
        
    end
    Model(:, end+1) = ones(size(Model, 1), 1);
    Model(isnan(Model)) = 0;
    x0(end+1) = 1;
    
    opts = optimset('Display','off');
    [f0, myF0] = myFun(x0);
    x = fminsearch(@myFun, x0);
    [HyperMat.ModelMZ_SSR(ii), HyperMat.ModelMZ{ii}] = myFun(x);
    
%     data2deconv = HyperMat.Data{ii};
%     deconMS = [];
%     for jj = 1:size(data2deconv, 3)
%         cData = squeeze(data2deconv(:, :, jj));
%         cData(isnan(cData)) = 0;
%         deconMS(:, :, jj) = cData'/HyperMat.ModelMZ{ii}';
%     end
    
end


    function [f, cModel] = myFun(x)
        
        cModel = zeros(size(Model));
        for k = 1:size(Model, 2)-1
            cModel(:, k) = interp1(tgtSpectra(:,1) + x(2*k-1), ...
                Model(:, k)*x(2*k), tgtSpectra(:,1));
        end
        cModel(:, end) = x(end);
        cModel(isnan(cModel)) = 0;
        
        f = sqrt(sum((tgtSpectra(:,2) - sum(cModel, 2)).^2));
    end

disp('dd')
end




