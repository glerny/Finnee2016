function [Model, Profiles, Data] = doDeconvMS(Data2Deconv, ModelSpectra, p)

%% Find model spectra and prediction for MS spectra
ThreCmpt = 95;
aveSpectra = [ModelSpectra(:, 1) mean( ModelSpectra(:,2:end), 2, 'omitnan')];

%%%%%

for ii = 1:size(HyperROIs.HyperROIs.Data, 1)
    ii
    FOM.DeconvMS.Id_HyperROIs(ii) = ii;
    tgtSpectra = HyperROIs.HyperROIs.axeMz{ii};
    lm = LocalMaxima(tgtSpectra(:, 1:2), 2, max(tgtSpectra(:,2))*0.01);
    FOM.DeconvMS.lm{ii} = lm;
    
    pca_Data = squeeze(sum(HyperROIs.HyperROIs.Data{ii}, 2, 'omitnan'));
    pca_Data(isnan(pca_Data)) = 0;
    warning off
    [~,~,~, ~, explained] = pca(pca_Data', 'algorithm','als');
    warning on
    NbrSpct = find(cumsum(explained) >= ThreCmpt, 1, 'first');
    if isempty(NbrSpct), NbrSpct = 1; end
    FOM.DeconvMS.NbrMSSpect_pca(ii) = NbrSpct;
    
    
    Model = [];
    step = (tgtSpectra(2,1) - tgtSpectra(1,1))/5;
    axeSpectra = tgtSpectra(1,1) - 5*step:step:tgtSpectra(end,1) + 5*step;
    
    clear x0
    for jj = 1:size(lm, 1)
        cModel = aveSpectra;
        cModel(:,1) = cModel(:,1)*sqrt(M2) + lm(jj,1);
        cModel(:,2) = cModel(:,2)* lm(jj,2);
        Model(:, end+1) = interp1(cModel(:,1), cModel(:,2), axeSpectra);
        x0(2*jj-1) = 0.00005;
        x0(2*jj) = 1;
    end
    
    opts = optimset('Display','off', 'MaxFunEvals' , 5000);
    x = fminsearch(@myFun, x0, opts);
    xold      = x;
    Model_old = Model;
    
    while size(Model, 2) < NbrSpct
        [Fs, myFs] = myFun(x);
        [IMax, IdAdd] = max(tgtSpectra(:,2) - sum(myFs, 2));
        cModel = aveSpectra;
        cModel(:,1) = cModel(:,1)*sqrt(M2) + tgtSpectra(IdAdd, 1);
        cModel(:,2) = cModel(:,2)* IMax;
        Model(:, end+1) = interp1(cModel(:,1), cModel(:,2), axeSpectra);
        x(end+1) = 0.00005;
        x(end+1) = 1;
        x = fminsearch(@myFun, x, opts);
        [Fe, myFe] = myFun(x);
        
        if Fs <= Fe
            break
            
        else
            xold      = x;
            Model_old = Model;
            
        end
    end
    x = xold;
    Model = Model_old;
    [F, Model] = myFun(x);
    FOM.DeconvMS.ModelMZ_SSR(ii) = F;
    FOM.DeconvMS.ModelMZ{ii} = Model./max(Model);
    
    data2deconv = HyperROIs.HyperROIs.Data{ii};
    Model = FOM.DeconvMS.ModelMZ{ii};
    Model(isnan(Model)) = 0;
    
    deconMS = [];
    for jj = 1:size(data2deconv, 3)
        cData = squeeze(data2deconv(:, :, jj));
        cData(isnan(cData)) = 0;
        
        % deconMS(:, :, jj) = cData'/Model';
        for kk = 1:size(cData, 2)
            x0 = [];
            for mdl = 1:size(Model, 2)
                x0(mdl) = cData(Model(:, mdl) == max(Model(:, mdl)), kk);
            end
            x0(x0 < 0) = 0;
            [x, ~, exitflag] = fminsearch(@deconv, x0, opts);
            if exitflag == 0
                count = 1;
                while 1
                    [x, ~, exitflag] = fminsearch(@deconv, x, opts);
                    if exitflag == 1, break, end
                    if count > 10, break, end
                    count = count+1;
                end
                
            end
            deconMS(kk, :, jj) = x;
        end
    end
    FOM.DeconvMS.axisMS{ii} = HyperROIs.HyperROIs.axeMz{ii};
    FOM.DeconvMS.axisTm{ii} = HyperROIs.HyperROIs.axeTm{ii};
    FOM.DeconvMS.deconMS{ii} = deconMS;
end

    function [f, cModel] = myFun(x)
        
        cModel = nan(size(tgtSpectra, 1), size(Model, 2));
        for k = 1:size(Model, 2)
            if x(2*k) >= 0
                cModel(:, k) = interp1(axeSpectra + x(2*k-1), ...
                    Model(:, k)*x(2*k), tgtSpectra(:,1));
                
            else
                cModel(:, k) = inf;
                
            end
        end
        % cModel(isnan(cModel)) = 0;
        
        f = sqrt(sum((tgtSpectra(:,2) - sum(cModel, 2, 'omitnan')).^2));
    end

    function [f, Y] = deconv(x)
        
        Y = zeros(size(Model));
        for k = 1:size(Model, 2)
            if x(k) > 0
                Y(:, k) = x(k)*Model(:, k);
            else
                Y(:, k) = inf;
            end
        end
        
        f = sqrt(sum((cData(:,kk) - sum(Y, 2)).^2));
    end

end




