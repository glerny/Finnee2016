function [FOM, HyperROIs] = doDeconv(HyperROIs)

%% Find model spectra and prediction for MS spectra

ModelSpectra(:, 1) = -6:0.1:6;
XY = [];
for ii = 1:size(HyperROIs, 1)
    cSpectra = [HyperROIs.axeMz{ii} sum(sum(HyperROIs.Data{ii}, 2, 'omitnan'), 3, 'omitnan')];
    
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
for ii = 1:size(HyperMat, 1)
    ii
    tgtSpectra = [HyperMat.axeMz{ii} ...
        sum(sum(HyperMat.Data{ii}, 2, 'omitnan'), 3, 'omitnan')];
    lm = LocalMaxima(tgtSpectra, 1, max(cSpectra(:,2))*0.01);
    HyperMat.LocalMaxima{ii} = lm;
    
    pca_Data = squeeze(sum(HyperMat.Data{ii}, 2, 'omitnan'));
    pca_Data(isnan(pca_Data)) = 0;
    warning off
    [~,~,~, ~, explained] = pca(pca_Data', 'algorithm','als');
    warning on
    NbrSpct = find(cumsum(explained) >= ThreCmpt, 1, 'first');
    if isempty(NbrSpct), NbrSpct = 1; end
    HyperMat.NbrMSSpect_pca{ii} = NbrSpct;
    
    
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
    
    % Model(isnan(Model)) = 0;
    
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
    HyperMat.ModelMZ_SSR(ii) = F;
    HyperMat.ModelMZ{ii} = Model./max(Model);
    
    data2deconv = HyperMat.Data{ii};
    Model = HyperMat.ModelMZ{ii};
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
    HyperMat.deconMS{ii} = deconMS;
    
    clear NbrPrf
    for jj = 1:size(deconMS, 2)
        nnmf_Data = squeeze(deconMS(:, jj, :));
        
        % find nbr of components
        [coeff,score,~, ~, explained] = pca(nnmf_Data');
        xNbr = find(cumsum(explained) >= ThreCmpt, 1, 'first');
        
        if isempty(xNbr)
            NbrPrf(jj) = 1;
        else
            NbrPrf(jj) = xNbr;
        end
        
        WO = coeff(:, 1: NbrPrf(jj));
        WO(WO < 0) = 0;
        H0 = score(1: NbrPrf(jj), :);
        H0(H0 <0) = 0;
        
        [W,H] = nnmf(nnmf_Data', NbrPrf(jj),...
            'H0', WO');
               %  'Options',opt,...
                % 'Algorithm','als');
             
    end
    HyperMat.NbrProf_pca{ii} = NbrPrf;
    
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




