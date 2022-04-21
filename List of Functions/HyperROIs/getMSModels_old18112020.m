function [Model, f] = getMSModels(ModelSpectra, tgtSpectra, p)

%% Find model spectra and prediction for MS spectra
ThreCmpt = 1;
ThresImprov = 0.95;
%%%%%

% 1. find Optimal number of spectra with PCA and local max
lm = LocalMaxima(tgtSpectra, 1, max(tgtSpectra(:,2))*(1-ThreCmpt)/100);

opts = optimset('Display','off', 'MaxFunEvals' , 5000);
[lm, f] = fminsearch(@myFun, lm, opts);

while 1
    [lm_c, f_c] = fminsearch(@myFun, lm, opts);
    if f_c <= 0.9995*f
        lm = lm_c;
        f = f_c;
    else
        break
    end
end

[f, Model] = myFun(lm);

while 1
    Residuals = tgtSpectra(:,2) - sum(Model, 2, 'omitnan');
    [Val2Res, Id2Res] = max(Residuals);
    lm(end+1, :) = [tgtSpectra(Id2Res,1), Val2Res];
    [lm_res, f_res] = fminsearch(@myFun, lm, opts);
    
    while 1
        [lm_c, f_c] = fminsearch(@myFun, lm_res, opts);
        if f_c <= 0.9995*f_res
            lm_res = lm_c;
            f_res = f_c;
        else
            break
        end
    end
    
    if f_res <= ThresImprov*f
        lm = lm_res;
        [f, Model] = myFun(lm);
    else
        lm(end, :) = [];
        lm = fminsearch(@myFun, lm, opts);
        [f, Model] = myFun(lm);
        break
    end
end
% Remove model with less than 3 non zeros pts
Model(isnan(Model)) = 0;
lm(sum(Model ~= 0) <= 3, :) = [];
lm = fminsearch(@myFun, lm, opts);
[f, Model] = myFun(lm);

Model = [tgtSpectra(:, 1), Model];
Model(isnan(Model)) = 0;

    function [f, Model] = myFun(lm)
        
        Model = [];
        for jj = 1:size(lm, 1)
            cModel = ModelSpectra;
            M2 = polyval(p, lm(jj,1));
            cModel(:,1) = cModel(:,1)*sqrt(M2) + lm(jj,1);
            if lm(jj, 1) >= 0 && lm(jj, 2) >= 0
                cModel(:,2) = cModel(:,2)* lm(jj,2);
                Model(:, end+1) = interp1(cModel(:,1), ...
                    cModel(:,2), tgtSpectra(:, 1));
            else
                Model(:, end+1) = inf(size(tgtSpectra(:, 1)));
            end
        end
        
        f = sqrt(sum((tgtSpectra(:,2) - sum(Model, 2, 'omitnan')).^2));
    end

end




