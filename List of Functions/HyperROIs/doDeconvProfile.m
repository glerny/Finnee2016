function FOM_2 = doDeconvProfile(FOM)

count = 0;
ThreCmpt = 95;
opts = optimset('MaxFunEvals', 1e5, 'MaxIter', 1e5);

count = 0;
for ii = 1:size(FOM.Slice.Profiles, 2)
    cData = [FOM.Slice.axisTm{ii} FOM.Slice.Profiles{ii}];
    alData = doAlignment_MinPearson(cData);
    
    %Model = fittEMGs(alData, FOM.Slice.NbrProf_pca{ii});
    Model = fittPMG2(alData, FOM.Slice.NbrProf_pca{ii});
    Model(:, sum(Model) == 0) = [];
    [~, Id0] = max(Model, [], 1);
    
    clear x0 x
    for jj = 1:size(alData, 2) - 1
        x0 = alData(Id0, jj+1);
        x(:, jj) = fminsearch(@deconv, x0, opts);
    end
    
    for jj = 1:size(x, 1)
        count = count + 1;
        
        FOM_2.Id_HyperROIs(count) = FOM.Slice.Id_HyperROIs(ii);
        FOM_2.Id_Slice(count) = ii;
        FOM_2.Id_Model(count) = jj;
        XY = [FOM.Slice.axisMS{ii} FOM.Slice.ModelMZ{ii}];
        XY(isnan(XY)) = 0;
        M   = ChrMoment(XY);
        FOM_2.M1_mz(count) = M(2);
        FOM_2.M2_mz(count) = M(3);
        FOM_2.M3_mz(count) = M(4);
        FOM_2.Ions{count} = XY;
        
        XY = [FOM.Slice.axisTm{ii} Model(:, jj)];
        XY(isnan(XY)) = 0;
        M   = ChrMoment(XY);
        FOM_2.M1_tm(count) = M(2);
        FOM_2.M2_tm(count) = M(3);
        FOM_2.M3_tm(count) = M(4);
        FOM_2.Profile{count} = XY;
        
        FOM_2.Areas(:, count) = x(jj, :);
        
    end
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
        
        f = sqrt(sum((alData(:,jj+1) - sum(Y, 2, 'omitnan')).^2, 'omitnan'));
        
        if isnan(f)
            disp('pp')
        end
    end
end