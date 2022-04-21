function [Model, FittData, SSR] = getMSModels_constrained_cleaned(ModelSpectra, tgtSpectra, Data, Thres1, Thres2, PtsPerPeaks, Options)

%% Find model spectra and prediction for MS spectra
ThresOpti                  = 0.95;
minSpectra                 = 3;
PearsonThresholdClustering = 0.5;
WeightNeg                  = 1;
%%%%%

warning off
%% Find spectra to fit
Dms = mean(tgtSpectra(2:end, 1) - tgtSpectra(1:end-1, 1));
minVar   = (Dms*PtsPerPeaks/30)^2; %Peakwith at base (6 sigma) => PtsPerPeaks/5
idealVar = (Dms*PtsPerPeaks/6)^2; %Peakwith at base (6 sigma) => PtsPerPeaks
maxVar   = (Dms*5*PtsPerPeaks/6)^2; %Peakwith at base (6 sigma) => 5*PtsPerPeaks


sqData = reshape(Data, size(Data, 1), []);
sqData(:, sum(sqData ~= 0, 1) <= 1) = [];

Z = linkage(sqData', 'complete', 'correlation');

Id2Cluster = cluster(Z, 'cutoff', PearsonThresholdClustering, 'Criterion','distance');
ClusterMS = tgtSpectra(:,1);
for ii = 1:max(Id2Cluster)
    IdX =  Id2Cluster == ii;
    if sum(IdX) >= minSpectra
        ClusterMS(:, end+1) = mean(sqData(:, IdX), 2);
    end
end
% remove the MS with less than 5 non negative points
ClusterMS(:, [false, sum(ClusterMS(:, 2:end) ~= 0) <= 2]) = [];

%% first Model
% find peaks for each ClusterMS'~
lm = LocalMaxima(tgtSpectra, 3, 0.05*max(tgtSpectra(:,2)));
if isempty(lm)
    Model = [];
    FittData = [];
    return;
end

lm = sortrows(lm, -2);
iDx = findCloser(lm(1, 1), ClusterMS(:,1));
x(1) = idealVar; 
x(2) = ClusterMS(iDx, 1);

opts = optimset('Display','off', 'MaxFunEvals' , 5000);

[x, f] = fminsearch(@myFun, x, opts);
while 1
    [x_c, f_c] = fminsearch(@myFun, x, opts);
    if f_c <= ThresOpti*f
        x = x_c;
        f = f_c;
    else
        break
    end
end

%% Recursivly improve the model
x_i = x;
while 1
    [f_i, ~, Fitted] = myFun(x_i);
    Residuals =  smoothdata(mean(ClusterMS(:, 2:end)- Fitted, 2), 'movmean', PtsPerPeaks);
    [~, IDNewPeak] = max(Residuals);
    x(end+1) = ClusterMS(IDNewPeak, 1);
    x(1) = idealVar; 
    
    [x, f] = fminsearch(@myFun, x, opts);
    while 1
        [x_c, f_c] = fminsearch(@myFun, x, opts);
        
        if f_c <= ThresOpti*f
            x = x_c;
            f = f_c;
        else
            break
        end
    end
    
    Resolution = getResolution([x(2:end)' sqrt(x(1))*ones(size(x(2:end)'))]);
    if f <= Thres2*f_i && ~any(Resolution < Thres1, 'all')
        x_i = x;
    else
        break
    end
end
[SSR, Model] = myFun(x_i);
FittData = [x_i(2:end)' sqrt(x_i(1))*ones(size(x_i(2:end)'))];
warning on

    function [f, fittedModel, Fitted] = myFun(lm)
        
        Sz = length(lm) - 1;
        fittedModel = zeros(size(ClusterMS, 1), Sz(1));
        for mF_jj = 1:Sz(1)
            if lm(1) < minVar ...
                    || lm(1) > maxVar ...
                    || lm(mF_jj+1) < ClusterMS(1, 1)   - 2*sqrt(lm(1)) ...
                    || lm(mF_jj+1) > ClusterMS(end, 1) + 2*sqrt(lm(1))
                fittedModel(:, mF_jj) = inf(size(ClusterMS(:, 1)));
                
            else
                cModel = ModelSpectra;
                cModel(:, 1) = cModel(:,1)*sqrt(lm(1)) + lm(mF_jj+1);
                fittedModel(:, mF_jj) = interp1(cModel(:,1), cModel(:,2), ClusterMS(:, 1));
            end
        end
        
        fittedModel(isnan(fittedModel)) = 0;
        fittedModel(:, sum(fittedModel ~= 0) <= 2) = inf(size( fittedModel(:, sum(fittedModel ~= 0) <= 2)));
        Fitted = ((ClusterMS(:, 2:end)'/fittedModel')*fittedModel')';
        Fitted(Fitted < 0) = Fitted(Fitted < 0)*WeightNeg;
        
        f = sum(sum((ClusterMS(:,2:end) - Fitted).^2));
        if ~isfinite(f)
            f = 3*sum(sum((ClusterMS(:,2:end)).^2));
        end
        
    end

end




