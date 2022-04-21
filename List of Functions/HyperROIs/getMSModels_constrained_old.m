function [Model, FittData] = getMSModels_constrained(ModelSpectra, tgtSpectra, Data, p)

%% Find model spectra and prediction for MS spectra
ThresOpti    = 0.95;
ThresOpti2   = 0.95;
PearsonThresholdClustering = 0.25; % Dissimilarity index (1-CC)
minSpectra = 5;
PtsPerPeaks  = 10;
ThresholdResolution = 0.5;
%%%%%

%% Find spectra to fit
minVar = (mean(tgtSpectra(2:end, 1) - tgtSpectra(1:end-1, 1))/2)^2;

sqData = reshape(Data, size(Data, 1), []);
sqData(:, sum(sqData ~= 0, 1) <= 1) = [];

Z = linkage(sqData', 'complete', 'correlation');

Id2Cluster = cluster(Z, 'cutoff', PearsonThresholdClustering, 'Criterion','distance');
ClusterMS = tgtSpectra(:,1);
for ii = 1:max(Id2Cluster)
    IdX =  Id2Cluster == ii;
    if sum(IdX) >= minSpectra
        ClusterMS(:, end+1) = sum(sqData(:, IdX), 2);
    end
end
% remove the MS with less than 5 non negative points
ClusterMS(:, [false, sum(ClusterMS(:, 2:end) ~= 0) <= 2]) = [];

%% first Model
% find peaks for each ClusterMS'~
lm = LocalMaxima(tgtSpectra, 2, 0);
if isempty(lm)
    Model = [];
    FittData = [];
    return;
end

lm = sortrows(lm, -2);
iDx = findCloser(lm(1, 1), ClusterMS(:,1));
x(1, 1) = ClusterMS(iDx, 1);
x(1, 2) = max(minVar, polyval(p{1},  x(1, 1), p{2}, p{3}));

opts = optimset('Display','off', 'MaxFunEvals' , 5000);

[x, f] = fminsearch(@myFun, x, opts);

while 1
    [x_c, f_c] = fminsearch(@myFun, x, opts);
    if any(x_c(:, 2) <= 0)
        x_c(x_c(:, 2) <= 0, :) = [];
        [x_c, f_c] = fminsearch(@myFun, x_c, opts);
    end
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
    Residuals =  smoothdata(sum(ClusterMS(:, 2:end)- Fitted, 2), 'movmean', PtsPerPeaks);
    [~, IDNewPeak] = max(Residuals);
    x(end+1, 1) = ClusterMS(IDNewPeak, 1);
    x(end, 2) = max(minVar, polyval(p{1},  x(end, 1), p{2}, p{3}));
    
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
    
    if f <= ThresOpti2*f_i %&& ~any(Resolution < ThresholdResolution, 'all')
        x_i = x;
    else
        break
    end
    
    if size(x, 1) > size(tgtSpectra, 1)/3
        break
    end
    
end
[~, Model] = myFun(x_i);

%% LAST STEP: CHECK FOR DEGREE OF SIMILARITY AND MERGE


isMerged = false(size(x, 1), 1);
FittData = [x'; isMerged'];

    function [f, fittedModel, Fitted] = myFun(lm)
        WeightNeg = 100;
        
        Sz = size(lm) - [0 2];
        fittedModel = zeros(size(ClusterMS, 1), Sz(1));
        Resolution = getResolution([lm(:,1) sqrt(lm(:,2))]);
        for mF_jj = 1:Sz(1)
            [Var_est, delta] = polyval(p{1}, lm(1, 1), p{2}, p{3});
            if lm(1, 2) < Var_est - 3*delta ...
                    || lm(1, 2) > Var_est + 3*delta ...
                    || lm(1, 2) < minVar ...
                    || lm(mF_jj, 1) < ClusterMS(1, 1) - 4*sqrt(lm(1, 2)) ...
                    || lm(mF_jj, 1) > ClusterMS(end, 1) + 4*sqrt(lm(1, 2)) ...
                    || any(Resolution < ThresholdResolution, 'all')
                fittedModel(:, mF_jj) = inf(size(ClusterMS(:, 1)));
                
            else
                cModel = ModelSpectra;
                cModel(:, 1) = cModel(:,1)*sqrt(lm(1,2)) + lm(mF_jj,1);
                fittedModel(:, mF_jj) = interp1(cModel(:,1), cModel(:,2), ClusterMS(:, 1));
            end
        end
        
        lm(:, 2) = lm(1, 2);
        fittedModel(isnan(fittedModel)) = 0;
        fittedModel(:, sum(fittedModel ~= 0) <= 2) = inf(size( fittedModel(:, sum(fittedModel ~= 0) <= 2)));
        Fitted = ((ClusterMS(:, 2:end)'/fittedModel')*fittedModel')';
        Fitted(Fitted < 0) = Fitted(Fitted < 0)*WeightNeg;
        
        f = sqrt(sum(sum((ClusterMS(:,2:end) - Fitted).^2)));
    end

end




