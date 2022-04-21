function [Model, FittData] = getMSModels_constrained(ModelSpectra, Data, MZAxis, p)

%% Find model spectra and prediction for MS spectra
ThresOpti                  = 0.95;
Thres1                     = 0.3;
Thres2                     = 0.95;
PtsPerPeaks                = 10;
%%%%%

warning off
%% Find spectra to fit
[y,delta] = polyval(p{1}, MZAxis(1), p{2}, p{3});
minVar   = max(0, y-delta);
[y,delta] = polyval(p{1}, MZAxis(end), p{2}, p{3});
maxVar   = y+delta;
idealVar = (minVar+maxVar)/2;
tgtSpectra = [MZAxis, mean(mean(Data, 2), 3)];

%% first Model
% find peaks for each ClusterMS'~
lm = LocalMaxima(tgtSpectra, 3, 0.05*max(tgtSpectra(:,2)));
if isempty(lm)
    Model = [];
    FittData = [];
    return;
end

lm = sortrows(lm, -2);
iDx = findCloser(lm(1, 1), tgtSpectra(:,1));
x(1) = idealVar; 
x(2) = tgtSpectra(iDx, 1);

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
    [f_i, fittedModel, Fitted] = myFun(x_i);
    Residuals =  smoothdata(mean(tgtSpectra(:, 2)- Fitted, 2), 'movmean', PtsPerPeaks);
    [~, IDNewPeak] = max(Residuals);
    x(end+1) = tgtSpectra(IDNewPeak, 1);
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
FittData = [x_i(2:end)' x_i(1)*ones(size(x_i(2:end)'))];
warning on

    function [f, fittedModel, Fitted] = myFun(lm)
        
        Sz = length(lm) - 1;
        fittedModel = zeros(size(tgtSpectra, 1), Sz(1));
        for mF_jj = 1:Sz(1)
            if lm(1) < minVar...
                    || lm(1) > maxVar ...
                    || lm(mF_jj+1) < tgtSpectra(1, 1)   - 3*sqrt(lm(1)) ...
                    || lm(mF_jj+1) > tgtSpectra(end, 1) + 3*sqrt(lm(1))
                fittedModel(:, mF_jj) = inf(size(tgtSpectra(:, 1)));
                
            else
                cModel = ModelSpectra;
                cModel(:, 1) = cModel(:,1)*sqrt(lm(1)) + lm(mF_jj+1);
                fittedModel(:, mF_jj) = interp1(cModel(:,1), cModel(:,2), tgtSpectra(:, 1));
            end
        end
        
        fittedModel(isnan(fittedModel)) = 0;
        fittedModel(:, sum(fittedModel ~= 0) <= 2) = inf(size( fittedModel(:, sum(fittedModel ~= 0) <= 2)));
        Fitted = ((tgtSpectra(:, 2)'/fittedModel')*fittedModel')';
        
        f = sum((tgtSpectra(:, 2) - Fitted).^2);
        if ~isfinite(f)
            f = 3*sum(tgtSpectra(:, 2).^2);
        end
        
    end

end