function [Model, PeakModels] = fittPMG_cleaned(Data, ThresholdResolution, PtsPerPeaks, ThresWhile, maxPeaks, Optimised, MaxTail, minMaxVar)

if nargin == 4
    maxPeaks = 10;
    Optimised = false;
    MaxTail = 0.1;
    minMaxVar = 100;
end

maxLoop = 5;
Data(isnan(Data)) = 0;
max_nnmf = 5;

% NNMF NOT CONCLUSIVE
% opt = statset('MaxIter',5,'Display','off');
% for ii = 1:max_nnmf
% [W0{ii}, H0{ii}, D0{ii}] = nnmf(Data(:, 2:end), ii, 'Replicates', 10,...
%                    'Options',opt,...
%                    'Algorithm','mult');
% end

% TrY with bottom up divisive HCA
% [IdSplit, PearSplit] = SplitMe(Data(:, 2:end))

opts = optimset('MaxFunEvals', 1e4, 'MaxIter', 1e4, ...
    'TolX', 1e-5, 'TolFun', 1e-5, 'Display','off');

% 1. Find Local maximum and peak variance
XY = [Data(:, 1), mean(Data(:, 2:end), 2, 'omitnan')];
XY(isnan(XY)) = 0;
x = [];
FittedModel = zeros(size(XY(:,2)));

Residuals = XY(:,2) - FittedModel;
[Val4Max, Id4Max] = max(Residuals);
x(1, end+1) =  XY(Id4Max, 1);
is = find(Residuals(1:Id4Max) <= Val4Max/2, 1, 'last');
if isempty(is), is = 1; end

ie = find(Residuals(Id4Max:end) <= Val4Max/2, 1, 'first');
if isempty(ie)
    ie = length(XY(:, 2));
else
    ie = min(ie + Id4Max, length(XY(:, 2)));
end

MinVar = (mean(XY(2:end, 1) - XY(1:end-1, 1)));
MaxVar = minMaxVar*MinVar;


x(2, end) =  (XY(ie, 1) - XY(is, 1))/2.355;
x(4, end) =  Val4Max;
x(5, end) =  prctile(XY(:,2), 10);

x(2, x(2, :) < MinVar) = MinVar;
[f_ini, FittedModel] = FitManyFunctions(x);
[x, fval, exitflag, PeakModels] = fminsearch(@FitManyFunctions, x, opts);

count = 1;
while exitflag == 0
    [x, fval, exitflag, PeakModels] = fminsearch(@FitManyFunctions, x, opts);
    if count > maxLoop, break; else, count = count +1; end
end
[f_ini, FittedModel] = FitManyFunctions(x);

x_ini = x;
while 1
    
    if size(x_ini, 2) >= maxPeaks
        break
    end
    
    Residuals = smoothdata(XY(:,2) - FittedModel, 'movmean', PtsPerPeaks);
    [Val4Max, Id4Max] = max(Residuals);
    x(1, end+1) =  XY(Id4Max, 1);
    is = find(Residuals(1:Id4Max) <= Val4Max/2, 1, 'last');
    if isempty(is), is = 1; end
    
    ie = find(Residuals(Id4Max:end) <= Val4Max/2, 1, 'first');
    if isempty(ie)
        ie = length(XY(:, 2));
    else
        ie = min(ie + Id4Max, length(XY(:, 2)));
    end
    
    x(2, end) =  (XY(ie, 1) - XY(is, 1))/2.355;
    x(2, x(2, :) < MinVar) = MinVar;
    x(2, x(2, :) > MaxVar) = MaxVar;
    x(4, end) =  Val4Max;
    x(5, end) =  0;
    
    [x, fval, exitflag, PeakModels] = fminsearch(@FitManyFunctions, x, opts);
    count = 1;
    while exitflag == 0
        [x, fval, exitflag, PeakModels] = fminsearch(@FitManyFunctions, x, opts);
        if count > maxLoop, break; else, count = count +1; end
    end
    
    [f, FittedModel] = FitManyFunctions(x);
    Resolution = getResolution([x(1, :)' x(2, :)']);
    if f <= ThresWhile*f_ini && ~any(Resolution < ThresholdResolution, 'all')
        f_ini = f;
        x_ini = x;
    else
        break
    end
end

if Optimised
    x_new = x_ini;
    x_new(2, :) = ones(size(x_new(2, :)))*mean(x_new(2, :));
    x_new(3, :) = ones(size(x_new(3, :)))*mean(x_new(3, :));
    x_new(5, :) = zeros(size(x_new(5, :)));
    
    [x_new, fval, exitflag, PeakModels] = fminsearch(@FitManyFunctions, x_new, opts);
    count = 1;
    while exitflag == 0
        [x_new, fval, exitflag, PeakModels] = fminsearch(@FitManyFunctions, x_new, opts);
        if count > maxLoop, break; else, count = count +1; end
    end
    [f_new, FittedModel] = FitManyFunctions(x_new);
    
    if f_new < f_ini
        x_ini = x_new;
    end
end

[f, FittedModel] = FitManyFunctions(x_ini);
Model = zeros(size(XY, 1), size(x_ini, 2));
PeakModels = {};
for ii = 1:size(x_ini, 2)
    Model(:, ii) = PMG1Peak(XY(:,1), x_ini(1:3, ii));
    PeakModels{ii}.Function = 'PMG1';
    PeakModels{ii}.FittingParameters = x_ini(1:3, ii);
end
Model(isnan(Model)) = 0;

Id = sum(Model ~= 0) <= 5;
PeakModels(Id) = [];
Model(:, Id) = [];

    function [f, FittedFunctions] = FitManyFunctions(x)
        FittedFunctions = zeros(size(XY(:,2)));
        for f_ii = 1:size(x, 2)
            if any(x([1:2 4], f_ii) < 0 ) ...
                    || x(2, f_ii) < MinVar ...
                    || x(2, f_ii) > MaxVar ...
                    || x(1, f_ii) < XY(1, 1) - 3*x(2, f_ii) ...
                    || x(1, f_ii) > XY(end, 1) + 3*x(2, f_ii) ...
                    || abs(x(3, f_ii)) > MaxTail
                FittedFunctions = inf(size(FittedFunctions));
                
            else
                FittedFunctions = FittedFunctions + x(5, f_ii) + (x(4, f_ii)*...
                    PMG1Peak(XY(:,1), x(1:3, f_ii)));
                FittedFunctions(FittedFunctions < 0) = 0;
            end
        end
        
        f = sqrt(sum((XY(:,2) - FittedFunctions).^2));
        if ~isfinite(f)
            f = 3*sqrt(sum((XY(:,2)).^2));
        end
    end
end
