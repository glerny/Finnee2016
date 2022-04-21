function [Model, PeakModels, isMerged] = fittPMG(Data, constrained)

ThresWhile          = 0.95;
ThresholdResolution = 1.5;
MaxTail             = 0.2;
PtsPerPeaks         = 20;

opts = optimset('MaxFunEvals', 1e4, 'MaxIter', 1e4, ...
    'TolX', 1e-6, 'TolFun', 1e-6, 'Display','off');

% 1. Find Local maximum and peak variance
XY = [Data(:, 1), mean(Data(:, 2:end), 2)];
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
x(2, end) =  (XY(ie, 1) - XY(is, 1))/2.355;
x(4, end) =  Val4Max;
x(2, x(2, :) < MinVar) = MinVar;

if constrained
    [x, f] = fminsearch(@FitManyFunctions_const, x, opts);
    while 1
        [xf, fval] = fminsearch(@FitManyFunctions_const, x, opts);
        if fval <= 0.95*f
            x = xf;
            f = fval;
        else
            break
        end
    end
    x(3, :) = ones(size(x(3, :))).*x(3, 1);
    
else
    [x, f] = fminsearch(@FitManyFunctions, x, opts);
    while 1
        [xf, fval] = fminsearch(@FitManyFunctions, x, opts);
        if fval <= 0.95*f
            x = xf;
            f = fval;
        else
            break
        end
    end
end

[f_ini, FittedModel] = FitManyFunctions_const(x);
x_ini = x;
    
while 1
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
    x(4, end) =  Val4Max;
    
    if constrained
        [x, f] = fminsearch(@FitManyFunctions_const, x, opts);
        while 1
            [xf, fval] = fminsearch(@FitManyFunctions_const, x, opts);
            if fval == inf
                break
            end
            
            if fval <= 0.95*f
                x = xf;
                f = fval;
            else
                break
            end
        end
        x(3, :) = ones(size(x(3, :))).*x(3, 1);
        
    else
        [x, f] = fminsearch(@FitManyFunctions, x, opts);
        while 1
            [xf, fval] = fminsearch(@FitManyFunctions, x, opts);
            if fval <= 0.95*f
                x = xf;
                f = fval;
            else
                break
            end
        end
    end
    
    [f, FittedModel] = FitManyFunctions_const(x);
    
    if f <= ThresWhile*f_ini 
        f_ini = f;
        x_ini = x;
    else
        break
    end
    
    if size(x_ini, 2) > 10
        break
    end
end

if constrained
    x_ini(3, :) = x_ini(3, 1);
end

isMerged = false(size(x_ini, 2), 1);
Model = zeros(size(XY, 1), size(x_ini, 2));
for ii = 1:size(x_ini, 2)
    Model(:, ii) = PMG1Peak(XY(:,1), x_ini(1:3, ii));
    PeakModels{ii}.Function = 'PMG1';
    PeakModels{ii}.FittingParameters = x_ini(1:3, ii);
end
Model(isnan(Model)) = 0;

Id = sum(Model ~= 0) <= 5;
PeakModels(Id) = [];
Model(:, Id) = [];
isMerged(Id) = [];

    function [f, FittedFunctions] = FitManyFunctions(x)
        FittedFunctions = zeros(size(XY(:,2))); 
        Resolution = getResolution([x(1, :)' x(2, :)']);
        for f_ii = 1:size(x, 2)
            if any(x([1:2 4], f_ii) < 0 ) ...
                || x(2, f_ii) < MinVar ...
                || x(1, f_ii) < XY(1, 1) - 3*x(2, f_ii) ...
                || x(1, f_ii) > XY(end, 1) + 3*x(2, f_ii) ...
                || abs(x(3, f_ii)) > MaxTail ...
                || any(Resolution < ThresholdResolution, 'all')
                FittedFunctions = inf(size(FittedFunctions));
                
            else
                FittedFunctions = FittedFunctions + (x(4, f_ii)*...
                    PMG1Peak(XY(:,1), x(1:3, f_ii)));
                FittedFunctions(FittedFunctions < 0) = 0;
            end
        end
        
        f = sqrt(sum((XY(:,2) - FittedFunctions).^2));
        if ~isfinite(f)
            f = sqrt(sum((XY(:,2)).^2));
        end
    end

    function [f, FittedFunctions] = FitManyFunctions_const(x)
        FittedFunctions = zeros(size(XY(:,2)));
        Resolution = getResolution([x(1, :)' x(2, :)']);
        for f_ii = 1:size(x, 2)
            if any(x([1:2 4], f_ii) < 0 ) ...
                || x(2, f_ii) < MinVar ...
                || x(1, f_ii) < XY(1, 1) - 3*x(2, f_ii) ...
                || x(1, f_ii) > XY(end, 1) + 3*x(2, f_ii) ...
                || abs(x(3, 1)) > MaxTail ...
                || any(Resolution < ThresholdResolution, 'all')
                FittedFunctions = inf(size(FittedFunctions));
                
            else
                FittedFunctions = FittedFunctions + (x(4, f_ii)*...
                    PMG1Peak(XY(:,1), [x(1:2, f_ii); x(3, 1)]));
                FittedFunctions(FittedFunctions < 0) = 0;
            end
        end
        
        f = sqrt(sum((XY(:,2) - FittedFunctions).^2));
        if ~isfinite(f)
            f = sqrt(sum((XY(:,2)).^2));
        end
    end
end
