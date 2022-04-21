function [Model, PeakModels] = fittPMG(Data, constrained)

ThresInt = 0.05;
opts = optimset('MaxFunEvals', 1e6, 'MaxIter', 1e6, ...
    'TolX', 1e-6, 'TolFun', 1e-6, 'Display','off');
minResolution = 0.5;
Wdw = 5;
%%% GAUSSIAN FITTING
% 1. Find Local maximum and peak variance
XY = [Data(:, 1), mean(Data(:, 2:end), 2)];
MaxInt = max(XY(:,2));
lm = LocalMaxima(XY, Wdw, ThresInt*MaxInt);
x = [];

for ii = 1:size(lm, 1)
    x(1, ii) =  lm(ii, 1);
    % Find peak variance
    Id4Max = findCloser(x(1, ii), XY(:,1));
    Val4Max = XY(Id4Max, 2);
    is = find(XY(1:Id4Max, 2) <= Val4Max/2, 1, 'last');
    if isempty(is), is = 1; end
    
    ie = find(XY(Id4Max:end, 2) <= Val4Max/2, 1, 'first');
    if isempty(ie)
        ie = length(XY(:, 2));
    else
        ie = min(ie + Id4Max, length(XY(:, 2)));
    end
    
    x(2, ii) =  (XY(ie, 1) - XY(is, 1))/2.355;
    x(4, ii) =  Val4Max;
end

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
    [x, f] = fminsearch(@FitManyFunctions, x_old, opts);
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
Model = zeros(size(XY, 1), size(x, 2));
for ii = 1:size(x, 2)
    Model(:, ii) = PMG1Peak(XY(:,1), x(1:3, ii));
    PeakModels{ii}.Function = 'PMG1';
    PeakModels{ii}.FittingParameters = x(1:3, ii);
end
Model = Model./max(Model);


    function [f, model] = FitSingleGauss(x)
        if any(x([1 2 4]) < 0)
            model = inf(size(XY, 1), 1);
        else
            model = x(4)*gaussPeak(XY(:,1), x(1:3));
        end
        f = sqrt(sum((model - XY(:,2)).^2));
    end

    function [f, FittedFunctions] = FitManyGauss(x)
        FittedFunctions = zeros(size(XY(:,2)));
        for f_ii = 1:size(x, 2)
            if any(x([1:2 4], f_ii) < 0 )
                FittedFunctions = inf(size(FittedFunctions));
                
            else
                FittedFunctions = FittedFunctions + (x(4, f_ii)*...
                    gaussPeak(XY(:,1), x(1:2, f_ii)));
                FittedFunctions(FittedFunctions < 0) = 0;
            end
        end
        
        f = sqrt(sum((XY(:,2) - FittedFunctions).^2));
    end

    function [f, FittedFunctions] = FitManyFunctions(x)
        FittedFunctions = zeros(size(XY(:,2)));
        for f_ii = 1:size(x, 2)
            if any(x([1:2 4], f_ii) < 0 )
                FittedFunctions = inf(size(FittedFunctions));
                
            else
                FittedFunctions = FittedFunctions + (x(4, f_ii)*...
                    PMG1Peak(XY(:,1), x(1:3, f_ii)));
                FittedFunctions(FittedFunctions < 0) = 0;
            end
        end
        
        f = sqrt(sum((XY(:,2) - FittedFunctions).^2));
    end

    function [f, FittedFunctions] = FitManyFunctions_const(x)
        FittedFunctions = zeros(size(XY(:,2)));
        for f_ii = 1:size(x, 2)
            if any(x([1:2 4], f_ii) < 0 )
                FittedFunctions = inf(size(FittedFunctions));
                
            else
                FittedFunctions = FittedFunctions + (x(4, f_ii)*...
                    PMG1Peak(XY(:,1), [x(1:2, f_ii); x(3, 1)]));
                FittedFunctions(FittedFunctions < 0) = 0;
            end
        end
        
        f = sqrt(sum((XY(:,2) - FittedFunctions).^2));
    end
end
