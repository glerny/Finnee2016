function [Model, PeakModels] = fittPMG(Data, constrained)

ThresPCA = 95;
maxPCAPeaks = 5;

% Find nbr of peaks using pca
[S, H, ~, ~, explained] = pca(Data(:, 2:end)', 'algorithm','als');
Npeaks =  min(5, find(cumsum(explained) >= ThresPCA, 1, 'first'));

optnmf = statset('MaxIter', 1000, 'Display','off', ...
    'TolFun', 1e-4, 'TolX', 1e-4);
[~, H_ini] = nnmf(Data(:, 2:end)', Npeaks, 'Replicates', 25,...
    'Options',optnmf,...
    'Algorithm','mult');

opts = optimset('MaxFunEvals', 1e6, 'MaxIter', 1e6, ...
    'TolX', 1e-6, 'TolFun', 1e-6, 'Display','off');

for ii = 1:size(H_ini, 1)
    XY = [Data(:,1), H_ini(ii, :)'];
    M  = ChrMoment(XY);
    X0(1) = M(2);
    X0(2) = sqrt(M(3));
    X0(3) = 0;
    X0(4) = max(XY(:,2));
    [X{ii}, f_0(ii)] = fminsearch(@FitSingleGauss, X0, opts);
    
    while 1
        [Xf, fval] = fminsearch(@FitSingleGauss, X{ii}, opts);
        if fval <= 0.95*f_0
            X{ii} = Xf;
            f_0(ii) = fval;
        else
            break
        end
    end
end
[~, Id] = min(f_0);
x0 = X{Id}';

x0 = [x0; zeros(size(H_ini, 1)-1, 1)];
[x, f] = fminsearch(@FitManyGauss, x0, opts);
while 1
    [xf, fval] = fminsearch(@FitManyGauss, x, opts);
    if fval <= 0.95*f
        x = xf;
        f = fval;
    else
        break
    end
end

[f_init, model] = FitManyGauss(x);
x_old = x;

while 1
    Residual = (H_ini-model);
    Residual(Residual <0) = 0;
    [~, I2] = max(max((Residual)));
    x(1, end+1) = XY(I2, 1);
    x(2:3, end) = X{Id}(2:3)';
    x(4:end, end) = Residual(:, I2);
    x(1:3, 1) = X{Id}(1:3)';
    [x, f] = fminsearch(@FitManyGauss, x, opts);
    while 1
        [xf, fval] = fminsearch(@FitManyGauss, x, opts);
        if fval <= 0.95*f
            x = xf;
            f = fval;
        else
            break
        end
    end
    
    if f <= f_init && ...
            any(x (4:end, end) >= 0.01*max(x(4:end, 1:end-1)), 'all') && ...
            all(abs(x(1, end) - x(1, 1:end-1))./(2*(x(2, end) + x(2, 1:end-1))) >= 0.5) && ...
            all(x(1, :) > XY(1, 1) - 3*x0(2)) && ...
            all(x(1, :) < XY(end, 1) + 3*x0(2)) && ...
            size(x, 2) < 10
             
            x = xf;
            f = fval;
        f_init = f;
        x_old = x;
        [~, model] = FitManyGauss(x);
    else
        x = x_old;
        break
    end
end

if constrained
    x_old(3, :) = zeros(size(x_old(3, :)));
    [x, f] = fminsearch(@FitManyFunctions_const, x_old, opts);
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
    x_old(3, :) = zeros(size(x_old(3, :)));
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
        FittedFunctions = zeros(size(H_ini));
        for f_ii = 1:size(H_ini, 1)
            for f_jj = 1:size(x, 2)
                if any([x(3+f_ii, f_jj);  x(1:2, f_jj)] < 0)
                    FittedFunctions(f_ii, :)  = inf(size(FittedFunctions(f_ii, :)));
                else
                    FittedFunctions(f_ii, :) = FittedFunctions(f_ii, :) + (x(3+f_ii, f_jj)*...
                        gaussPeak(XY(:,1), x(1:2, f_jj)))';
                    FittedFunctions(FittedFunctions < 0) = 0;
                end
            end
        end
        f = sqrt(sum(sum((H_ini-FittedFunctions)'.^2)));
    end

    function [f, FittedFunctions] = FitManyFunctions(x)
        FittedFunctions = zeros(size(H_ini));
        for f_ii = 1:size(H_ini, 1)
            for f_jj = 1:size(x, 2)
                if any([x(3+f_ii, f_jj);  x(1:2, f_jj)] < 0)
                    FittedFunctions(f_ii, :)  = inf(size(FittedFunctions(f_ii, :)));
                else
                    FittedFunctions(f_ii, :) = FittedFunctions(f_ii, :) + (x(3+f_ii, f_jj)*...
                        PMG1Peak(XY(:,1), x(1:3, f_jj)))';
                    FittedFunctions(FittedFunctions < 0) = 0;
                end
            end
        end
        f = sqrt(sum(sum((H_ini-FittedFunctions)'.^2)));
    end

    function [f, FittedFunctions] = FitManyFunctions_const(x)
        FittedFunctions = zeros(size(H_ini));
        for f_ii = 1:size(H_ini, 1)
            for f_jj = 1:size(x, 2)
                if any([x(3+f_ii, f_jj);  x(1:2, f_jj)] < 0)
                    FittedFunctions(f_ii, :)  = inf(size(FittedFunctions(f_ii, :)));
                else
                    FittedFunctions(f_ii, :) = FittedFunctions(f_ii, :) + (x(3+f_ii, f_jj)*...
                        PMG1Peak(XY(:,1), [x(1:2, f_jj); x(3, 1)]))';
                    FittedFunctions(FittedFunctions < 0) = 0;
                end
            end
        end
        f = sqrt(sum(sum((H_ini-FittedFunctions)'.^2)));
    end
end
