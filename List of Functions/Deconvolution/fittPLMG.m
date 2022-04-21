function Model = fittPLMG(Data)

ThresCC = 0.80;
n = 1;
minRs = 0.5;
optnmf = statset('MaxIter', 1000, 'Display','off', ...
    'TolFun', 1e-4, 'TolX', 1e-4);
[~, H_ini] = nnmf(Data(:, 2:end)', n, 'Replicates', 25,...
    'Options',optnmf,...
    'Algorithm','mult');
sumData = sum(Data, 2);
D_ini = sqrt(sum((sumData - ((sumData'/H_ini)*H_ini)').^2));
opts = optimset('MaxFunEvals', 1e6, 'MaxIter', 1e6, 'Display','iter',...
    'TolX',1e-15, 'TolFun', 1e-15);
while 1
    [~, H] = nnmf(Data(:, 2:end)', n+1, 'Replicates', 25,...
        'Options',optnmf,...
        'Algorithm','mult');
    D_test = sqrt(sum((sumData - ((sumData'/H)*H)').^2));
    
    if D_test <= 0.95*D_ini
        H_ini = H;
        D_ini = D_test;
        n = n+1;
    else
        break
    end
    
end

% Check correlation
if size(H_ini, 1) > 1
    CC = corrcoef(H_ini');
    while any(triu(CC, 1) >= ThresCC, 'all')
        [myLin, myCol] = find(triu(CC, 1) == max(max(triu(CC, 1))));
        H_ini(myLin(1), :) = H_ini(myLin(1), :) +  H_ini(myCol(1), :);
        H_ini(myCol(1), :) = [];
        
        CC = corrcoef(H_ini');
    end
end


if size(H_ini, 1) == 1
    model = H_ini';
    f = D_ini;
    return
end

x0 = zeros(6+size(H_ini, 1), size(H_ini, 1));
for ii = 1:size(H_ini, 1)
    XY = [Data(:,1), H_ini(ii, :)'];
    M  = ChrMoment(XY);
    X0(1) = M(2);
    X0(2) = sqrt(M(3));
    X0(3) = 0;
    X0(4) = 0;
    X0(5) = 0;
    X0(6) = 0;
    X0(7) = max(XY(:,2));
    [X, f_0(ii)] = fminsearch(@FitSingleFunction, X0, opts);
    
    while 1
        [Xf, fval] = fminsearch(@FitSingleFunction, X, opts);
        if fval <= 0.95*f_0
            X = Xf;
            f_0(ii) = fval;
        else
            break
        end
    end
    x0(1:6, ii) = X(1:6);
    x0(7+ii-1, ii) = X(7);
end

[~, Id] = min(f_0);
x0(2:6, :) = ones(size(x0(2:6, :))).*x0(2:6, Id);
[xf, fval] = fminsearch(@FitManyFunctions, x0, opts);
[f_init, model] = FitManyFunctions(xf);

[f_init, model] = FitManyFunctions(x0);
x0(3, :) = ones(size(x0(3, :)))*mean(x0(3, :));
x0(4, :) = ones(size(x0(4, :)))*mean(x0(4, :));
x0(5, :) = ones(size(x0(5, :)))*mean(x0(5, :));
x0(6, :) = ones(size(x0(6, :)))*mean(x0(6, :));
[x, f_0] = fminsearch(@FitManyFunctions, x0, opts);
while 1
    x(3, :) = ones(size(x(3, :)))*mean(x(3, :));
    x(4, :) = ones(size(x(4, :)))*mean(x(4, :));
    x(5, :) = ones(size(x(5, :)))*mean(x(5, :));
    x(6, :) = ones(size(x(6, :)))*mean(x(6, :));
    [xf, fval] = fminsearch(@FitManyFunctions, x, opts);
    if fval <= 0.95*f_0
        x = xf;
        f_0 = fval;
    else
        break
    end
end

[f_init, model] = FitManyFunctions(x);
if size(x, 2) > 1
    while 1
        clear f_all x_new
        
        for ii = 1:size(x, 2)
            [x_new{ii}, f_all(ii)] = fminsearch(@FitManyFunctions, x(:, 1:end~=ii), opts);
            
            while 1
                [xf, fval] = fminsearch(@FitManyFunctions, x_new{ii}, opts);
                if fval <= 0.95*f_all(ii)
                    x_new{ii} = xf;
                    f_all(ii) = fval;
                else
                    break
                end
            end
        end
        
        if any(f_all <= 1.10*f_init)
            [~, Id2Keep] = min(f_all);
            x = x_new{Id2Keep};
            f_init = f_all(Id2Keep);
            if size(x, 2) == 1, break; end
        else
            break
        end
    end
end

% check for minimum Rs
if size(x, 2) > 1
    while 1
        Model = zeros(size(XY, 1), size(x, 2));
        clear M
        for ii = 1:size(x, 2)
            Model(:, ii) = PMG2Peak(XY(:,1), x(1:4, ii));
            Model(Model < 0) = 0;
            M(ii, :) = ChrMoment([XY(:,1), Model(:, ii)]);
        end
        [M, Ix] = sortrows(M, 2);
        x = x(:, Ix);
        M(:, 5) = NaN;
        M(2:end, 5) = ...
            (M(2:end, 2) - M(1:end-1, 2))./(2*(sqrt(M(2:end, 3)) + sqrt(M(1:end-1, 3))));
        
        if any(M(:,5) < minRs)
            [~, Id2Rem] = min(M(:,5));
            x(:, Id2Rem) = [];
            [x, f_0] = fminsearch(@FitManyFunctions, x, opts);
            while 1
                [xf, fval] = fminsearch(@FitManyFunctions, x, opts);
                if fval <= 0.95*f_0
                    x = xf;
                    f_0 = fval;
                else
                    break
                end
            end
            if size(x, 2) == 1, break; end
        else
            break
        end
    end
end

Model = zeros(size(XY, 1), size(x, 2));
for ii = 1:size(x, 2)
    Model(:, ii) = PMG2Peak(XY(:,1), x(1:4, ii));
end
Model = Model./max(Model);


    function [f, model] = FitSingleFunction(x)
        if any(x([1 2 7]) < 0) || x(4) < (x(3)^2)/4
            model = inf(size(XY, 1), 1);
        else
            model = x(7)*PLMGPeak(XY(:,1), x(1:6));
        end
        f = sqrt(sum((model - XY(:,2)).^2));
    end

    function [f, FittedFunctions] = FitManyFunctions(x)
        FittedFunctions = zeros(size(H_ini));
        for f_ii = 1:size(H_ini, 1)
            for f_jj = 1:size(x, 2)
                if any([x(6+f_ii, f_jj);  x(1:2, f_jj)] < 0) ||...
                        x(4, f_jj) < (x(3, f_jj)^2)/4
                    FittedFunctions(f_ii, :)  = inf(size(FittedFunctions(f_ii, :)));
                else
                    FittedFunctions(f_ii, :) = FittedFunctions(f_ii, :) + (x(6+f_ii, f_jj)*...
                        PLMGPeak(XY(:,1), x(1:6, f_jj)))';
                    FittedFunctions(FittedFunctions < 0) = 0;
                end
            end
        end
        f = sqrt(sum(sum((H_ini-FittedFunctions)'.^2)));
    end
end
