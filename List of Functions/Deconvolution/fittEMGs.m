function model = fittEMGs(Data, n)

while 1
    [~, H] = nnmf(Data(:, 2:end)', n, 'Replicates', 10);
    if any(sum(H, 2) == 0)
        if n - 1 == 0
            break
        else
            n = n - 1;
        end
    else
        break
    end
end

XY = Data(:, 1);
XY(:, 2) = mean(Data(:, 2:end), 2);
opts = optimset('MaxFunEvals', 1e5, 'MaxIter', 1e5);


for ii = 1:n
    cXY = [XY(:,1), H(ii, :)'];
    
    M   = ChrMoment(cXY);
    [~, IdMax] = max(cXY(:,2));
    x0(3*ii-2) = XY(IdMax, 2);
    x0(3*ii-1) = XY(IdMax, 1);
    x0(3*ii) = sqrt(M(3))/n;
end

SumModel = [];
[f_Gau, model] = myfunGauss(x0);
x = fminsearch(@myfunGauss, x0);
[f_Gau, model] = myfunGauss(x);

x_emg_plus = reshape(x, [3 n]);
if n > 1
    x_emg_plus(3, 2:end) = ones(1, n-1)*x_emg_plus(3, 1);
end
x_emg_plus(end+1, :) = 0;
x_emg_plus = x_emg_plus(:);
[f_emg_plus, model_emg_plus] = myfun(x_emg_plus);
x_emg_plus = fminunc(@myfun, x_emg_plus);
[f_emg_plus, model_emg_plus] = myfun(x_emg_plus);


SumModel = [];
[f, model] = myfun(x0);
x = fminsearch(@myfun, x0);
[f, model] = myfun(x);

if n > 1
    x = reshape(x, 4, n);
    x(3, :) = ones(1, n)*mean(x(3,:));
    x(4, :) = ones(1, n)*mean(x(4,:));
    x = x(:);
    x = fminsearch(@myfun, x0);
    [f, model] = myfun(x);
end

    function [f, model] = myfunGauss(x)
        nfct = length(x)/3;
        SumModel = zeros(size(XY, 1), 1);
        for jj = 1:nfct
            if any(x(3*jj-2:3*jj-1) < 0)
                model(:,jj) = zeros(size(XY, 1), 1);
            else
                model(:,jj) = gaussPeak(XY(:,1), x([3*jj-1 3]));
            end
            SumModel = SumModel + x(3*jj-2)*model(:,jj);
        end
        
        f = sqrt(sum((SumModel - XY(:,2)).^2));
    end

    function [f, model] = myfun(x)
        nfct = length(x)/4;
        SumModel = zeros(size(XY, 1), 1);
        for jj = 1:nfct
            if any(x(4*jj-3:4*jj-1) < 0)
                model(:,jj) = zeros(size(XY, 1), 1);
            else
                model(:,jj) = EMGPeak(XY(:,1), x([4*jj-2 3 4]));
            end
            SumModel = SumModel + x(4*jj-3)*model(:,jj);
        end
        
        f = sqrt(sum((SumModel - XY(:,2)).^2));
    end
end
