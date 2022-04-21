function [Aligned_Data, Master, Shift, Coef] = doAlignX2Deconv(Data, leeway, options)

Thres4Ali    = 0.01;
Master       = mean(Data, 1);
[Lin, Col]   = size(Data);
AxeX         = 1:Col;
Aligned_Data = []; 
opts = optimset('Display', 'off');

for ii = 1:Lin
    [Shift(ii, 1), Coef(ii, 1)] = fmincon(@fun, 0.1, 0, 0, [], [], ...
        -leeway, leeway, [], opts);
    [~, vq] = fun(Shift(ii, 1));
    Aligned_Data(ii, :) = vq;
end

Master       = mean(Aligned_Data, 1);
for ii = 1:Lin
    [Shift(ii, 2), Coef(ii, 2)] = fmincon(@fun, Shift(ii, 1), 0, 0, [], [], ...
        -leeway, leeway, [], opts);
    [~, vq] = fun(Shift(ii, 2));
    Aligned_Data(ii, :) = vq;
end

    function [f, vq] = fun(x)
        vq = interp1(AxeX + x, Data(ii, :), AxeX);
        vq(isnan(vq)) = 0;
        f = pdist2( Master, vq, 'correlation');
    end

end

