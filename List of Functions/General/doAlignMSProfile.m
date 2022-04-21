function [shift, Coef] = doAlignMSProfile(XY4ref, XY2alig, Thrs_mz, options)

if nargin == 3
    options = optimset('Display', 'off');
end

if size(XY4ref, 1) < size(XY4ref, 2)
    XY4ref = XY4ref';
end

if size(XY2alig, 1) < size(XY2alig, 2)
    XY2alig = XY2alig';
end

warning off

x0 = XY4ref(XY4ref(:,2) == max(XY4ref(:,2)), 1) - ...
    XY2alig(XY2alig(:,2) == max(XY2alig(:,2)), 1);
mzAve = sum(XY4ref(:,1).*XY4ref(:,2))/sum(XY4ref(:,2));

[shift, Coef] = fmincon(@fun, x0, 0, 0, [], [], ...
    -mzAve*Thrs_mz/1000000, mzAve*Thrs_mz/1000000, [],  options);

warning on

    function f = fun(x)
        vq = interp1(XY2alig(:,1) + x, XY2alig(:,2), XY4ref(:,1));
        vq(isnan(vq)) = 0;
        f = pdist2( XY4ref(:,2)', vq', 'correlation');
    end
end
