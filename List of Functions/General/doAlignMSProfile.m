function [shift, Coef] = doAlignMSProfile(XY4ref, XY2alig, options)
if nargin == 2
    options = optimset('Display', 'off');
end

warning off

x0 = XY4ref(XY4ref(:,2) == max(XY4ref(:,2)), 1) - ...
    XY2alig(XY2alig(:,2) == max(XY2alig(:,2)), 1);
mzAve = sum(XY4ref(:,1).*XY4ref(:,2))/sum(XY4ref(:,2));

[shift, Coef, exitflag, output] = fmincon(@fun, x0, 0, 0, [], [], ...
    -mzAve*5/1000000, mzAve*5/1000000, [],  options);

warning on

    function f = fun(x)
        vq = interp1(XY2alig(:,1) + x, XY2alig(:,2), XY4ref(:,1));
        vq(isnan(vq)) = 0;
        f = pdist2( XY4ref(:,2)', vq', 'correlation');
    end
end
