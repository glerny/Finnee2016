function [shift, Coef, Ratio] = doAlignAllMSProfile(XY4ref, MSlist, Axis, Rp)

options = optimset('Display', 'off');
options = optimset('Display','iter','PlotFcns',@optimplotfval);;
x0(1) = 0;
if length(MSlist) > 1
    disp('dtt')
    for ii = 1:length(MSlist)
        IdM(ii) = MSlist{ii}(MSlist{1}(:,2) == max(MSlist{ii}(:,2)), 1);
    end
    
    I1 = XY4ref(findCloser(IdM(1), XY4ref(:,1)), 2);
    for ii = 2:length(IdM)
        x0(ii)  = XY4ref(findCloser(IdM(ii), XY4ref(:,1)), 2); 
        x0(ii)  = x0(ii)/I1;
    end
end

[shift, Coef, exitflag, output] = fminsearch(@fun, x0, options);

Ratio = shift(2:end);
shift = shift(1);

    function f = fun(x)
        XY     = Axis;
        for jj = 1:length(MSlist)
            MSc = MSlist{jj};
            MSc(:,1) = MSc(:,1) + x(1);
            if jj > 1
                xy = centr2prof(MSc, Rp, Axis);
                XY(:,2) = XY(:,2)+ x(jj)*xy(:, 2);
            else
                XY = centr2prof(MSc, Rp, Axis);
            end
        end
        x
        f = pdist2( XY4ref(:,2)', XY(:,2)' , 'correlation');
    end
end
