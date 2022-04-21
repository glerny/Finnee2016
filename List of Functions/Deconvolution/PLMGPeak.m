function dataOut = PLMGPeak(axe, a)
%http://dx.doi.org/10.1016/j.aca.2012.10.035
% a(1) : tm
% a(2) : sigma0
% a(3) : w
% a(4) : z
% a(5) : b
% a(6) : c
% Conditions:
% a(2) > 0; 
% a(4) > a(3)^2/4

tc = (axe-a(1));
num = 1 + 100*a(3)*tc + 100*a(4)*(tc.^2);
den = a(2)^2+a(5)*tc + 100*a(6)*(tc.^2);

dataOut = exp(-0.5*(num./den).*tc.^2);
end

