function dataOut = PMG1Peak(axe, a)
%http://dx.doi.org/10.1016/j.aca.2012.10.035

dataOut = exp(-0.5*(((axe-a(1)).^2)./((a(2) + a(3)*(axe-a(1))).^2)));
dataOut(dataOut < 0) = 0;

end

