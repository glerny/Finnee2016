function dataOut = gaussPeak(axe, a)

dataOut = exp(-(axe-a(1)).^2/(2*a(2)^2));


end

