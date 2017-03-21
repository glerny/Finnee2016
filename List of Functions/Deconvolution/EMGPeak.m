function dataOut = EMGPeak(x, a)


% FROM: http://link.springer.com/article/10.1007%2FBF02268168#page-3

if a(3) > 0
    dataOut = exp(0.5*(a(2)/a(3))^2 - (x-a(1))/a(3)).*...
        (erf(1/sqrt(2)*(a(1)/a(2)+a(2)/a(3))) + ...
        erf(1/sqrt(2).*((x-a(1))/a(2)-a(2)/a(3))));
else
    dataOut = exp(0.5*(a(2)/a(3))^2 - (x-a(1))/a(3)).*...
        (erfc(1/sqrt(2)*(a(1)/a(2)+a(2)/a(3))) + ...
        erfc(1/sqrt(2).*((x-a(1))/a(2)-a(2)/a(3))));
end

if sum( ~isfinite( dataOut ) ) ~= 0 % the distribution is very close to a Gaussian, replace with a Gaussian
    dataOut = exp(-(x-a(1)).^2/(2*a(2)^2));

end

