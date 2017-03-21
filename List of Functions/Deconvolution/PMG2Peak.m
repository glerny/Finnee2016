function dataOut = PMG2Peak(axe, a)
% https://www.academia.edu/2160514/On_the_equations_describing_chromatographic_peaks_and_the_problem_of_the_deconvolution_of_overlapped_peaks

s = a(2) + a(3)*(axe-a(1)) + a(4)*(axe-a(1)).^2;
dataOut = a(2)./s.*exp(-(axe-a(1)).^2./(2*s.^2));

end

