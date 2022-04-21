function M = ChrMoment3D(X, Y, XY)

Y2X = (trapz(X, XY'))';
X2Y = (trapz(Y, XY))';

M2X = ChrMoment([X, X2Y]);
M2Y = ChrMoment([Y, Y2X]);

M = array2table([M2X(1:4), M2Y(1:4)]);
M.Properties.VariableNames = {'M0_X', 'M1_X', 'M2_X', 'M3_X', ...
    'M0_Y', 'M1_Y', 'M2_Y', 'M3_Y'};