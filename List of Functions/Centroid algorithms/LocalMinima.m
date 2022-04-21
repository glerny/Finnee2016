%% DESCRIPTION
% 1. INTRODUCTION
% LOCAL MAXIMA find all local maxima within a profile, where each local 
% maximum is any poinst that his higher or equal to its wdz closest 
% neighbourgs. The accurate position whith the X axis is determined by 
% fitting each local maximum and it two neighbourg points with a polynomial 
% of degree 2. The correspondinh intensity is also calculate as well as the 
% peak width at half height. 
%
% 2. PARAMETERS
%   .input: MS is a 2xn array that contain the profile spectrum with in the 
%   first column the m/z increments and in the m/z columns the frequency of 
%   ions detected.
%   .output: redMS id the reduced spectrum with in the first column the
%   accurate mass, in the second the corresponding intensity and in the
%   third the peak width at half height.
%
% 3. EXAMPLES
%       
% 4. COPYRIGHT
% Copyright 2014-2015 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal

function [redMS, I4LM] = LocalMaxima(XY, wdz, thrI)

%% CORE OF THE FUNCTION
% Finding non nul local maxima
yy = XY(:,2);
A = zeros(length(yy), 2*wdz + 1);
for ii = 1:2*wdz+1
    ind1 = max(1, wdz+2-ii);
    ind2 = min(length(yy), length(yy) + wdz+1-ii);
    A(length(yy)-ind2+1:length(yy)-ind1+1,ii) = yy(ind1:ind2);
end

lmax = yy <= min(A, [], 2) & yy <= thrI;
% calculate accurate masses with the three most intense
% points
jj   = 1;
a    = [];
b    = [];
c    = [];
I4LM = []; % to record the indices of the peak maxima.
lmax(1) = false;
lmax(end) = false;

while 1
    if lmax(jj)
        % We used the analytical solution for the polynomial of degree 2 as
        % with 3 points there is only one solution. This is much faster than
        % a fitting procedure.
        a(end+1) = ((XY(jj-1,2)-XY(jj,2))/(XY(jj-1,1)-XY(jj,1))...
            -(XY(jj-1,2)-XY(jj+1,2))/(XY(jj-1,1)-XY(jj+1,1)))...
            /(XY(jj,1)-XY(jj+1,1));
        b(end+1) = ((XY(jj-1,2)-XY(jj+1,2))-a(end)...
            *(XY(jj-1,1)^2-XY(jj+1,1)^2))...
            /(XY(jj-1,1)-XY(jj+1,1));
        c(end+1) = XY(jj, 2) - a(end)*XY(jj, 1)^2 - b(end)*XY(jj, 1);
        I4LM(end+1) = jj;
        if lmax(jj+1)   % Check if the following data points is also a 
                        % local maxima. If yes, skip it.
            jj = jj + 1;
            lmax(jj+1) = false;
        end
    end
    jj = jj + 1;
    if jj >= length(lmax), break; end
end

if ~isempty(a)
    redMS      = (-b./(2*a))';
    redMS(:,2) =(a'.*redMS(:,1).^2) + (b'.*redMS(:,1)) + c';
    
    % redMS is nan when the peak contain only one point
    redMS(isnan(redMS(:,1)) | isnan(redMS(:,2)), :) = [];
else
    redMS = [];
end