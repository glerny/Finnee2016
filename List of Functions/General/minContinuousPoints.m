function mCP = minContinuousPoints(data, spacing, threshold)

narginchk(1, 3)
if nargin == 1
    spacing = 1;
    threshold = 0;
elseif nargin == 2
    threshold = 0;
end
    
   
Vector = data(:,1);
for ii = 1:spacing
    Vector(:, ii+1) = circshift(Vector(:,1) , ii);
end

Vector(Vector <= threshold) = 0;

mCP =  max(diff([0; find(sum(Vector, 2, 'omitnan') == 0); size(Vector, 1)+1]))-1-spacing;
end

