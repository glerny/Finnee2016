function out = centre3D(X, Y, Area)

if sum(sum(Area == 0)) > 0.5*numel(Area)
    out = [NaN, NaN, NaN];
else
    tot_mass = sum(Area(:));
    [ii,jj] = ndgrid(1:size(Area,1),1:size(Area,2));
    R = sum(ii(:).*Area(:))/tot_mass;
    C = sum(jj(:).*Area(:))/tot_mass;
    
    [TX, TY] = meshgrid(X, Y);
    out(1) = tot_mass/numel(Area);
    p = polyfit(1:size(X), X', 1);
    out(2) = polyval(p, R);
    p = polyfit(1:size(Y), Y', 1);
    out(3) = polyval(p, C);
end
end

