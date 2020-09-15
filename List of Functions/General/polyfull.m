function [p, IOutl, XY] = polyfull(x, y, n, w)
    XY = (1:length(x))';
    XY(:, 2) = x;
    XY(:, 3) = y;
    if isempty(w)
        XY(:, 4) = 1;
    else
        XY(:, 4) = w;
    end
    XY(isnan(XY(:,2)) | isnan(XY(:,3)), :) = [];
    p  = polyfitweighted(XY(:,2), XY(:,3), n,  1./XY(:,4));
    IOutl = isoutlier( XY(:,3) - polyval(p, XY(:,2)));
    XY(IOutl, :) = [];
    p  = polyfitweighted(XY(:,2), XY(:,3), n,  1./XY(:,4));
end