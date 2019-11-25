function PL = checkCo(PL)
Dtm     = 0.15;
Dmz     = 40e-06;
tM = GetDistance(PL.FOM{1}.Data(:,5), 1);
tM = tM + diag(inf(size(tM, 1), 1));
MZ = GetDistance(PL.FOM{1}.Data(:,10), 2);
MZ = MZ + diag(inf(size(MZ, 1), 1));
    
while 1
    [x, y] = find(tM <= Dtm & MZ <= Dmz);
    if isempty(x)
        break
    end
    lin = min(x(1), y(1));
    col = max(x(1), y(1));
    
    if isempty(x)
        break
    end
    cData = [PL.LstPIP{1}{lin}.Data; PL.LstPIP{1}{col}.Data];
    cData = sortrows(cData, 3);
    nPIP  = PL.LstPIP{1}{lin};
    nPIP.Data = cData;
    nPIP.x = PL.AxisX{1}.Data(min(cData(:,3)):max(cData(:,3)));
    nFOM = max([PL.FOM{1}.Data(lin, :); PL.FOM{1}.Data(col, :)]);
    nFOM(2:end-1) = nPIP.FOM;
    PL.FOM{1}.Data(lin, :) = nFOM;
    PL.LstPIP{1}{lin}      = nPIP;
    PL.FOM{1}.Data(col, :) = [];
    PL.LstPIP{1}(col)      = [];
    tM(col, :) = [];
    tM(:, col) = [];
    l2r        = sqrt((PL.FOM{1}.Data(lin, 5) -PL.FOM{1}.Data(:, 5)).^2);
    l2r(lin)   = inf;
    tM(lin, :) = l2r;
    tM(:, lin) = l2r;
    MZ(col, :) = [];
    MZ(:, col) = [];
    l2r = sqrt((PL.FOM{1}.Data(lin, 10) -PL.FOM{1}.Data(:, 10)).^2)/...
        PL.FOM{1}.Data(lin, 10);
    l2r(lin) = inf;
    MZ(lin, :) = l2r;
    MZ(:, lin) = l2r;
    
end
PL.FOM{1, 1}.Data(:,1) = 1:size(PL.FOM{1, 1}.Data(:,1), 1);