function BCKG = MADPeaks(XY)

sumOut = 0;
outl = false(size(XY, 1), 1);
S0 = sum(outl);
BCKG = XY;

while 1
    io = isoutlier(BCKG(:,2) - mean(BCKG(:,2)));
    if sum(io) == 0
        break
    end
    BCKG(io, :) = [];

end


