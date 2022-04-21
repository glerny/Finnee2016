function [S2N, maxPW, flag, BckgPts, XY] = getSignalToNoise_md1(XY)


LOQ = 10;
nPoly = 1;
minBckg = 0.25;
sg_order = 2;
sg_window = 5;
%%STEP 1: ALIGNED QC

%% Get the second derivative
dt = mean(XY(2:end, 1) - XY(1:end-1, 1));
[~,g] = sgolay(sg_order, sg_window);

XY(:,3) = conv(XY(:, 2), factorial(2)/(-dt)^2 * g(:,3), 'same'); 
flag = -1;
isBckg(:, 1) = 1:size(XY, 1);

%% step1. Use second derivatine
while 1
    
    Out = isoutlier(XY(isBckg, 3));
    if any(Out)
        isBckg(Out) = [];
        if size(isBckg, 1) <= minBckg*size(XY, 1)
            
            flag = 10;
            break
        else
            
            flag = 1;
        end
    else
        
        break
    end
end


%% step2. Use raw data
while 1
    p = polyfit(XY(isBckg, 1), XY(isBckg, 2), nPoly);
    y_ = polyval(p, XY(isBckg, 1));
    
    Out = isoutlier(XY(isBckg, 2) - y_);
    if any(Out)
        isBckg(Out) = [];
        if size(isBckg, 1) <= minBckg*size(XY, 1)
            
            flag = 20;
            break
        else
            
            flag = 2;
        end
    else
        
        break
    end
end

%% step3. calculate signal to noize

if flag == -1
    y_ = polyval(p, XY(:, 1));
    XY(:, 4) = XY(:, 2) - y_;
    BckgPts = false(size(XY, 1), 1);
    BckgPts(isBckg) = true;
    signal = max(XY(~BckgPts, 4));
    Noise  = 2*std(XY(BckgPts, 4));
    S2N    = NaN;
    maxPW  = NaN;
    
else
    y_ = polyval(p, XY(:, 1));
    XY(:, 4) = XY(:, 2) - y_;
    BckgPts = false(size(XY, 1), 1);
    BckgPts(isBckg) = true;
    signal = max(XY(~BckgPts, 4));
    Noise  = 2*std(XY(BckgPts, 4));
    S2N    = signal/Noise;
    maxPW  = max(diff(find(BckgPts)));
    
end


end

