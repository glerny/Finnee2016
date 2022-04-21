function [Aligned_Data, Opt, leeway] = doAlignment_MinPearson(Data, IdQC, leeway, robust, ThresholdFactor)

if nargin == 2
    leeway = 1;
    robust = false;
elseif nargin == 3
    robust = false;
end

if nargin == 4 & robust
    ThresholdFactor = 5;
end
minLeeway = leeway/100;
Edge = 5;

%%STEP 1: ALIGNED QC
IdQC = find(IdQC);
XY = Data(:, 1);
XY(:,2) = Data(:, IdQC(1)+1);
Opt = zeros(size(IdQC, 1), 2);
% 1st Loop

[~, IdM] = max(XY(:, 2));
f = (-leeway:leeway/100:leeway)';

for ii = 2:size(IdQC, 1)
    
    TobeAligned = Data(:, [1 IdQC(ii)+1]);
    for x = 1:size(f, 1)
        [f(x, 2), model] = align(f(x, 1));
    end
    
    f(:,2) = 1-f(:,2);
    [~, IdM] = max(f(:,2), [], 'omitnan');
    [Opt(ii, 1), Opt(ii, 2)] = fminsearch(@align,f(IdM,1));
end

if robust
    io = ~isoutlier(Opt(:,1), 'ThresholdFactor', ThresholdFactor);
else
    io = true(size(IdQC));
end

Opt(io, 1) = Opt(io, 1) -  mean(Opt(io, 1), 'omitnan');
clear Aligned_Data
for ii = 1:size(io, 1)
    TobeAligned = Data(:, [1 IdQC(ii)+1]);
    [~, Aligned_Data(:, ii)] = align(Opt(ii, 1));
end

%%STEP 2: ALIGNED All
XY = Data(:, 1);
XY(:,2) = mean(Aligned_Data, 2, 'omitnan');
Aligned_Data = XY(:, 1);
leeway = max(5*std(Opt(io, 1)), minLeeway);

Opt = zeros(size(Data, 2)-1, 2);
[~, IdM] = max(XY(:, 2));
f = (-leeway:leeway/100:leeway)';
for ii = 2:size(Data, 2)
    TobeAligned = Data(:, [1 ii]);
    for x = 1:size(f, 1)
        f(x, 2) = align(f(x, 1));
    end
    
    f(:,2) = 1-f(:,2);
    [~, IdM] = max(f(:,2), [], 'omitnan');
    [Opt(ii-1, 1), Opt(ii-1, 2)] = fminsearch(@align,f(IdM,1));
    Aligned_Data(:, ii) = ...
        interp1(TobeAligned(:, 1)-Opt(ii-1, 1), TobeAligned(:, 2), XY(:, 1));
end
Opt(:, 2) = 1-Opt(:, 2);


    function [f, model] = align(x)
        
        model = interp1(TobeAligned(:, 1)-x, TobeAligned(:, 2), XY(:, 1));
        if abs(x) > leeway, model = model*0; end
        model(isnan(model)) = 0;
        C = corrcoef(XY(:,2), model);
        f = 1 - abs(C(1, 2));
        
        if ~isfinite(f)
            f = 1;
        end
        
    end

end

