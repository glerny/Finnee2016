function [Aligned_Data, x] = doAlignment_MinPearson(Data, tmWdw)

step     = 1000;
threshCC = 0.75;
leeway   = inf;

%%STEP 1: ALIGNED QC
XY = Data(:, 1);
XY(:,2) = Data(:, 2);

Opt = zeros(size(Data, 2), 2);
% 1st Loop
for ii = 2:size(Data, 2)
    count = 1;
    TobeAligned = Data(:, ii);
    for x = -tmWdw:2*tmWdw/step:tmWdw
        [f(count, 2), model] = align(x);
        f(count, 1) = x;
        count = count+1;
    end
    
    [~, IdM] = max(f(:,1));
end

% 2nd Loop
Opt(:, 1) = Opt(:, 1) -  mean(Opt(2:end, 1), 'omitnan');
for ii = 1:length(Id4QC)
    TobeAligned = Data(:, [1 Id4QC(ii)]);
    [~, model(:, ii)] = align(Opt(ii, 2));
end

Mall = ChrMoment(XY);

Aligned_Data = Data(:,1);
for ii = 1:length(Id4QC)
    if nnz(Data(:, ii+1)) < 4
        x(ii) = NaN;
        model = Data(:, ii+1);
        
    else
        
        Mcur = ChrMoment([Data(:,1), Data(:, ii+1)]);
        x(ii) = fminsearch(@align, Mcur(2)-Mall(2));
        [~, model] = align(x(ii));
    end
    
    Aligned_Data(:, ii+1) = model;
end

% 2nd Loop
XY(:,2) =  mean(Aligned_Data(:, 2:end), 2);
Aligned_Data = Data(:,1);

for ii = 1:length(Id4QC)
    if nnz(Data(:, ii+1)) < 4
        x(ii) = NaN;
        model = Data(:, ii+1);
        
    else
        x(ii) = fminsearch(@align, x(ii));
        [~, model] = align(x(ii));
    end
    
    Aligned_Data(:, ii+1) = model;
end

%%STEP 2: ALIGNED All
XY(:,2) = mean(Aligned_Data(:, 2:end), 2);
Aligned_Data = Data(:,1);
Mall = ChrMoment(XY);

for ii = 1:size(Data, 2)-1
    if nnz(Data(:, ii+1)) < 4
        x(ii) = NaN;
        model = Data(:, ii+1);
        
    else
        
        Mcur = ChrMoment([Data(:,1), Data(:, ii+1)]);
        x(ii) = fminsearch(@align, Mcur(2)-Mall(2));
        [~, model] = align(x(ii));
    end
    
    Aligned_Data(:, ii+1) = model;
end

    function [f, model] = align(x)
        
        model = interp1(TobeAligned(:, 1)-x, TobeAligned(:, 2), XY(:, 1));
        if abs(x) > leeway, model = model*0; end
        model(isnan(model)) = 0;
        C = corrcoef(XY(:,2), model);
        f = C(1, 2);
        
        if ~isfinite(f)
            f = 0;
        end
        
    end

end

