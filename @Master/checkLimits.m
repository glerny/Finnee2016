function [newLimits, TIPs, TimeAxis, Options, iscomplete] = checkLimits(obj, dts, Limits, Options)

% TODO: Pass variable as options

if nargin == 3
    Options.minMzPts   = 3;
    Options.WdwLength  = 25;
    Options.minPPP     = 25;
    Options.nbrQc      = 10;
    Options.LoopMax    = 2;
    Options.Thres      = 10;
    Options.maxPeakLng = 2; %min
    Options.stepLng    = 0.25;
    Options.cut4Bsl    = 1;
    Options.spkSz      = 3;
    Options.AddLength  = 0.3;
end
ML = floor(Options.maxPeakLng/Options.stepLng);
iscomplete = true;

myCount   = 0;
tic;

QCFiles = obj.QC.Files;
% 2. Setting Limits for ROIs
% 2.1. Open first QC
MF = load(fullfile(QCFiles{1}, 'myFinnee.mat'));
myFinnee = MF.myFinnee; clear MF;

fprintf('\n\t Initialisation with %s\n', myFinnee.FileID)
TimeAxis = myFinnee.Datasets{dts}.AxisX;


fprintf('\t Calculating the TIP\n')
% TODO: Change .MkMnROI_2 to .MkMnROI_det for improved speed.
TIP = myFinnee.Datasets{dts}.multiTIC(Limits, Options.minMzPts);

count = 1;

% 2.3. Add remaining QCs
for cF = 2:length(QCFiles)
    count = count+1;
    MF = load(fullfile(QCFiles{cF}, 'myFinnee.mat'));
    myFinnee = MF.myFinnee; clear MF;
    fprintf('\n\t Adding %s\n', myFinnee.FileID)
    fprintf('\t Calculating the TIP\n')
    newTIP = myFinnee.Datasets{dts}.multiTIC(Limits, Options.minMzPts);
    
    fprintf('\t Interpreting to common time axis\n')
    for ii = 1:size(newTIP, 1)
        vq = interp1(newTIP{ii}(:, 1), newTIP{ii}(:, 2), TIP{ii}(:, 1));
        vq(isnan(vq)) = 0;
        TIP{ii}(:, count+1) = vq;
    end
end

if Options.nbrQc > count, Options.nbrQc = count; end

fprintf('\n\n\t Splitting the TIP\n')

for ii = 1:size(Limits, 1)
    cX = TIP{ii}(:,1);
    XY = TIP{ii}(:,2:end);
    for jj = 1:size(XY, 2)
        XY(:, jj) = spikesRemoval_2(XY(:, jj), Options.spkSz);
    end
    
    XY = smoothdata(XY, 'movmean', Options.WdwLength);
    XY(XY < Options.Thres) = 0;
    xy = sum(XY ~= 0, 2);
    curLoop = 1;
    
    while 1
        Id = find(xy >= Options.nbrQc, 1, "first");
        if isempty(Id)
            if curLoop > Options.LoopMax, break; end
            
            XY = smoothdata(XY, 'movmean', Options.WdwLength);
            XY(XY < Options.Thres) = 0;
            xy = sum(XY ~= 0, 2);
            curLoop = curLoop + 1;
            
        else
            Id_start = find(xy(1:Id) <= Options.cut4Bsl, 1, "last");
            if isempty(Id_start), Id_start = 1; end
            
            Id_end = find(xy(Id:end) <= Options.cut4Bsl, 1, "first");
            if isempty(Id_end)
                Id_end = size(xy, 1);
            else
                Id_end = Id + Id_end - 1;
            end
            
            if Id_end - Id_start >= Options.minPPP + curLoop*Options.WdwLength
                myCount                        = myCount + 1;
                newLimits.ID(myCount, 1)       = nan;
                newLimits.mz_min(myCount, 1)   = Limits.mz_min(ii);
                newLimits.mz_max(myCount, 1)   = Limits.mz_max(ii);
                if Id_start == 1 && mean(XY(Id_start, :)) >= Options.Thres
                    newLimits.Time_min(myCount, 1) = ...
                        max(TimeAxis.Data(1), cX(1) - Options.AddLength);
                else
                    newLimits.Time_min(myCount, 1) = cX(Id_start);
                end
                
                
                if Id_end == size(xy, 1) && mean(XY(Id_end, :)) >= Options.Thres
                    newLimits.Time_max(myCount, 1) = ...
                        min(TimeAxis.Data(end), cX(end) + Options.AddLength);
                else
                    newLimits.Time_max(myCount, 1) = cX(Id_end);
                end
                
                newLimits.ID_TIPs(myCount, 1)  = ii;
                
                Ys = TIP{ii}(Id_start:Id_end, 2:end); 
                X  = TIP{ii}(Id_start:Id_end, 1);
                Moment = [];
                
                for jj = 1:size(Ys, 2)
                    Moment = [Moment; ChrMoment([X, Ys(:, jj)])];
                end
                newLimits.Moment{myCount, 1}  = Moment;
                TIPs{myCount, 1} = TIP{ii}(Id_start:Id_end, :);
            end
            
            XY(Id_start:Id_end, :) = 0;
            xy = sum(XY ~= 0, 2);
            
        end
    end
end
newLimits = struct2table(newLimits);
newLimits.ID = (1:size(newLimits, 1))'; 

Dm = newLimits.Time_max - newLimits.Time_min;
ID2Cut = find(Dm > Options.maxPeakLng);
myCount = 0;

for ii = 1:length(ID2Cut)
    ID2TIPs = ID2Cut(ii);
    
    XY = TIPs{ID2TIPs, 1}(:, 2:end);
    Tm = TIPs{ID2TIPs, 1}(:, 1);
    
    for jj = 1:size(XY, 2)
        
        xy = [Tm, XY(:, jj)];
        
        pst = Tm(1) -(ML-1)*Options.stepLng;
        I1 = findCloser(Tm,  pst);
        I2 = findCloser(Tm,  pst+Options.maxPeakLng);
        while 1
            [~, w] = doPF(xy(I1:I2, 1:2), 2);
            xy(I1:I2, end+1) = w;
            
            pst = pst + Options.stepLng;
            if  pst > Tm(end) - Options.stepLng
                break
            end
            I1 = findCloser(Tm,  pst);
            I2 = findCloser(Tm,  pst+Options.maxPeakLng);
            if I2 == 1; I2 = length(Tm); end 
        end
        xy(:, 3) = sum(xy(:, 3:end), 2);
    end
    XY(xy(:, 3) > ML/2, :) = 0;
    
    
    for jj = 1:size(XY, 2)
        XY(:, jj) = spikesRemoval_2(XY(:, jj), Options.spkSz);
    end
    XY = smoothdata(XY, 'movmean', Options.WdwLength);
    XY(XY < Options.Thres) = 0;
    xy = sum(XY ~= 0, 2);
    curLoop = 1;
    
    while 1
        Id = find(xy >= Options.nbrQc, 1, "first");
        if isempty(Id)
            if curLoop > Options.LoopMax, break; end
            
            XY = smoothdata(XY, 'movmean', Options.WdwLength);
            XY(XY < Options.Thres) = 0;
            xy = sum(XY ~= 0, 2);
            curLoop = curLoop + 1;
            
        else
            Id_start = find(xy(1:Id) <= Options.cut4Bsl, 1, "last");
            if isempty(Id_start), Id_start = 1; end
            
            Id_end = find(xy(Id:end) <= Options.cut4Bsl, 1, "first");
            if isempty(Id_end)
                Id_end = size(xy, 1);
            else
                Id_end = Id + Id_end - 1;
            end
            
            if Id_end - Id_start >= Options.minPPP + curLoop*Options.WdwLength
                myCount                        = myCount + 1;
                addLimits.ID(myCount, 1)       = nan;
                addLimits.mz_min(myCount, 1)   = newLimits.mz_min(ID2Cut(ii));
                addLimits.mz_max(myCount, 1)   = newLimits.mz_max(ID2Cut(ii));
                addLimits.Time_min(myCount, 1) = Tm(Id_start);
                addLimits.Time_max(myCount, 1) = Tm(Id_end);
                addLimits.ID_TIPs(myCount, 1)  = ID2TIPs;
                
                Ys = TIPs{ID2TIPs}(Id_start:Id_end, 2:end); 
                X  = TIPs{ID2TIPs}(Id_start:Id_end, 1);
                Moment = [];
                
                for jj = 1:size(Ys, 2)
                    Moment = [Moment; ChrMoment([X, Ys(:, jj)])];
                end
                addLimits.Moment{myCount, 1}  = Moment;
                addTIPs{myCount, 1} = TIPs{ID2TIPs, 1}(Id_start:Id_end, :);
            end
            
            XY(Id_start:Id_end, :) = 0;
            xy = sum(XY ~= 0, 2);
        end
    end
    
end

try
    if ~isempty(ID2Cut)
        newLimits(ID2Cut, :) = [];
        newLimits = [newLimits; struct2table(addLimits)];
        newLimits.ID = (1:size(newLimits, 1))';
        TIPs(ID2Cut, :) = [];
        TIPs = [TIPs; addTIPs];
    end
catch
    disp("pp")
end

Options.TOC     = toc;