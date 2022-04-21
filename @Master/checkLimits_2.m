function [newLimits, TISs, MzAxis, Options] = checkLimits_2(obj, dts, Limits, Options)

% TODO: Pass variable as options

if nargin == 3
    Options.minMzPts   = 3;
    Options.nbrQc      = 10;
    Options.Thres      = 0;
    Options.minPPP     = 5;
end

myCount   = 0;
tic;

QCFiles = obj.QC.Files;
% 2. Setting Limits for ROIs
% 2.1. Open first QC
MF = load(fullfile(QCFiles{1}, 'myFinnee.mat'));
myFinnee = MF.myFinnee; clear MF;

fprintf('\n\t Initialisation with %s\n', myFinnee.FileID)
MzAxis = myFinnee.Datasets{dts}.AxisY;


fprintf('\t Calculating the TIS\n')
% TODO: Change .MkMnROI_2 to .MkMnROI_det for improved speed.
TIS = myFinnee.Datasets{dts}.multiTIS(Limits, Options.minMzPts);

count = 1;

% 2.3. Add remaining QCs
for cF = 2:length(QCFiles)
    count = count+1;
    MF = load(fullfile(QCFiles{cF}, 'myFinnee.mat'));
    myFinnee = MF.myFinnee; clear MF;
    fprintf('\n\t Adding %s\n', myFinnee.FileID)
    fprintf('\t Calculating the TIP\n')
    newTIS = myFinnee.Datasets{dts}.multiTIS(Limits, Options.minMzPts);
    
    fprintf('\t Interpreting to common time axis\n')
    for ii = 1:size(newTIS, 1)
        vq = interp1(newTIS{ii}(:, 1), newTIS{ii}(:, 2), TIS{ii}(:, 1));
        vq(isnan(vq)) = 0;
        TIS{ii}(:, count+1) = vq;
    end
end

if Options.nbrQc > count, Options.nbrQc = count; end

fprintf('\n\n\t Splitting the TIP\n')

for ii = 1:size(Limits, 1)
    XY = TIS{ii};
    filt_XY = sum(sum(XY(:, 2:end) > Options.Thres, 2) >= Options.nbrQc, 2);
    filt_XY(:,2) = mean(XY(:, 2:end), 2);
    filt_XY(:,3) = filt_XY(:,1).*filt_XY(:,2);
    
    while 1
        [myMax, myId] = max(filt_XY(:,3));
        if myMax == 0
            break
            
        else
            Id_start = find(filt_XY(1:myId, 2)<= Options.Thres, 1, "last");
            if isempty(Id_start), Id_start = 1; end
            
            Id_end = find(filt_XY(myId:end, 2)<= Options.Thres, 1, "first");
            if isempty(Id_end)
                Id_end = size(filt_XY, 1);
            else
                Id_end = myId + Id_end - 1;
            end
            
            if Id_end - Id_start >= Options.minPPP
                myCount                        = myCount + 1;
                newLimits.ID(myCount, 1)       = nan;
                
                if filt_XY(Id_start, 2) <= Options.Thres
                    newLimits.mz_min(myCount, 1) = XY(Id_start, 1);
                else
                    Id = find(MzAxis.Data == XY(Id_start, 1));
                    Id = max(1, Id-1);
                    newLimits.mz_min(myCount, 1) = ...
                        MzAxis.Data(Id);
                end
                
                if filt_XY(Id_end, 2) <= Options.Thres
                    newLimits.mz_max(myCount, 1) = XY(Id_end, 1);
                else
                    
                    Id = find(MzAxis.Data == XY(Id_end, 1));
                    Id = min(size(MzAxis.Data, 1), Id+1);
                    newLimits.mz_max(myCount, 1) = ...
                        MzAxis.Data(Id);
                    
                end
                
                
                TISs{myCount} = XY(Id_start:Id_end, :);
                newLimits.Time_min(myCount, 1) = Limits.Time_min(ii);
                newLimits.Time_max(myCount, 1) = Limits.Time_max(ii);
                newLimits.ID_TIS(myCount, 1)  = ii;
            end
            
            filt_XY(Id_start:Id_end, 1) = 0;
            filt_XY(:,3) = filt_XY(:,1).*filt_XY(:,2);
        end
    end
end
newLimits = struct2table(newLimits);
newLimits.ID = (1:size(newLimits, 1))'; 

Options.TOC     = toc;