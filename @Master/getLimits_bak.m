function Limits = getLimits(obj, methodId, FilterOut, mzWdW, tmWdW)

% 1- Initialisation and options
narginchk(2, 5);
if nargin == 3
    mzWdW = 5;
    tmWdW = 0.5;
end
if nargin == 4
    tmWdW = 0.5;
end

switch methodId
    case 1
        FOM = obj.QC.Method1.FOM;
        
    case 2
        error('not implemented yet')
        
    case 3
        error('not implemented yet')
        
end

% 2- Step 1. Find Limits
FOM(FilterOut, :) = [];
ML = load(fullfile(obj.QC.Files{1}, 'myFinnee.mat'));
AxisMZ = ML.myFinnee.Datasets{1, 2}.AxisY.Data;
AxisTm = ML.myFinnee.Datasets{1, 2}.AxisX.Data; clear ML;
Lim = obj.QC.Method1.Limits;

for ii = 1:size(FOM, 1)
    Limits.ID(ii, 1) = FOM.ID(ii);
    TgtInd = findCloser(FOM.AccurateMass(ii), AxisMZ);
    if isempty(TgtInd)
        disp('pp')
    end
    
    IdS = max(1, TgtInd - mzWdW);
    IdE = min(size(AxisMZ, 1), TgtInd + mzWdW);
    Limits.mz_min(ii, 1) = AxisMZ(IdS);
    Limits.mz_max(ii, 1) = AxisMZ(IdE);
    Limits.Time_min(ii, 1) = max(AxisTm(1)...
        , FOM.Time_Start(ii) - tmWdW);
    Limits.Time_max(ii, 1) = min(AxisTm(end)...
        , FOM.Time_End(ii) + tmWdW);
end
Limits = struct2table(Limits);

Limits = sortrows(Limits,'Time_max','ascend');
Limits.ID = (1:size(Limits, 1))';

while 1
    LoopMe = true;
    myCount = 1;
    
    while 1
        id = Limits.Time_min <= Limits.Time_max(myCount) & (...
            (Limits.mz_max >= Limits.mz_min(myCount) & Limits.mz_max <= Limits.mz_max(myCount)) |...
            (Limits.mz_min >= Limits.mz_min(myCount) & Limits.mz_min <= Limits.mz_max(myCount)) |...
            (Limits.mz_min <= Limits.mz_min(myCount) & Limits.mz_max >= Limits.mz_max(myCount)));
        
        if sum(id) > 1
            Tbl = Limits(id, :);
            Limits(id, :) = [];
            NewLine.ID = nan;
            NewLine.mz_min = min(Tbl.mz_min);
            NewLine.mz_max = max(Tbl.mz_max);
            NewLine.Time_min = min(Tbl.Time_min);
            NewLine.Time_max = max(Tbl.Time_max);
            Limits = [Limits; struct2table(NewLine)];
            Limits = sortrows(Limits,'Time_min','ascend');
            Limits.ID = (1:size(Limits, 1))';
            LoopMe = false;
            if myCount > size(Limits, 1)
                break
            end
            
        else
            myCount = myCount+1;
            if myCount > size(Limits, 1)
                break
            end
        end
    end
    if LoopMe, break; end
end
