function Limits = getLimits_full(obj, methodId, mzWdW)

% 1- Initialisation and options
narginchk(2, 4);
if nargin == 2
    mzWdW = 5;
    tmWdW = 0.5;
end
if nargin == 3
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

ML = load(fullfile(obj.QC.Files{1}, 'myFinnee.mat'));
AxisMZ = ML.myFinnee.Datasets{1, 2}.AxisY.Data;
AxisTm = ML.myFinnee.Datasets{1, 2}.AxisX.Data; clear ML;

for ii = 1:size(FOM, 1)
    TgtInd(ii) = findCloser(FOM.AccurateMass(ii), AxisMZ);
end
TGI(:, 1) = unique(TgtInd);
TGI(1:end-1, 2) = diff(TGI(:, 1));

count = 1;
TgtInMz = [];

while 1
    if TGI(count, 2) == 1
        IdS = count;
        IdE = count;
        
        while TGI(count, 2) == 1
            IdE = IdE + 1;
            count = count+1;
        end
        TgtInMz(end+1, 1:2) = [ceil(mean(TGI(IdS:IdE,1))), IdE-IdS+1];
        count = count+1;
    else
        TgtInMz(end+1, 1:2) = [TGI(count,1), 1];
        count = count+1;
    end
    
    if count > size(TGI, 1)
        break
    end
    
    ii = size(TgtInMz, 1);
    Limits.ID(ii, 1) = FOM.ID(ii);
    IdS = max(1, TgtInMz(end, 1) - mzWdW - floor(TgtInMz(end, 2)/2));
    IdE = min(size(AxisMZ, 1), TgtInMz(end, 1) + mzWdW + floor(TgtInMz(end, 2)/2));
    Limits.mz_min(ii, 1) = AxisMZ(IdS);
    Limits.mz_max(ii, 1) = AxisMZ(IdE);
    Limits.Time_min(ii, 1) = AxisTm(1);
    Limits.Time_max(ii, 1) = AxisTm(end);
    
end
Limits = struct2table(Limits);