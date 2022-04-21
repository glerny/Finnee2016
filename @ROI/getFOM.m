function FOM = getFOM(obj, nGol, wdwGol, wdwLM, Lim)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

minPts = 4;

%% 1. Initialisation
FOM.M0  = [];
FOM.M1  = [];
FOM.M2  = [];
FOM.M3  = [];
FOM.TimeAtPeakMax = [];
FOM.IntAtPeakMax = [];
FOM.nbrPts = [];
FOM.Lm1 = [];
FOM.Lm2 = [];
FOM.mean_mz = [];
FOM.std_mz = [];
FOM.accuMass = [];

if isempty(obj.StoredData), FOM = struct2table(FOM); return; end

%% 2. Getting trace and filtering
XYZ = obj.AxisTm.Data;
XYZ(:,2) = trapz(obj.AxisMZ.Data, obj.StoredData);
XYZ(:,3) = trapz(obj.AxisMZ.Data, obj.StoredData.*obj.AxisMZ.Data);
XYZ(:,3) = XYZ(:,3)./XYZ(:,2);
I2rem = XYZ(:,1) < obj.TgtTm - obj.tmWdW | XYZ(:,1) > obj.TgtTm + obj.tmWdW;
XYZ(I2rem, :) = [];
XYZ(isnan(XYZ)) = 0;

%% 3. Measuring Figures of Merits
if isempty(Lim)
    XYZ(:, 4) = sgolayfilt(XYZ(:,2), nGol, wdwGol);
    [~, I4LM] = LocalMaxima(XYZ(:, [1 4]), wdwLM, 0);
    
    while ~isempty(I4LM)
        ix1 = find(XYZ(1:I4LM(1), 4) <= 0, 1, 'last');
        ix2 = find(XYZ(I4LM(1):end, 4) <= 0 , 1, 'first');
        if isempty(ix1), ix1 = 1; end
        if isempty(ix2), ix2 = size(XYZ, 1); else, ix2 = ix2 + I4LM(1) -1; end
        
        if(ix2-ix1) < minPts
            [~,ia] = intersect(I4LM, ix1:ix2);
            I4LM(ia) = [];
            [~, ia] = intersect(I4LM, ix1:ix2);
            I4LM(ia) = [];
            continue
        end
        
        CM = ChrMoment(XYZ(ix1:ix2, 1:2));
        FOM.M0(end+1, 1) = CM(1);
        FOM.M1(end+1, 1) = CM(2);
        FOM.M2(end+1, 1) = CM(3);
        FOM.M3(end+1, 1) = CM(4);
        [~, iMax] = max(XYZ(ix1:ix2, 2)); iMax = iMax + ix1 -1;
        FOM.TimeAtPeakMax(end+1, 1) = XYZ(iMax, 1);
        FOM.IntAtPeakMax(end+1, 1) = XYZ(iMax, 2);
        FOM.nbrPts(end+1, 1) = ix2-ix1-1;
        FOM.Lm1(end+1, 1) = XYZ(ix1, 1);
        FOM.Lm2(end+1, 1) = XYZ(ix2, 1);
        FOM.mean_mz(end+1, 1) = mean(nonzeros(XYZ(ix1:ix2, 3)));
        FOM.std_mz(end+1, 1) = std(nonzeros(XYZ(ix1:ix2, 3)));
        FOM.accuMass(end+1, 1) = sum(XYZ(ix1:ix2, 2).*XYZ(ix1:ix2, 3))/sum(XYZ(ix1:ix2, 2));
        [~,ia] = intersect(I4LM, ix1:ix2);
        I4LM(ia) = [];
        [~, ia] = intersect(I4LM, ix1:ix2);
        I4LM(ia) = [];
    end
else
    ix1 = find(XYZ(:, 1) <= Lim(1), 1, 'last');
    ix2 = find(XYZ(:, 1) >= Lim(2), 1, 'first');
    if isempty(ix1), ix1 = 1; end
    if isempty(ix2), ix2 = size(XYZ, 1); end
    
    if nnz(XYZ(ix1:ix2, 2)) < minPts
        FOM = struct2table(FOM); return;
    end
    
    CM = ChrMoment(XYZ(ix1:ix2, 1:2));
    FOM.M0(end+1, 1) = CM(1);
    FOM.M1(end+1, 1) = CM(2);
    FOM.M2(end+1, 1) = CM(3);
    FOM.M3(end+1, 1) = CM(4);
    [~, iMax] = max(XYZ(ix1:ix2, 2)); iMax = iMax + ix1 -1;
    FOM.TimeAtPeakMax(end+1, 1) = XYZ(iMax, 1);
    FOM.IntAtPeakMax(end+1, 1) = XYZ(iMax, 2);
    FOM.nbrPts(end+1, 1) = ix2-ix1-1;
    FOM.Lm1(end+1, 1) = XYZ(ix1, 1);
    FOM.Lm2(end+1, 1) = XYZ(ix2, 1);
    FOM.mean_mz(end+1, 1) = mean(nonzeros(XYZ(ix1:ix2, 3)));
    FOM.std_mz(end+1, 1) = std(nonzeros(XYZ(ix1:ix2, 3)));
    FOM.accuMass(end+1, 1) = sum(XYZ(ix1:ix2, 2).*XYZ(ix1:ix2, 3))/sum(XYZ(ix1:ix2, 2));
end

FOM = struct2table(FOM);