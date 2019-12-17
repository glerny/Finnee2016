function FOM = getFOM(obj, tgtTm, eTm, tgtMz, eMz, optInt)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


FOM = [];
XYZ = obj.AxisTm.Data;

XYZ(:,2) = trapz(obj.AxisMZ.Data, obj.StoredData);
XYZ(:,3) = trapz(obj.AxisMZ.Data, obj.StoredData.*obj.AxisMZ.Data);
XYZ(:,3) = XYZ(:,3)./XYZ(:,2);
XYZ(isnan(XYZ)) = 0;

if nnz(XYZ(:,2)) >= 0.9*length(XYZ(:,2))
    % do baseline correction
    [z, ~] = doPF(XYZ(:, 1:2), 2);
    XYZ(:,2) = XYZ(:,2) - z;
    XYZ(XYZ(:, 2) < 0, 2) = 0;
end

[~, I4LM] = LocalMaxima(XYZ(:, 1:2), 5, 0);

if any(I4LM)
    
    ix = abs(XYZ(I4LM, 3) - tgtMz) <= eMz & abs(XYZ(I4LM, 1) - tgtTm) <= eTm;
    if sum(ix) > 1
        [~, id] = max( XYZ(I4LM, 2).*ix);
        ix = false(size(ix));
        ix(id) = true;
    end
    
    if sum(ix) == 1
        
        if nargin == 5
            iLm = I4LM(ix);
            IdS = find(XYZ(1:iLm(1), 2) <= 0, 1, 'last');
            if isempty(IdS)
                IdS = 1;
            end
            
            IdE = find(XYZ(iLm(end):end, 2) <= 0, 1, 'first');
            if isempty(IdE)
                IdE = size(XYZ, 1);
            else
                IdE = IdE + iLm(end)-1;
            end
            
            if (IdE - IdS) >= 5
                
                FOM(end+1, :) = ChrMoment(XYZ(IdS:IdE, 1:2));
                XYZ(XYZ==0) = nan;
                FOM(end, 5) =  XYZ(I4LM(ix), 1);
                FOM(end, 6) =  XYZ(I4LM(ix), 2);
                FOM(end, 7) =  XYZ(I4LM(ix), 3);
                FOM(end, 8) = (IdE - IdS);
                FOM(end, 9) = obj.AxisTm.Data(IdS);
                FOM(end, 10) = obj.AxisTm.Data(IdE);
                FOM(end, 11) = sum(XYZ(IdS:IdE, 3).*XYZ(IdS:IdE, 2), 'omitnan')...
                    /sum(XYZ(IdS:IdE, 2), 'omitnan');
                FOM(end, 12) = mean(XYZ(IdS:IdE, 3), 'omitnan');
                FOM(end, 13) = std(XYZ(IdS:IdE, 3), 'omitnan');
            else
                FOM = NaN(1, 13);
                
            end
        else
            if any(isnan(optInt))
                FOM = NaN(1, 13);
            else
                IdS = find(XYZ(:, 1) <= optInt(1), 1, 'last');
                if isempty(IdS)
                    IdS = 1;
                end
                
                IdE = find(XYZ(:, 1) >= optInt(2), 1, 'first');
                if isempty(IdE)
                    IdE = size(XYZ, 1);
                end
                
                if sum(XYZ(IdS:IdE, 2) > 0) >= 5
                    M = ChrMoment(XYZ(IdS:IdE, 1:2));
                    if M(1) > 0
                        FOM(end+1, :) = M;
                        XYZ(XYZ==0) = nan;
                        FOM(end, 5) =  XYZ(I4LM(ix), 1);
                        FOM(end, 6) =  XYZ(I4LM(ix), 2);
                        FOM(end, 7) =  XYZ(I4LM(ix), 3);
                        FOM(end, 8) = (IdE - IdS);
                        FOM(end, 9) = obj.AxisTm.Data(IdS);
                        FOM(end, 10) = obj.AxisTm.Data(IdE);
                        FOM(end, 11) = sum(XYZ(IdS:IdE, 3).*XYZ(IdS:IdE, 2), 'omitnan')...
                            /sum(XYZ(IdS:IdE, 2), 'omitnan');
                        FOM(end, 12) = mean(XYZ(IdS:IdE, 3), 'omitnan');
                        FOM(end, 13) = std(XYZ(IdS:IdE, 3), 'omitnan');
                    else
                        FOM = NaN(1, 13);
                    end
                else
                    FOM = NaN(1, 13);
                    
                end
            end
        end
    else
        FOM = NaN(1, 13);
    end
else
    FOM = NaN(1, 13);
end