function FOM = getFOM(obj, tgtTm, eTm, tgtMz, eMz, optInt)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

IThres = 0.05;
minRes = 1.5;

[out, lmax, Nz, obj] = obj.findPeaks(5, 15);
if ~isempty(out)
    ix = abs(out(:,1) - tgtMz) <= eMz & abs(out(:,2) - tgtTm) <= eTm;
    FOM = [];
    
    if any(ix)
        
        if nargin == 5
            XY = sum(obj.Smoothed);
            iLm = lmax(ix);
            IdS = find(XY(1:iLm(1)) <= Nz, 1, 'last');
            if isempty(IdS)
                IdS = 1;
            end
            
            IdE = find(XY(iLm(end):end) <= Nz, 1, 'first');
            if isempty(IdE)
                IdE = size(XY, 2);
            else
                IdE = IdE + iLm(end)-1;
            end
            
            LF = find(lmax >= IdS &  lmax <= IdE & out(:,3)' >= IThres*max(out(ix,3)));
            
            if isempty(LF)
                disp('no can do')
            elseif (IdE - IdS) >= 10
                x  = obj.AxisTm.Data(IdS:IdE);
                y  = obj.AxisMZ.Data;
                xy = obj.StoredData(:, IdS:IdE);
                
                FOM(end+1, :) = ChrMoment3D(x, y, xy);
                FOM(end, 6) = length(LF);
                FOM(end, 7) = obj.AxisTm.Data(IdS);
                FOM(end, 8) = obj.AxisTm.Data(IdE);
            else
                FOM = NaN(1, 8);
                
            end
        elseif nargin == 6
            
            X = obj.AxisTm.Data;
            IdS = find(X < optInt(1), 1, 'last');
            if isempty(IdS)
                IdS = 1;
            end
            
            IdE = find(X > optInt(2), 1, 'first');
            if isempty(IdE)
                IdE = size(X, 2);
            end
            
            LF = find(lmax >= IdS &  lmax <= IdE & out(:,3)' >= IThres*max(out(ix,3)));
            
            if isempty(LF)
                FOM = NaN(1, 8);
            elseif (IdE - IdS) >= 10
                x  = obj.AxisTm.Data(IdS:IdE);
                y  = obj.AxisMZ.Data;
                xy = obj.StoredData(:, IdS:IdE);
                
                FOM(end+1, :) = ChrMoment3D(x, y, xy);
                FOM(end, 6) = length(LF);
                FOM(end, 7) = obj.AxisTm.Data(IdS);
                FOM(end, 8) = obj.AxisTm.Data(IdE);
            else
                FOM = NaN(1, 8);
                
            end
        end
    else
        FOM = NaN(1, 8);
    end
else
    FOM = NaN(1, 8);
end






