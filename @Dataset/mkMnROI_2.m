%% DESCRIPTION
%
%% Copyright
% BSD 3-Clause License
% Copyright 2016-2017 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Data4ROI,X, Y] = mkMnROI_2(obj, mz_start, mz_end, tm_start, tm_end)

AxisX     = obj.AxisX.Data;
AxisY     = obj.AxisY.Data;

for ii = 1:size(mz_start, 1)
    ix = findCloser(mz_start(ii), AxisY);
    indMZ(ii, 1) = max(1, ix);
    ix = findCloser(mz_end(ii), AxisY);
    indMZ(ii, 2) = min(ix, length(AxisY));
    Data4ROI{ii} = [];
    Y{ii} = [];
    
    timesInt(ii, :) = [tm_start(ii) tm_end(ii)];
    X{ii} = AxisY(indMZ(ii, 1):indMZ(ii, 2));
    
end


TLim      = [min(timesInt(:,1)), max(timesInt(:,2))];
IdTLim1 = find(AxisX < TLim(1), 1, 'last');
if isempty(IdTLim1), IdTLim1 = 1; end

IdTLim2 = find(AxisX >  TLim(2), 1, 'first');
if isempty(IdTLim2), IdTLim2 = length(AxisX); end

IdTLim = [IdTLim1 IdTLim2];

switch obj.Format
    case 'profile'
        if ~isempty(obj.AxisY.Data)
            
            h = waitbar(0,'Making ROI, please wait');
            for ii = IdTLim(1):IdTLim(2)
                waitbar(ii/length(AxisX(:,1)))
                
                XMS = xpend(obj, obj.ListOfScans{ii});
                
                id = find(timesInt(:,1) <= AxisX(ii) & timesInt(:, 2) >= AxisX(ii));
                for jj = 1:length(id)
                    Data4ROI{id(jj)}(:, end+1) = XMS.Data(indMZ(id(jj),1):indMZ(id(jj),2), 2);
                    Y{id(jj)}(end+1) =  AxisX(ii);
                end
            end
            
            
        end
        
    case 'centroid'
        error('L37')
end
if ishandle(h), close(h); end

end
