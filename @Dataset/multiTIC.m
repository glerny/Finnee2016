%% DESCRIPTION
%
%% Copyright
% BSD 3-Clause License
% Copyright 2016-2017 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function TIC = multiTIC(obj, Limits, k)
%TODO: Check if the dataset is mz normalised

AxisX     = obj.AxisX.Data;
AxisY     = obj.AxisY.Data;

for ii = 1:size(Limits, 1)
   
    ix = findCloser(Limits.mz_min(ii), AxisY);
    indMZ(ii, 1) = max(1, ix);
    ix = findCloser(Limits.mz_max(ii), AxisY);
    indMZ(ii, 2) = min(ix, length(AxisY));
    timesInt(ii, :) = [Limits.Time_min(ii) Limits.Time_max(ii)];
    
    TIC{ii, 1}(:, 1) = AxisX(AxisX >= timesInt(ii, 1) & AxisX <= timesInt(ii, 2));
    TIC{ii, 1}(:, 2) = 0;
    
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
            
            h = waitbar(0,'Bulding Total Ion Profiles, please wait');
            for ii = IdTLim(1):IdTLim(2)
                waitbar(ii/length(AxisX(:,1)))
                
                XMS = xpend(obj, obj.ListOfScans{ii});
                
                id = find(timesInt(:,1) <= AxisX(ii) & timesInt(:, 2) >= AxisX(ii));
                for jj = 1:length(id)
                    if nnz(XMS.Data(indMZ(id(jj),1):indMZ(id(jj),2), 2)) >= k
                        TIC{id(jj)}(TIC{id(jj)}(:, 1) == AxisX(ii), 2) = ...
                            sum(XMS.Data(indMZ(id(jj),1):indMZ(id(jj),2), 2));
                    end
                end
            end
            
            
        end
        
    case 'centroid'
        error('L37')
end
if ishandle(h), close(h); end

end
