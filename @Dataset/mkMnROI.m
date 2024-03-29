%% DESCRIPTION
%
%% Copyright
% BSD 3-Clause License
% Copyright 2016-2017 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tgtROIs = mkMnROI(obj, mzList, mzWdw, tmList, tmWdW, name, varList)

if isfile(fullfile(obj.Path2Fin, [name, '.mat']))
    tgtROIs = load(fullfile(obj.Path2Fin, name), 'tgtROIs');
    tgtROIs = tgtROIs.tgtROIs;
    return
end

AxisX     = obj.AxisX.Data;
AxisY     = obj.AxisY.Data;

if size(mzWdw, 1) == 1
    mzWdw = ones(size(mzList))*mzWdw;
end

if size(tmWdW, 1) == 1
    tmWdW = ones(size(mzList))*tmWdW;
end

varList = ones(size(mzList));

for ii = 1:length(mzList)
    ix = findCloser(mzList(ii), AxisY);
    indMZ(ii, 1) = max(1, ix-mzWdw(ii));
    indMZ(ii, 2) = min(ix+mzWdw(ii), length(AxisY));
    Data4ROI{ii} = [];
    
    timesInt(ii, 1) = max(AxisX(1), tmList(ii)-tmWdW(ii));
    timesInt(ii, 2) = min(tmList(ii)+tmWdW(ii), AxisX(end));
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
                end
            end
            
            
        end
        
    case 'centroid'
        error('L37')
end
if ishandle(h), close(h); end

h = waitbar(0,'Saving ROI, please wait');
for jj = 1:length(mzList)
    waitbar(jj/length(mzList))
    [InfoROI.Title, errmsg] = ...
        sprintf('ROI with tm from %.2f:%.2f %s and mz=%.4f +/- %u datapoints.', ...
        timesInt(jj, 1), timesInt(jj, 2), obj.AxisX.Unit, mzList(jj), mzWdw(jj));
    InfoROI.TagOfDts = obj.Log;
    InfoROI.TgtMz    = mzList(jj);
    InfoROI.MzWdw    = mzWdw(jj);
    InfoROI.TgtTm    = tmList(jj);
    InfoROI.tmWdW    = tmWdW(jj);
    InfoROI.TgtVarnc = varList(jj);
    
    AxisXc     = AxisX(timesInt(jj, 1) <= AxisX & timesInt(jj, 2) >= AxisX);
    AxisTm     = Axis(obj.AxisX.InfoAxis, AxisXc);
    InfoROI.AxisTm   = AxisTm;
    
    AxisYc  = AxisY(indMZ(jj,1):indMZ(jj,2));
    AxisMz = Axis(obj.AxisY.InfoAxis, AxisYc);
    InfoROI.AxisMZ   = AxisMz;
    InfoROI.Path2Fin = obj.Path2Fin;
    InfoROI.StoredData = Data4ROI{jj};
    
    s{jj} = ROI(InfoROI);
end
if ishandle(h), close(h); end

tgtROIs.mzList = mzList;
tgtROIs.mzWdw  = mzWdw;
tgtROIs.tmList = tmList;
tgtROIs.tmWdW  = tmWdW;
tgtROIs.ROI    = s;

if ~isempty(name)
    save(fullfile(obj.Path2Fin, name), 'tgtROIs')
end
end
