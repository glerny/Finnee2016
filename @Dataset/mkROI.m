%% DESCRIPTION
%
%% Copyright 
% BSD 3-Clause License
% Copyright 2016-2017 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function s = mkROI(obj, tgtMz, mzWdw, timeInt)
if length(timeInt) == 1
    [timeStart, timeEnd] = deal(timeInt(1));
else
    timeStart = min(timeInt);
    timeEnd   = max(timeInt);
end
AxisX      = obj.AxisX.Data;
AxisTm     = obj.AxisX;
indTimeStt = findCloser(timeStart, AxisX);
indTimeEnd = findCloser(timeEnd, AxisX);
AxisX      = AxisX(indTimeStt:indTimeEnd);
AxisTm     = Axis(AxisTm.InfoAxis, AxisX);

switch obj.Format
    case 'profile'
        if ~isempty(obj.AxisY.Data)
            AxisY  = obj.AxisY.Data;
            AxisMz = obj.AxisY;
            indMZ  = findCloser(tgtMz, AxisY);
            AxisY  = AxisY(indMZ-mzWdw:indMZ+mzWdw);
            AxisMz = Axis(AxisMz.InfoAxis, AxisY);
            Data   = zeros(2*mzWdw+1, indTimeEnd-indTimeStt+1);
            
            for ii = indTimeStt:indTimeEnd
                XMS = xpend(obj, obj.ListOfScans{ii});
                Data(:, ii-indTimeStt+1) = XMS.Data(indMZ-mzWdw:indMZ+mzWdw,2);
            end
        else
            error('L32')
        end
        
    case 'centroid'
            error('L37')
end

[InfoROI.Title, errmsg] = ...
    sprintf('ROI with tm from %.2f:%.2f %s and mz=%.4f +/- %u datapoints.', ...
    timeInt(1), timeInt(2), AxisTm.Unit, tgtMz, mzWdw);
InfoROI.TagOfDts = obj.Log;
InfoROI.TgtMz    = tgtMz;
InfoROI.MzWdw    = mzWdw;
InfoROI.TInt     = timeInt;
InfoROI.AxisTm   = AxisTm;
InfoROI.AxisMZ   = AxisMz;
InfoROI.Path2Fin = obj.Path2Fin;
InfoROI.StoredData = Data;

s = ROI(InfoROI);

end