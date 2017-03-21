%% DESCRIPTION
%
%% Copyright 
% BSD 3-Clause License
% Copyright 2016-2017 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function s = getSpectra(obj, timeInt)
if length(timeInt) == 1
    [timeStart, timeEnd] = deal(timeInt(1));
else
    timeStart = min(timeInt);
    timeEnd   = max(timeInt);
end
AxisX       = obj.AxisX.Data;
indTimeStt = findCloser(timeStart, AxisX);
indTimeEnd = findCloser(timeEnd, AxisX);

switch obj.Format
    case 'profile'
        data(:,1) = obj.AxisY.Data;
        data(:,2) = 0;
        for ii = indTimeStt:indTimeEnd
            XMS = xpend(obj, obj.ListOfScans{ii});
            data(:,2) = data(:,2) + XMS.Data(:,2);
        end
        InfoTrc.TT = 'PRF';
        
    case 'centroid'
        data = [];
        for ii = indTimeStt:indTimeEnd
            MS = obj.ListOfScans{ii}.Data;
            data = [data; MS];
            data = sortrows(data, 1);
        end
        InfoTrc.TT = 'CTR';
end

InfoTrc.Title  = ['MS scan from ', num2str(AxisX(indTimeStt), obj.AxisX.fo),...
    ' to ', num2str(AxisX(indTimeEnd), obj.AxisX.fo), ' ', obj.AxisX.Unit];
strLog         = decipherLog(obj.Log, 1);
InfoTrc.FT     = strLog{1};
InfoTrc.AxisX   = Axis(obj.AxisY.InfoAxis);
InfoTrc.AxisY   = Axis(obj.AxisZ.InfoAxis);
InfoTrc.Loc    = 'inTrace';
InfoTrc.AdiPrm = {};

s = Trace(InfoTrc, data);
end