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
AxisX      = obj.AxisX.Data;
indTimeStt = findCloser(timeStart, AxisX);
indTimeEnd = findCloser(timeEnd, AxisX);
InfoTrc    = obj.ListOfScans{indTimeStt}.InfoTrc;

switch obj.Format
    case 'profile'
        if ~isempty(obj.AxisY.Data)
            data(:,1) = obj.AxisY.Data;
            data(:,2) = 0;
            for ii = indTimeStt:indTimeEnd
                XMS = xpend(obj, obj.ListOfScans{ii});
                data(:,2) = data(:,2) + XMS.Data(:,2);
            end
            InfoTrc.TT = 'PRF';
        else
            if indTimeStt == indTimeEnd
                data = obj.ListOfScans{indTimeStt}.Data;
            else
                txt = ('warning. \nThere is no master mz axis in dataset %i. \nA provisory one will be created.\nIt may take time\n');
                Log = decipherLog(obj.Log);
                warning(txt,  Log{1}.dtsId)
                data(:,1) = obj.ListOfScans{indTimeStt}.extrapolMZ;
                data(:,2) = 0;

                for ii = indTimeStt:indTimeEnd
                    XMS = extrapol2axis(obj.ListOfScans{ii}.Data, data(:,1));
                    data(:,2) = data(:,2) + XMS(:,2);
                end
            end
            InfoTrc.TT = 'PRF';
        end
        
    case 'centroid'
        data = [];
        for ii = indTimeStt:indTimeEnd
            MS = obj.ListOfScans{ii}.Data;
            data = [data; MS]; %#ok<AGROW>
            data = sortrows(data, 1);
        end
        InfoTrc.TT = 'CTR';
end

InfoTrc.Title  = ['MS scan from ', num2str(AxisX(indTimeStt), obj.AxisX.fo),...
    ' to ', num2str(AxisX(indTimeEnd), obj.AxisX.fo), ' ', obj.AxisX.Unit];
[~, partial] = decipherLog(obj.Log);
InfoTrc.FT     = partial{1};
InfoTrc.Loc    = 'inTrace';
InfoTrc.AdiPrm = {};

s = Trace(InfoTrc, data);
end