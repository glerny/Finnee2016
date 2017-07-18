%% DESCRIPTION
%
%% Copyright 
% BSD 3-Clause License
% Copyright 2016-2017 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function s = getNoise(obj, timeInt)
if length(timeInt) == 1
    error('Please choose and interval')
else
    timeStart = min(timeInt);
    timeEnd   = max(timeInt);
end

AxisX      = obj.AxisX.Data;
indTimeStt = findCloser(timeStart, AxisX);
indTimeEnd = findCloser(timeEnd, AxisX);
InfoTrc    = obj.ListOfScans{indTimeStt}.InfoTrc;
if indTimeEnd-indTimeStt+1 < 10
    error('Please choose and larger interval')
end

switch obj.Format
    case 'profile'
        if ~isempty(obj.AxisY.Data)
            AxisY = obj.AxisY.Data;
            wrkMtx = zeros(size(AxisY,1), indTimeEnd-indTimeStt+1);
            for ii = indTimeStt:indTimeEnd
                XMS = xpend(obj, obj.ListOfScans{ii});
                wrkMtx( :,ii-indTimeStt+1) = XMS.Data(:,2);
                
            end
            InfoTrc.TT = 'OTR';
            wrkMtx(wrkMtx == 0) = NaN;
            Pk2PkNoise = max(wrkMtx, [], 2) - min(wrkMtx, [], 2);
            Pk2PkNoise(sum(isnan(wrkMtx), 2) > (indTimeEnd - indTimeStt)) = nan;
        else
            error('You should first normalize all scan to a commun mz axis')
        end
        
    otherwise
        error('the dataset should be in profile mode')
end

InfoTrc.Title  = ['Peak to peak noise caluclated between', num2str(AxisX(indTimeStt), obj.AxisX.fo),...
    ' to ', num2str(AxisX(indTimeEnd), obj.AxisX.fo), ' ', obj.AxisX.Unit];
[~, partial] = decipherLog(obj.Log);
InfoTrc.FT     = partial{1};
InfoTrc.Loc    = 'inTrace';
InfoTrc.AdiPrm = {};

s = Trace(InfoTrc, [AxisY, Pk2PkNoise]);
end