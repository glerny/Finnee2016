%% DESCRIPTION
%
%% Copyright
% BSD 3-Clause License
% Copyright 2016-2017 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function s = getProfile(obj, mzInt)
if length(mzInt) == 1
    [mzStart, mzEnd] = deal(mzInt(1));
else
    mzStart = min(mzInt);
    mzEnd   = max(mzInt);
end
AxisX     = obj.AxisY.Data;

prof(:,1) = obj.AxisX.Data;
prof(:,2) = 0;

switch obj.Format
    case 'profile'
        indMZStt = findCloser(mzStart, AxisX);
        indMZEnd = findCloser(mzEnd, AxisX);

        h = waitbar(0,'Calculating profile, please wait');
        for ii = 1:length(prof(:,1))
            waitbar(ii/length(prof(:,1)))
            XMS = xpend(obj, obj.ListOfScans{ii});
            prof(ii,2) = sum(XMS.Data(indMZStt:indMZEnd, 2));
        end
        InfoTrc.Title = ['Extracted ion Profiles from ', ...
            num2str(AxisX(indMZStt), obj.AxisY.fo),...
            ' to ', num2str(AxisX(indMZEnd), obj.AxisY.fo), ...
            ' ', obj.AxisY.Unit];
        
    case 'centroid'
        h = waitbar(0,'Calculating profile, please wait');
        for ii = 1:length(prof(:,1))
            waitbar(ii/length(prof(:,1)))
            MS = obj.ListOfScans{ii}.Data;
            InfoTrc.Title = ['Extracted ion Profiles from ', ...
                num2str(mzStart, obj.AxisY.fo),...
                ' to ', num2str(mzEnd, obj.AxisY.fo),...
                ' ', obj.AxisY.Unit];
            if isempty(MS)
                prof(ii,2) = 0;
                continue
            end
                
            ind2keep = MS(:,1) >= mzStart & MS(:,1) <= mzEnd;
            if any(ind2keep)
                prof(ii,2) = sum(MS(ind2keep, 2));
            else
                prof(ii,2) = 0;
            end
        end
end

try close(h); catch, end
InfoTrc.TT     = 'SEP';
strLog         = decipherLog(obj.Log, 1);
InfoTrc.FT     = strLog{1};
InfoTrc.AxisX   = Axis(obj.AxisX.InfoAxis);
InfoTrc.AxisY   = Axis(obj.AxisZ.InfoAxis);
InfoTrc.Loc    = 'inTrace';
InfoTrc.AdiPrm = {};

s = Trace(InfoTrc, prof);
end
