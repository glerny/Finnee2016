%% DESCRIPTION
%
%% Copyright
% BSD 3-Clause License
% Copyright 2016-2017 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot(obj)
% Function dealing in the plotting of the trace

fig = figure('Name', obj.FigureTitle);
dcm_obj = datacursormode(fig);
set(dcm_obj,'UpdateFcn',@myupdatefcn);
c = uicontextmenu; 
uimenu(c, 'Label','ExportObj','Callback',@doExport);
fig.UIContextMenu = c;

% Main plot
infoX = obj.AxisX.InfoAxis;
foX   = obj.AxisX.fo;
infoY = obj.AxisY.InfoAxis;
foY   = obj.AxisY.fo;

switch upper(obj.TraceType)
    case {'SEP', 'PRF', 'OTR'} 
        plot(obj.Data(:,1), obj.Data(:,2));
        title(obj.Title);
        xlabel([infoX.Label, ' / ', infoX.Unit]);
        ylabel([infoY.Label, ' / ', infoY.Unit]);
        
    case {'CTR'}
        stem(obj.Data(:,1), obj.Data(:,2), 'Marker', 'none');
        title(obj.Title);
        xlabel([infoX.Label, ' / ', infoX.Unit]);
        ylabel([infoY.Label, ' / ', infoY.Unit]);
end

    function doExport(~, ~)
        
        assignin('base', 'currentTrace', obj)
    end

    function txt = myupdatefcn(~,event_obj)
        % New data cursor update function
        
        pos = get(event_obj,'Position');
        xString = [infoX.Label, ' = ', num2str(pos(1), foX), ' ', infoX.Unit];
        yString = [infoY.Label, ' = ', num2str(pos(2), foY), ' ', infoY.Unit];
        txt = {xString, yString};
    end
end