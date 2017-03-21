function plot(obj)
figure
XY = obj.TraceIn.Data;
s  = subplot(2,1, 1);
s.Parent.Name = [obj.TraceIn.Title, ' - Results AnalyzeThis'];
IdS = find(XY(:,1) <=  obj.Options.XLim(1), 1, 'last');
if isempty(IdS), IdS = 1; end
IdE = find(XY(:,1) >=  obj.Options.XLim(2), 1, 'first');
if isempty(IdE), IdE = size(XY,1); end

subplot(2, 1, 1);
XYori = obj.TraceIn.Data;
scatter(XYori(IdS:IdE,1), XYori(IdS:IdE,2), '.k');
title('Smoothing and Baseline')
XYsmo  = obj.SmoothData;
hold on
plot(XYsmo(IdS:IdE,1), XYsmo(IdS:IdE,2), 'k');
plot(XYsmo(IdS:IdE,1), obj.Baseline.vals, 'r');
hold off
xlabel([obj.TraceIn.AxisX.Label, ' / ', obj.TraceIn.AxisX.Unit]);
ylabel([obj.TraceIn.AxisY.Label, ' / ', obj.TraceIn.AxisY.Unit]);


sb = subplot(2, 1, 2);
XY = obj.XY;
plotline = plot(XY(:,1), XY(:,2), 'k');
title('Peak picking and figure of merits')
PL  = obj.PeakList;
hold on
stem(PL.Data(:,1), PL.Data(:,2), 'r')
hold off

xlabel([obj.TraceIn.AxisX.Label, ' / ', obj.TraceIn.AxisX.Unit]);
ylabel([obj.TraceIn.AxisY.Label, ' / ', obj.TraceIn.AxisY.Unit]);

c = uicontextmenu;
plotline.UIContextMenu = c;

uD = plotline.UserData;
if ~isempty(uD)
    dothat
end
% Set c to be the plot line's UIContextMenu
% Create menu items for the uicontextmenu
uimenu(c,'Label','FiguresOfMerits' ,'Callback',@getFOM);
    function getFOM(~,~)
        h = msgbox('Select peak start', '','modal');
        uiwait(h)
        axis(sb);
        [TmS, ~] = ginput(1);
        if isempty(TmS)
            return
        end
        h = msgbox('Select peak end', '','modal');
        uiwait(h)
        axis(sb);
        [TmE, ~] = ginput(1);
        if isempty(TmE)
            return
        end
        XLim(1) = min(TmS, TmE);
        XLim(2) = max(TmS, TmE);
        IdX(1)  = find(XY(:,1) <= XLim(1), 1, 'last');
        IdX(2)  = find(XY(:,1) >= XLim(2), 1, 'first');
        
        XY = obj.XY;
        uD = plotline.UserData;
        if isempty(uD)
            Id = 1;
            [Imax, idImax] = max(XY(IdX(1):IdX(2), 2));
            tmax = XY(idImax+IdX(1)-1, 1);
            x = XY(IdX(1):IdX(2),1); y = XY(IdX(1):IdX(2),2);
            M0 = trapz(x, y);
            M1 = trapz(x, x.*y)/M0;
            M2 = trapz(x, (x-M1).^2.*y)/M0;
            T     = table;
            T.Id  = Id;
            T.Imax = Imax;
            T.tmax = tmax;
            T.M0   = M0;
            T.M1   = M1;
            T.M2   = M2;
            disp(T)
        end
    end

end
