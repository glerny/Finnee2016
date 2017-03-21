function plot(obj)

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


subplot(2, 1, 2);
XY = obj.XY;
plot(XY(:,1), XY(:,2), '--k');
title('Peak picking and peak fitting')
PL  = obj.PeakList;
hold on
stem(PL.Data(:,1), PL.Data(:,2), 'r')
hold off
dcP = obj.Fitting;
if ~isempty(dcP)
    hold on
    for ii = 1:length(dcP.Int)
        plot(XY(:,1), dcP.Int(ii)* dcP.model(ii, :), 'b');
    end
    hold off
end
    
xlabel([obj.TraceIn.AxisX.Label, ' / ', obj.TraceIn.AxisX.Unit]);
ylabel([obj.TraceIn.AxisY.Label, ' / ', obj.TraceIn.AxisY.Unit]);
end
