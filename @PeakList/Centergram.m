function Centergram(obj)

XLim    = obj.options.XLim;
X       = obj.AxisX.Data;
infoX   = obj.AxisX.InfoAxis;
foX     = obj.AxisX.fo;
infoY   = obj.AxisY.InfoAxis;
foY     = obj.AxisY.fo;
infoZ   = obj.AxisZ.InfoAxis;
foZ     = obj.AxisZ.fo;
IdS     = find( X >= XLim(1), 1, 'first');
IdE     = find( X <= XLim(2), 1, 'last');

fig     = figure;
dcm_obj = datacursormode(fig);
set(dcm_obj,'UpdateFcn',@myupdatefcn);
pH      = plot(obj.BPP.Data(IdS:IdE,1), obj.BPP.Data(IdS:IdE,2));
pH.Tag  = 'pH';
hold on
sH      = stem(obj.FOM.Data(:,5), obj.FOM.Data(:,2));
sH.Tag  = 'sH';
title('Centergram plot')
xlabel([infoX.Label, ' / ', infoX.Unit]);
ylabel([infoZ.Label, ' / ', infoZ.Unit]);

    function txt = myupdatefcn(~,event_obj)
        % New data cursor update function
        
        tgt = get(event_obj,'Target');
        pos = get(event_obj,'Position');
        tgt.Tag
        switch tgt.Tag;
            case 'pH'
                xString = [infoX.Label, ' = ', num2str(pos(1), foX), ...
                    ' ', infoX.Unit];
                yString = [infoZ.Label, ' = ', num2str(pos(2), foZ), ...
                    ' ', infoZ.Unit];
                txt = {xString, yString};
            case 'sH'
                M1   = pos(1);
                IMax = pos(2);
                Id   = find(obj.FOM.Data(:,5) == M1 & obj.FOM.Data(:,2) == IMax);
                cFOM = obj.FOM.Data(Id,:);
                zrtL = ['PIP # ', num2str(cFOM(1))];
                fstl = ['Time @ pk max: ',  num2str(cFOM(3), foX), ...
                    ' ', infoX.Unit];
                sndl = ['Pk Area: ',  num2str(cFOM(4), foZ), ...
                    ' ', infoZ.Unit];               
                cPIP = obj.LstPIP{Id};
                AM   = sum(cPIP.Data(:,1).*cPIP.Data(:,2))/...
                    sum(cPIP.Data(:,2));
                thrl = ['Acc. Mass: ', num2str(AM, foY)];
                txt = {zrtL, fstl, sndl, thrl};

                cPIP.plot(sprintf('PIP #%u', Id));
        end
    end
end

