function Centergram(obj)

if strcmp(obj.Type, 'singleton')
    XLim    = [0 inf];
    X       = obj.AxisX{1}.Data;
    infoX   = obj.AxisX{1}.InfoAxis;
    foX     = obj.AxisX{1}.fo;
    infoY   = obj.AxisY{1}.InfoAxis;
    foY     = obj.AxisY{1}.fo;
    infoZ   = obj.AxisZ{1}.InfoAxis;
    foZ     = obj.AxisZ{1}.fo;
    IdS     = find( X >= XLim(1), 1, 'first');
    IdE     = find( X <= XLim(2), 1, 'last');
    
    % Drawing the figure
    fig     = figure(...
        'Name' , obj.Path2Fin{1},...
        'Units', 'normalized');
    
    dcm_obj = datacursormode(fig);
    set(dcm_obj,'UpdateFcn',@myupdatefcn4s);
    pH      = plot(obj.BPP{1}.Data(IdS:IdE,1), obj.BPP{1}.Data(IdS:IdE,2), 'k');
    pH.Tag  = 'pH';
    hold on
    sH      = stem(obj.FOM{1}.Data(:,5), obj.FOM{1}.Data(:,2), 'r');
    sH.Tag  = 'sH';
    title('Centergram plot')
    xlabel([infoX.Label, ' / ', infoX.Unit]);
    ylabel([infoZ.Label, ' / ', infoZ.Unit]);
    
    ToglH1 = uicontrol(fig ,...
        'Style'   , 'togglebutton', ...
        'Visible' , 'on'        , ...
        'Units'   , 'normalized', ...
        'String'  , 'popup PIP on',...
        'Tag'     , 'ToglH1    ',...
        'Callback', @pushMe     , ...
        'Position', [2/3 0 1/3 1]);
    wPB1            = ToglH1.Extent(3);
    hPB1            = ToglH1.Extent(4);
    ToglH1.Position = [1-1.3*wPB1 1-1.3*hPB1  1.2*wPB1  1.2*hPB1];
    
elseif strcmp(obj.Type, 'replicates')
    lineFrm ={'-k'; '--k'; ':k'; '-.k'; '-b'; '--b'; ':b'; '-.b';...
        '-m'; '--m'; ':m'; '-.m'};
    XLim    = [0 inf];
    rpts    = length(obj.Path2Fin);
    X       = obj.AxisX{1}.Data;
    infoX   = obj.AxisX{1}.InfoAxis;
    foX     = obj.AxisX{1}.fo;
    infoY   = obj.AxisY{1}.InfoAxis;
    foY     = obj.AxisY{1}.fo;
    infoZ   = obj.AxisZ{1}.InfoAxis;
    foZ     = obj.AxisZ{1}.fo;
    IdS     = find( X >= XLim(1), 1, 'first');
    IdE     = find( X <= XLim(2), 1, 'last');
    
    stringName = '';
    Label      = {};
    for ii = 1:rpts-1
        stringName = [stringName, obj.Path2Fin{ii}, ' & '];
        Label{ii}  = ['Replicate: #', num2str(ii)];
    end
    
    stringName = [stringName, obj.Path2Fin{rpts}];
    Label{rpts}  = ['Replicate: #', num2str(rpts)];
    
    fig     = figure(...
        'Name', stringName,...
        'Units', 'normalized');
    
    ToglH1 = uicontrol(fig ,...
        'Style'   , 'togglebutton', ...
        'Visible' , 'on'        , ...
        'Units'   , 'normalized', ...
        'String'  , 'popup PIP on',...
        'Tag'     , 'ToglH1    ',...
        'Callback', @pushMe     , ...
        'Position', [0 0 1 1]);
    wPB1            = ToglH1.Extent(3);
    hPB1            = ToglH1.Extent(4);
    ToglH1.Position = [1-1.3*wPB1 1-1.3*hPB1  wPB1  hPB1];
    
    dcm_obj = datacursormode(fig);
    set(dcm_obj,'UpdateFcn',@myupdatefcn4r);
    hold on
    for ii = 1:rpts
        pH{ii} = plot(obj.BPP{ii}.Data(:,1), obj.BPP{ii}.Data(:,2), lineFrm{ii});
        pH{ii}.Tag  = 'pH';
    end
    
    data4ctg.intMax = [];
    data4ctg.centre = [];
    for ii = 1:rpts
        data4ctg.intMax = [data4ctg.intMax, obj.FOM{ii}.Data(:,2)];
        data4ctg.centre = [data4ctg.centre, obj.FOM{ii}.Data(:,5)];
    end
    legend(Label)
    
    sH      = stem(mean(data4ctg.centre,2), mean(data4ctg.intMax,2), 'r');
    sH.Tag  = 'sH';
    title('Centergram plot')
    xlabel([infoX.Label, ' / ', infoX.Unit]);
    ylabel([infoZ.Label, ' / ', infoZ.Unit]);
    
    
    
end

    function txt = myupdatefcn4s(~,event_obj)
        % New data cursor update function
        
        tgt = get(event_obj,'Target');
        pos = get(event_obj,'Position');
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
                Id   = find(obj.FOM{1}.Data(:,5) == M1 & obj.FOM{1}.Data(:,2) == IMax);
                cFOM = obj.FOM{1}.Data(Id,:);
                zrtL = ['PIP # ', num2str(cFOM(1))];
                fstl = ['Time @ pk max: ',  num2str(cFOM(3), foX), ...
                    ' ', infoX.Unit];
                sndl = ['Pk Area: ',  num2str(cFOM(4), foZ), ...
                    ' ', infoZ.Unit];
                thrl = ['Acc. Mass: ',  num2str(cFOM(10), foY)];
                txt = {zrtL, fstl, sndl, thrl};
                
                if ToglH1.Value == 0
                    cPIP = obj.LstPIP{1}{Id};
                    cPIP.plot(sprintf('PIP #%u', Id));
                end
        end
    end


    function txt = myupdatefcn4r(~,event_obj)
        % New data cursor update function
        
        tgt = get(event_obj,'Target');
        pos = get(event_obj,'Position');
        switch tgt.Tag;
            case 'pH'
                Title   = tgt.DisplayName;
                xString = [infoX.Label, ' = ', num2str(pos(1), foX), ...
                    ' ', infoX.Unit];
                yString = [infoZ.Label, ' = ', num2str(pos(2), foZ), ...
                    ' ', infoZ.Unit];
                txt = {Title, xString, yString};
            case 'sH'
                M1   = pos(1);
                IMax = pos(2);
                Id   = find(mean(data4ctg.centre,2) == M1 & mean(data4ctg.intMax,2) == IMax);
                cFOM = [];
                for ii = 1:rpts
                    cFOM = [cFOM; obj.FOM{ii}.Data(Id,:)];
                end
                
                zrtL = ['PIP # ', num2str(cFOM(1))];
                fo = ['%.', num2str(signFig(std(cFOM(:,5)))), 'f'];
                fstl = ['Centers: ',  num2str(mean(cFOM(:,5)), fo), ...
                    ' (', num2str(std(cFOM(:,5)), fo), ...
                    ' n=', num2str(rpts),') ', infoX.Unit];
                fo = ['%.', num2str(signFig(std(cFOM(:,4)))), 'f'];
                sndl = ['Pk Area: ',  num2str(mean(cFOM(:,4)), fo), ...
                    ' (', num2str(std(cFOM(:,4)), fo), ...
                    ' n=', num2str(rpts),') ', infoZ.Unit];
                fo = ['%.', num2str(signFig(std(cFOM(:,10)))), 'f'];
                thrl = ['Acc. Mass: ',  num2str(mean(cFOM(:,10)), fo), ...
                    ' (', num2str(std(cFOM(:,10)), fo), ...
                    ' n=', num2str(rpts),') '];
                txt = {zrtL, fstl, sndl, thrl};
                
                
                if ToglH1.Value == 0
                    for ii = 1:rpts
                        cPIP{ii} = obj.LstPIP{ii}{Id};
                    end
                    assignin('base', 'cPIP', cPIP)
                    
                    mergePIP(cPIP)
                end
        end
    end

    function pushMe(source, ~)
        if source.Value == 1
            source.String = 'popup PIP off';
        else
            source.String = 'popup PIP on';
        end
    end

end

