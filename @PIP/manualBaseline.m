function PIPs = manualBaseline(obj, stringTitle)
if nargin == 1, stringTitle = ''; end
InputFig = figure(                          ...
    'Visible'          , 'on'             , ...
    'Name'             , stringTitle      , ...
    'Toolbar'          , 'none'           , ...
    'MenuBar'          , 'none'           , ...
    'Units'            , 'normalized'     , ...
    'WindowStyle'      , 'modal'         , ...
    'UserData'         , obj              , ...
    'Resize'           , 'on');

AxisH1 = axes(InputFig ,...
    'Units'            , 'normalized'     , ...
    'OuterPosition'    , [0 0 2/3 1]);

PushH{1} = uicontrol(InputFig ,...
    'Style'   , 'pushbutton', ...
    'Visible' , 'on'        , ...
    'Units'   , 'normalized', ...
    'String'  , 'Select peaks pts',...
    'Tag'     , 'PushH1',...
    'Callback', @pushMe     , ...
    'Position', [2/3 0 1/3 1]);
wPB1            = 1.1*PushH{1}.Extent(3);
hPB1            = 1.2*PushH{1}.Extent(4);

x = 2/3;
y = 1-3*hPB1;

PushH{2} = uicontrol(InputFig ,...
    'Style'   , 'pushbutton', ...
    'Visible' , 'on'        , ...
    'Units'   , 'normalized', ...
    'String'  , 'Reset',...
    'Tag'     , 'PushH2',...
    'Callback', @pushMe     , ...
    'Position', [x y wPB1 hPB1]);
y = y-1.2*hPB1;

PushH{3} = uicontrol(InputFig ,...
    'Style'   , 'pushbutton', ...
    'Visible' , 'on'        , ...
    'Units'   , 'normalized', ...
    'String'  , 'Crop',...
    'Tag'     , 'PushH3',...
    'Callback', @pushMe     , ...
    'Position', [x y wPB1 hPB1]);
y = y-1.2*hPB1;


PushH{4} = uicontrol(InputFig ,...
    'Style'   , 'pushbutton', ...
    'Visible' , 'on'        , ...
    'Units'   , 'normalized', ...
    'String'  , 'Done',...
    'Tag'     , 'PushH4',...
    'Callback', @pushMe     , ...
    'Position', [x y wPB1 hPB1]);
y = y-1.2*hPB1;

PushH{1}.Position = [x y wPB1 hPB1];

y = y-1.2*hPB1;

uicontrol(InputFig ,...
    'Style'              , 'Text'      , ...
    'Visible'            , 'on'        , ...
    'Units'              , 'normalized', ...
    'String'             , 'Degree for the polynomial:' ,...
    'HorizontalAlignment', 'left', ...
    'Position'           , [x y 2*wPB1 hPB1]);

y = y-1.05*hPB1;

EditH{1} = uicontrol(InputFig ,...
    'Style'              , 'Edit'      , ...
    'Visible'            , 'on'        , ...
    'Units'              , 'normalized', ...
    'String'             , '3' ,...
    'HorizontalAlignment', 'left', ...
    'Callback', @changeMe, ...
    'Position'           , [x y wPB1 hPB1]);



IdMin    = 1;
IdMax    = length(obj.x);
Id4Bas   = true(size(obj.x));
n        = 3;
plotPIP

if ishghandle(InputFig)
    uiwait(InputFig);
end

if ~ishghandle(InputFig)
    disp('done'); % if no handles and answer 0 either esc or stop the fuck
end

    function plotPIP
        obj.Data = sortrows(obj.Data, 3);
        Id2Rem   = obj.Data(:,3) < IdMin - 1+ obj.IdS | obj.Data(:,3) > IdMax - 1+ obj.IdS;
        IdMin    = 1;
        IdMax    = length(obj.x);
        obj.Data(Id2Rem, :) = [];
        obj.IdS  = min(obj.Data(:,3));
        obj.x    = obj.x(min(obj.Data(:,3)) - obj.IdS +1:max(obj.Data(:,3)) - obj.IdS +1);
        Id4Bas   = Id4Bas(min(obj.Data(:,3)) - obj.IdS +1:max(obj.Data(:,3)) - obj.IdS +1);
        
        p = polyfitweighted(obj.x(Id4Bas),obj.y(Id4Bas),n);
        plot(AxisH1, obj.x, obj.y, 'k');
        
        hold on
        plot(AxisH1, obj.x(Id4Bas), obj.y(Id4Bas), 'ko');
        plot(AxisH1, obj.x, polyval(p, obj.x), 'r');
        
        
        title('Manual baseline correction')
        xlabel([obj.AxisX.Label, ' / ', obj.AxisX.Unit]);
        ylabel([obj.AxisZ.Label, ' / ', obj.AxisZ.Unit]);
        
        hold off
    end

    function pushMe(source, ~)
        switch source.Tag
            case 'PushH1'
                
                h = msgbox('Select the peak start', 'Correct','modal');
                uiwait(h)
                figure(InputFig);
                [fmin, ~] = ginput(1);
                if isempty(fmin)
                    return
                end
                Id1 = find(obj.x <= fmin, 1, 'last');
                if isempty(Id1)
                    Id1 = 1;
                end
                h = msgbox('Select the peak end', 'Correct','modal');
                uiwait(h)
                figure(InputFig);
                [fmax, ~] = ginput(1);
                if isempty(fmax)
                    Id1 = 1;
                    return
                end
                Id2 = find(obj.x >= fmax, 1, 'first');
                if isempty(Id2)
                    Id2 = length(obj.x);
                end

                Id4Bas(Id1:Id2) = false;
                plotPIP
                
            case 'PushH2'
                obj = InputFig.UserData;
                IdMin    = 1;
                IdMax    = length(obj.x);
                Id4Bas   = true(size(obj.x));
                plotPIP
                
            case 'PushH3'
                h = msgbox('Select the starting points', 'Correct','modal');
                uiwait(h)
                figure(InputFig);
                [fmin, ~] = ginput(1);
                if isempty(fmin)
                    return
                end
                IdMin = find(obj.x <= fmin, 1, 'last');
                if isempty(IdMin)
                    IdMin = 1;
                end
                h = msgbox('Select the ending points', 'Correct','modal');
                uiwait(h)
                figure(InputFig);
                [fmax, ~] = ginput(1);
                if isempty(fmax)
                    IdMin = 1;
                    return
                end
                IdMax = find(obj.x >= fmax, 1, 'first');
                if isempty(IdMax)
                    IdMax = length(obj.x);
                end
                plotPIP
                
            case 'PushH4'
                p = polyfitweighted(obj.x(Id4Bas),obj.y(Id4Bas),n);
                for ii = 1:size(obj.Data, 1)
                    x_ = obj.x(obj.Data(ii,3)-obj.IdS+1);
                    obj.Data(ii,2) = obj.Data(ii,2) - polyval(p, x_);
                end
                
                [~, w] = polyfitweighted(obj.x, obj.y, 0);
                Id4Pks = obj.Data(:,2) > 3*std(obj.Data(:,2), w);
                PIPs = {};
                
                while 1
                    if ~any(Id4Pks)
                        break
                    end
                    
                    IdP = find(Id4Pks, 1, 'first');
                    IdS = find(obj.Data(1:IdP-1,2) <= 0, 1, 'last');
                    if isempty(IdS)
                        IdS = 1;
                    end
                    IdE = find(obj.Data(IdP+1:end,2) <= 0, 1, 'first');
                    if isempty(IdE)
                        IdE = size(obj.Data, 1);
                    else
                        IdE = IdE + IdP;
                    end
                    Id4Pks(IdS:IdE) = false;       
                    
                    if IdE - IdS >= 5
                        nPIP = obj;
                        nPIP.Data = obj.Data(IdS:IdE, :);
                        nPIP.Data(:,2) = round(nPIP.Data(:,2));
                        nPIP.Data(nPIP.Data(:,2) < 0,2) = 0;
                        nPIP.x = obj.x(IdS:IdE);
                        nPIP.IdS = obj.IdS + IdS -1;
                    end
                    PIPs{end+1} = nPIP;
                end
                
                close(InputFig)
        end
    end

    function changeMe(source, ~)
        n = round(str2double(source.String));
        if isnan(n)
            n = 3;
        end
        
        if n > 3
            n = 3;
        end
        source.String = num2str(n);
        plotPIP
    end

end
