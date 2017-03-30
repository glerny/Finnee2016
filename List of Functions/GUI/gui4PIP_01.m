%% DESCRIPTION
% GUI for @PeakList
%% Copyright
% BSD 3-Clause License
% Copyright 2016-2017 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PIPs, ExitFlag, stDim] = gui4PIP_01(cPIP, str2fig, stDim)

ExitFlag  = -2; %-2 Abord without saving; -1 stop and save ; 0 reset last PIP; 1 discard PIP; use PIPs
PIPs{1}   = cPIP;

%% 1. Filling the GUI
% 1.1 Dimension
cor.baseline = [];
cor.split    = [];
cor.crop     = []; 
X            = cPIP.AxisX.Data;

InputFig = figure(                          ...
    'Visible'          , 'on'             , ...
    'Name'             , str2fig{1}       , ...
    'Toolbar'          , 'none'           , ...
    'MenuBar'          , 'none'           , ...
    'Units'            , 'normalized'     , ...
    'Position'         , stDim            , ...
    'Tag'              , 'tagFig'         , ...
    'HandleVisibility' , 'callback'       , ...
    'WindowStyle'      , 'normal'         , ...
    'UserData'         , cor              , ...
    'Resize'           , 'on');


InputFig.Units = 'pixels';
insDim = InputFig.Position; %[x y width height]

queryH = uicontrol (InputFig, ...
    'Style'  , 'Text'   , ...
    'Units'  , 'pixels' , ...
    'Visible', 'off'    , ...
    'String' ,'Average PIPs'' Variance (min^2)') ;
dftTX(3) = queryH.Extent(3);
dftTX(4) = queryH.Extent(4);

% Query longer PB
queryH = uicontrol (InputFig, ...
    'Style'  , 'pushbutton' , ...
    'Units'  , 'pixels'     , ...
    'Visible', 'off'        , ...
    'String' , 'Discard');
dftPB(3) = queryH.Extent(3);
dftPB(4) = queryH.Extent(4);



% Query longer text
queryH = uicontrol (InputFig, ...
    'Style'  , 'Text'   , ...
    'Units'  , 'pixels' , ...
    'Visible', 'off'    , ...
    'String' ,[str2fig{2}, '; ', str2fig{3}, ...
    ';  ', str2fig{4}]);
minWdh  = queryH.Extent(3);

% 1.2. Place the elements
needHeight   = 20*dftPB(4);
if insDim(4) < needHeight, insDim(4) = needHeight; end
needWidth    = minWdh;
if insDim(3) < needWidth, insDim(3) = needWidth; end
InputFig.Position = insDim;

% text and slider
x = dftTX(4);
y = insDim(4) - 1.5*dftTX(4);

uicontrol (InputFig               , ...
    'Style'   , 'text'                      , ...
    'Position', [x y  insDim(3)  dftTX(4)]  , ...
    'Tag'     , 'textH1'                    , ...
    'FontWeight', 'bold'                    , ...
    'String'  ,[str2fig{2}, '; ', str2fig{3}, ...
    ';  ', str2fig{4}]);

x = insDim(3) - 3*dftPB(3);
y = insDim(4) - 3*dftTX(4);
uicontrol (InputFig             , ...
    'Style'   , 'pushbutton'             , ...
    'Position', [x y  dftPB(3)  dftPB(4)], ...
    'Tag'     , 'PsBt1'                  , ...
    'Callback', @pushMe                  , ...
    'String'  , 'Split');
x = x + 1.2*dftPB(3);

uicontrol (InputFig             , ...
    'Style'   , 'pushbutton'             , ...
    'Position', [x y  dftPB(3)  dftPB(4)], ...
    'Tag'     , 'PsBt2'                  , ...
    'Callback', @pushMe                  , ...
    'String'  , 'Crop');
y = y - 1.2*dftPB(4);
x = x - 1.2*dftPB(3);

uicontrol (InputFig                 , ...
    'Style'   , 'pushbutton'                 , ...
    'Position', [x y  2.2*dftPB(3)  dftPB(4)], ...
    'Tag'     , 'PsBt3'                      , ...
    'Callback', @pushMe                      , ...
    'String'  , 'Baseline cor.');
y = y - 1.2*dftPB(4);
x = x + 0.6*dftPB(3);

uicontrol (InputFig             , ...
    'Style'   , 'pushbutton'             , ...
    'Position', [x y  dftPB(3)  dftPB(4)], ...
    'Tag'     , 'PsBt4'                  , ...
    'Callback', @pushMe                  , ...
    'String'  , 'Reset');
y = y - 3.6*dftPB(4);
x = x - 0.6*dftPB(3);

uicontrol (InputFig             , ...
    'Style'   , 'pushbutton'             , ...
    'Position', [x y  dftPB(3)  dftPB(4)], ...
    'Tag'     , 'PsBt5'                  , ...
    'Callback', @pushMe                  , ...
    'String'  , 'Keep');
x = x + 1.2*dftPB(3);

uicontrol (InputFig             , ...
    'Style'   , 'pushbutton'             , ...
    'Position', [x y  dftPB(3)  dftPB(4)], ...
    'Tag'     , 'PsBt6'                  , ...
    'Callback', @pushMe                  , ...
    'String'  , 'Discard');
y = y - 1.2*dftPB(4);
x = x - 1.2*dftPB(3);

uicontrol (InputFig                 , ...
    'Style'   , 'pushbutton'                 , ...
    'Position', [x y  2.2*dftPB(3)  dftPB(4)], ...
    'Tag'     , 'PsBt7'                      , ...
    'Callback', @pushMe                      , ...
    'String'  , 'Peak deconvol.');
y = y - 3.6*dftPB(4);
x = x + 0.6*dftPB(3);

uicontrol (InputFig                     , ...
    'Style'   , 'pushbutton'                     , ...
    'Position', [x y  1.2*dftPB(3)  1.2*dftPB(4)], ...
    'Tag'     , 'PsBt8'                          , ...
    'Callback', @pushMe                          , ...
    'FontWeight', 'bold'                         , ...
    'ForegroundColor', 'r'                       ,...
    'String'  , 'Finish');
y = y - 1.5*dftPB(4);

uicontrol (InputFig                     , ...
    'Style'   , 'pushbutton'                     , ...
    'Position', [x y  1.2*dftPB(3)  1.2*dftPB(4)], ...
    'Tag'     , 'PsBt9'                          , ...
    'Callback', @pushMe                          , ...
    'FontWeight', 'bold'                         , ...
    'ForegroundColor', 'r'                       ,...
    'String'  , 'Abord');

wRem = insDim(3) - 3*dftPB(3);
hRem = insDim(4) - 2*dftTX(4);

% Place axes
AxisH1 = axes(InputFig        , ...
    'Units'        , 'pixels' , ...
    'Tag'          , 'AxisH1' , ...
    'OuterPosition'     , [0 0 wRem hRem]); % main axis


% 1.3. activate the figure and resize function

InputFig.Visible        = 'on';
InputFig.SizeChangedFcn =  @scf;
plotPIP(cPIP);
if ishghandle(InputFig)
    uiwait(InputFig);
end

if ~ishghandle(InputFig)
    disp('done'); % if no handles and answer 0 either esc or stop the fuck
end
    function plotPIP(cPIP)
        yyaxis(AxisH1, 'left')
        plot(AxisH1, cPIP.x, cPIP.y);
        
        yyaxis(AxisH1, 'right')
        scatter(AxisH1, cPIP.x(cPIP.Data(:,3) - cPIP.IdS+1), cPIP.Data(:,1), 'o');
        %AxisH{1}.Title.String = sprintf('PIP #%u', FOM(stID, 1));
        AxisH1.XAxis.Label.String  = sprintf('%s / %s', ...
            cPIP.AxisX.Label, cPIP.AxisX.Unit);
        AxisH1.YAxis(1).Label.String  = sprintf('%s / %s', ...
            cPIP.AxisZ.Label, cPIP.AxisZ.Unit);
        AxisH1.YAxis(2).Label.String  = sprintf('%s / %s', ...
            cPIP.AxisY.Label, cPIP.AxisY.Unit);
        drawnow()
        
    end

    function scf(~,~)
        fig = gcbo;
        insDim = fig.Position;
        if insDim(4) < needHeight, insDim(4) = needHeight; end
        if insDim(3) < needWidth, insDim(3) = needWidth; end
        fig.Position = insDim;
        
        % text and slider
        x = dftTX(4);
        y = insDim(4) - 1.5*dftTX(4);
        u = findobj(gcbo,'Tag','textH1');
        u.Position =  [x y  insDim(3)  dftTX(4)];
        
        x = insDim(3) - 3*dftPB(3);
        y = insDim(4) - 3*dftTX(4);
        
        u = findobj(gcbo,'Tag','PsBt1');
        u.Position =  [x y  dftPB(3)  dftPB(4)];
        x = x + 1.2*dftPB(3);
        
        
        u = findobj(gcbo,'Tag','PsBt2');
        u.Position =  [x y  dftPB(3)  dftPB(4)];
        y = y - 1.2*dftPB(4);
        x = x - 1.2*dftPB(3);
        
        u = findobj(gcbo,'Tag','PsBt3');
        u.Position =  [x y  2.2*dftPB(3)  dftPB(4)];
        y = y - 1.2*dftPB(4);
        x = x + 0.6*dftPB(3);
        
        u = findobj(gcbo,'Tag','PsBt4');
        u.Position =  [x y  dftPB(3)  dftPB(4)];
        y = y - 3.6*dftPB(4);
        x = x - 0.6*dftPB(3);
        
        u = findobj(gcbo,'Tag','PsBt5');
        u.Position =  [x y  dftPB(3)  dftPB(4)];
        x = x + 1.2*dftPB(3);
        
        u = findobj(gcbo,'Tag','PsBt6');
        u.Position =  [x y  dftPB(3)  dftPB(4)];
        y = y - 1.2*dftPB(4);
        x = x - 1.2*dftPB(3);
        
        u = findobj(gcbo,'Tag','PsBt7');
        u.Position =  [x y  2.2*dftPB(3)  dftPB(4)];
        y = y - 3.6*dftPB(4);
        x = x + 0.6*dftPB(3);
        
        u = findobj(gcbo,'Tag','PsBt9');
        u.Position =  [x y  1.2*dftPB(3)  1.2*dftPB(4)];
        y = y - 1.5*dftPB(4);
        
        u = findobj(gcbo,'Tag','PsBt8');
        u.Position =  [x y  1.2*dftPB(3)  1.2*dftPB(4)];
        y = y - 1.5*dftPB(4);
        
        wRem = insDim(3) - 3*dftPB(3);
        hRem = insDim(4) - 2*dftTX(4);
        
        u = findobj(gcbo,'Tag','AxisH1');
        if ~ isempty(u)
            u.OuterPosition = [0 0 wRem hRem]; % main axis
        end
        
        InputFig.Units = 'normalized';
        stDim = InputFig.Position; %[x y width height]
        InputFig.Units = 'pixels';
    end

    function pushMe(source, ~)
        cor = InputFig.UserData;
        
        switch source.Tag
            case 'PsBt1'
                % split
                cor.split = []; 
                figure(InputFig);
                [TgX, ~] = ginput(1);
                if ~isempty(TgX)
                    Id = find(X >= TgX, 1, 'first');
                    if ~isempty(Id)
                        cor.split = Id; 
                    end 
                end
                
            case 'PsBt2'
                % crop
                cor.crop = [];
                figure(InputFig);
                [TgX, ~] = ginput(2);
                if length(TgX) == 2
                    TgX = sort(TgX);
                    Id(1) = find(X >= TgX(1), 1, 'first');
                    Id(2) = find(X <= TgX(2), 1, 'last');
                    cor.crop = [Id(1) Id(2)];
                end
                
            case 'PsBt3'
                % baseline correctiom
                cor.baseline = [];            
                figure(InputFig);
                [TgX, TgY] = ginput;
                if length(TgX) > 2
                end
                
            case 'PsBt4'
                delete(InputFig)
            case 'PsBt5'
                delete(InputFig)
            case 'PsBt6'
                delete(InputFig)
            case 'PsBt7'
                delete(InputFig)
            case 'PsBt8'
                delete(InputFig)
            case 'PsBt9'
                ExitFlag = -1;
                delete(InputFig)
        end
    end
end
























