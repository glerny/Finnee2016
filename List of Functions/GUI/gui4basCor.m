%% DESCRIPTION
% GUI for @Finnee\BaselineCorrection
%
%% Copyright 
% BSD 3-Clause License
% Copyright 2016-2017 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function [method, ExitFlag] = gui4basCor( AxisX, AxisZ, prof2cor, XLim)

%%%%%% LIST OF METHODS AND INITIAL METHOD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
listString = {'PolyFit'; 'ArPLS'; 'ArPLS2'};
method     = 'ArPLS:10E4:10E-2';
lastLS     = {'PolyFit:2'; 'ArPLS:10E4:10E-2'; 'ArPLS2:10E4'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ExitFlag = -1; %-1 Abord, 0, finished

%% 1. Filling the GUI
% 1.1 Dimension

InputFig = figure(                          ...
    'Visible'          , 'on'             , ...
    'Name'             , 'Baseline method', ...
    'Toolbar'          , 'none'           , ...
    'MenuBar'          , 'none'           , ...
    'Units'            , 'normalized'     , ...
    'Position'         , [0, 0, 0.9, 0.9] ,...
    'UserData'         , false            , ...
    'Tag'              , 'tagFig'         , ...
    'HandleVisibility' , 'callback'       , ...
    'WindowStyle'      , 'normal'         , ...
    'Resize'           , 'on');

InputFig.Units = 'pixels';
insDim = InputFig.Position; %[x y width height]
posAx1 = []; posAx2 = []; posAx3 = [];

% blck push buttin to get the offset
queryH = uicontrol (InputFig, ...
    'Style'  , 'pushbutton' , ...
    'Units'  , 'pixels'     , ...
    'Visible', 'off'        , ...
    'String' , ' '         );
mgWidth  = (queryH.Extent(3)/4);
mgHeight = (queryH.Extent(4)/4);
dftPB    = [0 0 0 0]; % default position for push button
dftTX    = [0 0 0 0]; % default position for text

% Query longer text
queryH = uicontrol (InputFig, ...
    'Style'  , 'Text'   , ...
    'Units'  , 'pixels' , ...
    'Visible', 'off'    , ...
    'String' ,'Order (0- cst, 1- linear,...):') ;
dftTX(3) = max(dftTX(3), queryH.Extent(3) + 2*mgWidth);
dftTX(4) = max(dftTX(4), queryH.Extent(4) + 2*mgHeight);

% Query longer PB
queryH = uicontrol (InputFig, ...
    'Style'  , 'pushbutton' , ...
    'Units'  , 'pixels'     , ...
    'Visible', 'off'        , ...
    'String' , 'Next Profile');
dftPB(3) = max(dftPB(3), queryH.Extent(3) + 2*mgWidth);
dftPB(4) = max(dftPB(4), queryH.Extent(4) + 2*mgHeight);


% 1.2. Place the elements
needHeight   = 9*dftPB(4) + 14*mgHeight;
if insDim(4) < needHeight, insDim(4) = needHeight; end
needWidth    = dftTX(3) + 4*mgWidth;
if insDim(3) < needWidth, insDim(3) = needWidth; end
InputFig.Position = insDim;

% text and slider
x = insDim(3) - dftTX(3) - 10*mgWidth;
y = insDim(4) - dftPB(4) - 4*mgHeight ;

textH{1} =  uicontrol (InputFig          , ...
    'Style'   , 'text'                   , ...
    'Position', [x y  dftTX(3)  dftTX(4)], ...
    'Tag'     , 'tagTextH1'              , ...
    'HorizontalAlignment', 'left'        , ...
    'String'  , 'Choose your baseline');
y = y - 3.8*dftPB(4);

textLB{1} =  uicontrol (InputFig          , ...
    'Style'   , 'listbox'                   , ...
    'Position', [x y  dftTX(3)  4*dftTX(4)], ...
    'Tag'     , 'tagLstStr'              , ...
    'Callback', @changeE                  , ...
    'Value'   , 2, ...
    'String'  , listString);
y = y - 1.5*dftPB(4);

textH{2} =  uicontrol (InputFig          , ...
    'Style'   , 'text'                   , ...
    'Position', [x y  dftTX(3)  dftTX(4)], ...
    'Tag'     , 'tagTextH2'              , ...
    'HorizontalAlignment', 'left'        , ...
    'String'  , 'Lambda (10E2 - 10E9): ');
y = y - 0.6*dftTX(4);

editE{1} =  uicontrol (InputFig          , ...
    'Style'   , 'edit'                   , ...
    'Position', [x y  dftPB(3)  dftPB(4)], ...
    'Tag'     , 'tagEditE1'              , ...
    'Callback', @changeE                  , ...
    'String'  , '10E4');
y = y - 1.5*dftPB(4);

textH{3} =  uicontrol (InputFig          , ...
    'Style'   , 'text'                   , ...
    'Position', [x y  dftTX(3)  dftTX(4)], ...
    'Tag'     , 'tagTextH3'              , ...
    'HorizontalAlignment', 'left'        , ...
    'String'  , 'Alpha (0 - 0.5): ');
y = y - 0.6*dftTX(4);

editE{2} =  uicontrol (InputFig          , ...
    'Style'   , 'edit'                   , ...
    'Position', [x y  dftPB(3)  dftPB(4)], ...
    'Tag'     , 'tagEditE2'              , ...
    'Callback', @changeE                  , ...
    'String'  , '10E-2');
y = y - 1.5*dftPB(4);

editPB{1} =  uicontrol (InputFig          , ...
    'Style'   , 'pushbutton'                   , ...
    'Position', [x y  dftTX(3)  dftPB(4)], ...
    'Tag'     , 'tagEditPB1'              , ...
    'Callback', @pressB                  , ...
    'String'  , 'Change profile');
y = y - 1.2*dftPB(4);

editPB{2} =  uicontrol (InputFig          , ...
    'Style'   , 'pushbutton'                   , ...
    'Position', [x y  dftPB(3)  dftPB(4)], ...
    'Tag'     , 'tagEditPB2'              , ...
    'Callback', @pressB                  , ...
    'String'  , 'Finished');
y = y - 1.2*dftPB(4);

editPB{3} =  uicontrol (InputFig          , ...
    'Style'   , 'pushbutton'                   , ...
    'Position', [x y  dftPB(3)  dftPB(4)], ...
    'Tag'     , 'tagEditPB3'              , ...
    'Callback', @pressB                  , ...
    'String'  , 'Abord');
y = y - 11.3*dftPB(4);

wRem = insDim(3) - max(dftPB(3), dftTX(3));
hRem = insDim(4);

% Place axes
AxisH{1} = axes(InputFig                      , ...
    'Units'        , 'pixels'                , ...
    'Tag'          , 'tagAxisH1'              , ...
    'OuterPosition'     , [0 0 wRem 1/2*hRem]); % main axis
posAx1 = AxisH{1}.OuterPosition;

AxisH{2} = axes(InputFig                      , ...
    'Units'        , 'pixels'                , ...
    'Tag'          , 'tagAxisH2'              , ...
    'OuterPosition'     , [0 1/2*hRem wRem 1/2*hRem]); % profile
posAx2 = AxisH{2}.OuterPosition;


% 1.3. activate the figure and resize function

InputFig.Visible        = 'on';
InputFig.SizeChangedFcn =  @scf;

% 2. creating random trace and plotting
trc    = randomTrc( AxisX, AxisZ, prof2cor, XLim);


if ishghandle(InputFig)
    uiwait(InputFig);
end

if ~ishghandle(InputFig)
    disp('done'); % if no handles and answer 0 either esc or stop the fuck
end

    function method = plotTrc(trc, saveM)
        % 1. First method name and parameter
        if saveM
            newMetd = lastLS{textLB{1}.Value};
            myAT    = AnalyzeThis(trc, ...
                'BaseMethod'  , newMetd, ...
                'SmoothMethod', 'None' , ...
                'PeakPicking' , 'None' , ...
                'Func4Deconv' , 'None');
        else
            newMetd = listString{textLB{1}.Value};
            pm1     = editE{1}.String;
            pm2     = editE{2}.String;
            myAT    = AnalyzeThis(trc, ...
                'BaseMethod'  , [newMetd, ':', pm1, ':', pm2],...
                'SmoothMethod', 'None',                       ...
                'PeakPicking' , 'None',                       ...
                'Func4Deconv' , 'None');
            
        end
        
        method = myAT.BaseMethod;
        lastLS{textLB{1}.Value} = method;
        
        
        plot(AxisH{2}, trc.Data(:,1), trc.Data(:,2), 'k')
        hold on
        plot(AxisH{2}, trc.Data(:,1), myAT.Baseline.vals, 'r')
        hold off
        title('Original profile')
        xlabel([trc.InfoTrc.AxisX.Label, ' / ', trc.InfoTrc.AxisX.Unit]);
        ylabel([trc.InfoTrc.AxisY.Label, ' / ', trc.InfoTrc.AxisY.Unit]);
        
        Yc          = round(trc.Data(:,2)- myAT.Baseline.vals);
        Yc(Yc < 0) = 0;
        plot(AxisH{1}, trc.Data(:,1), Yc, 'k')
        title('Corrected profile')
        xlabel([trc.InfoTrc.AxisX.Label, ' / ', trc.InfoTrc.AxisX.Unit]);
        ylabel([trc.InfoTrc.AxisY.Label, ' / ', trc.InfoTrc.AxisY.Unit]);
        
        drawnow()
        
        DM     = strsplit(method, ':');
        switch DM{1}
            
            case 'PolyFit'
                textLB{1}.Value  = 1;
                editE{1}.Visible = 'on';
                editE{1}.String   = DM{2};
                editE{2}.Visible = 'off';
                textH{2}.Visible = 'on';
                textH{2}.String  = 'Order (0- cst, 1- linear,...):';
                textH{3}.Visible = 'off';
                
            case 'ArPLS'
                textLB{1}.Value  = 2;
                editE{1}.Visible = 'on';
                editE{1}.String   = DM{2};
                editE{2}.Visible = 'on';
                editE{2}.String   = DM{3};
                textH{2}.Visible = 'on';
                textH{2}.String  = 'Lambda (10E2 - 10E9): ';
                textH{3}.Visible = 'on';
                textH{3}.String  = 'Alpha (0 - 0.5): ';
                
            case 'ArPLS2'
                textLB{1}.Value  = 3;
                editE{1}.Visible = 'on';
                editE{1}.String   = DM{2};
                editE{2}.Visible = 'off';
                textH{2}.Visible = 'on';
                textH{2}.String  = 'Lambda (10E2 - 10E9): ';
                textH{3}.Visible = 'off';
        end
        
    end

    function changeE(source,~)
        if strcmpi(source.Tag, 'tagLstStr')
            method = plotTrc(trc, true);
        else
            method = plotTrc(trc, false);
        end
    end

    function trc = randomTrc( AxisX, AxisZ, prof2cor, XLim, Id)
        if nargin == 4
            prf      = AxisX.Data;
            prf(:,2) = prof2cor(int16(rand*size(prof2cor,1)), :);
        else
            prf      = AxisX.Data;
            prf(:,2) = prof2cor(Id, :);
        end
        
        ind2keep = prf(:,1) >= XLim(1) & prf(:,1) <= XLim(2);
        infoTrc.Title = '';
        
        infoTrc.FT     = '';
        infoTrc.TT     = 'PRF';
        infoTrc.AxisX  = Axis(AxisX.InfoAxis);
        infoTrc.AxisY  = Axis(AxisZ.InfoAxis);
        infoTrc.Loc    = 'inTrace';
        infoTrc.AdiPrm = {};
        infoTrc.P2Fin  = '';
        trc = Trace(infoTrc, prf(ind2keep,:));
    end

    function scf(src,callbackdata)
        fig = gcbo;
        insDim = fig.Position;
        if insDim(4) < needHeight, insDim(4) = needHeight; end
        if insDim(3) < needWidth, insDim(3) = needWidth; end
        fig.Position = insDim;
        
        % text and slider
        x = insDim(3) - dftTX(3) - 10*mgWidth;
        y = insDim(4) - dftPB(4) - 4*mgHeight ;
        
        u = findobj(gcbo,'Tag','tagTextH1');
        u.Position = [x y  dftTX(3)  dftTX(4)];
        y = y - 3.8*dftPB(4);
        
        u = findobj(gcbo,'Tag','tagLstStr');
        u.Position =  [x y  dftTX(3)  4*dftTX(4)];
        y = y - 1.5*dftPB(4);
        
        u = findobj(gcbo,'Tag', 'tagTextH2');
        u.Position =  [x y  dftTX(3)  dftTX(4)];
        y = y - 0.6*dftTX(4);
        
        u = findobj(gcbo,'Tag','tagEditE1');
        u.Position = [x y  dftPB(3)  dftPB(4)];
        y = y - 1.5*dftPB(4);
        
        u = findobj(gcbo,'Tag','tagTextH3');
        u.Position = [x y  dftTX(3)  dftTX(4)];
        y = y - 0.6*dftTX(4);
        
        u = findobj(gcbo,'Tag','tagEditE2');
        u.Position = [x y  dftPB(3)  dftPB(4)];
        y = y - 1.5*dftPB(4);
        
        u = findobj(gcbo,'Tag','tagEditPB1');
        u.Position = [x y  dftTX(3)  dftPB(4)];
        y = y - 1.3*dftPB(4);
        
        u = findobj(gcbo,'Tag','tagEditPB2');
        u.Position = [x y  dftPB(3)  dftPB(4)];
        y = y - 1.3*dftPB(4);
        
        u = findobj(gcbo,'Tag','tagEditPB3');
        u.Position = [x y  dftPB(3)  dftPB(4)];
        y = y - 1.3*dftPB(4);
        
        wRem = insDim(3) - max(dftPB(3), dftTX(3));
        hRem = insDim(4);
        
        u = findobj(gcbo,'Tag','tagAxisH1');
        if ~ isempty(u)
            u.OuterPosition = [0 0 wRem 1/2*hRem]; % main axis
            posAx1 = u.OuterPosition;
        end
        
        u = findobj(gcbo,'Tag','tagAxisH2');
        if ~ isempty(u)
            u.OuterPosition = [0 1/2*hRem wRem 1/2*hRem];
            posAx2 = u.OuterPosition;
        end
        
        u = findobj(gcbo,'Tag','tagAxisH3');
        if~isempty(u)
            u.OuterPosition = [0 3*dftPB(4) 1/4*wRem 3/4*hRem];
            posAx3 = u.OuterPosition;
        end
    end

    function pressB(source, ~)
        switch source.Tag
            case 'tagEditPB1'
                trc    = randomTrc( AxisX, AxisZ, prof2cor, XLim);
                plotTrc(trc, false);
                
            case 'tagEditPB2'
                ExitFlag = 0;
                close(InputFig)
                
            case 'tagEditPB3'
                ExitFlag = -1;
                close(InputFig)
        end
    end
end

