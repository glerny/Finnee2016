function newPeakList = Gui4PeakListReview(myPeakList, name)

newPeakList = myPeakList;
options.CI = 1.96;
CI = options.CI;
myGuiFig = figure(                ...
    'Name'         , name        ,...
    'Units'        , 'normalized',...
    'Outerposition',    [0 0 1 1],...
    'UserData'     , myPeakList,...
    'Resize'       , 'off');
DefPos = myGuiFig.Position;

% Place axes
AxisH{1} = axes(myGuiFig                      , ...
    'Units'        , 'normalized'                , ...
    'Tag'          , 'tagAxisH1'              , ...
    'OuterPosition'     , [DefPos(1) DefPos(2) DefPos(3)*2/3-DefPos(1) DefPos(4)*1/3-DefPos(2)]); % main axis

AxisH{2} = axes(myGuiFig                      , ...
    'Units'        , 'normalized'                , ...
    'Tag'          , 'tagAxisH2'              , ...
    'OuterPosition'     , [DefPos(1) +DefPos(4)*1/3 + DefPos(2) DefPos(3)*2/3-DefPos(1) DefPos(4)*1/3-DefPos(2)]); % main axis

AxisH{3} = axes(myGuiFig                      , ...
    'Units'        , 'normalized'                , ...
    'Tag'          , 'tagAxisH2'              , ...
    'OuterPosition'     , [DefPos(1) +DefPos(4)*2/3 + DefPos(2) DefPos(3)*2/3-DefPos(1) DefPos(4)*1/3-DefPos(2)]); % main axis

testPb    = uicontrol (myGuiFig          , ...
    'Style'   , 'pushbutton'                   , ...
    'Units'   , 'normalized',...
    'Visible' , 'off',...
    'String'   , 'X');
DSC = testPb.Extent;

x = DefPos(1)+DefPos(3)*2/3;
y = DefPos(2)+3*DSC(4);

AxisH{4} = axes(myGuiFig                      , ...
    'Units'        , 'normalized'                , ...
    'Tag'          , 'tagAxisH2'              , ...
    'OuterPosition', [x y DefPos(3)*1/3-DefPos(1) DefPos(4)*1/2-DefPos(2)]); % main axis

pos = AxisH{4}.Position;
x = pos(1) + DSC(3);
y = y - 2*DSC(4);

PushH{1} = uicontrol (myGuiFig , ...
    'Style'   , 'pushbutton'   , ...
    'Units'   , 'normalized'   , ...
    'Tag'     , 'PB1'          , ...
    'Callback' , @pushMe        , ...
    'String'  , '<'            , ...
    'Position', [x y 3*DSC(3) 1.2*DSC(4)]);

x = x + 4*DSC(3);
PushH{2} = uicontrol (myGuiFig , ...
    'Style'   , 'pushbutton'   , ...
    'Units'   , 'normalized'   , ...
    'Tag'     , 'PB2'          , ...
    'Callback' , @pushMe        , ...
    'String'  , '>'            , ...
    'Position', [x y 3*DSC(3) 1.2*DSC(4)]);

x = x + 4*DSC(3);
PushH{3} = uicontrol (myGuiFig , ...
    'Style'   , 'pushbutton'   , ...
    'Units'   , 'normalized'   , ...
    'Tag'     , 'PB3'          , ...
    'Callback' , @pushMe        , ...
    'String'  , 'X'            , ...
    'Position', [x y 3*DSC(3) 1.2*DSC(4)]);

x = x + 4*DSC(3);
PushH{4} = uicontrol (myGuiFig , ...
    'Style'   , 'pushbutton'   , ...
    'Units'   , 'normalized'   , ...
    'Tag'     , 'PB4'          , ...
    'Callback' , @pushMe        , ...
    'String'  , 'reprocess'            , ...
    'Position', [x y 12*DSC(3) 1.2*DSC(4)]);

x =  pos(1) + DSC(3);
y = 1 - 3*DSC(4);
PushH{5} = uicontrol (myGuiFig , ...
    'Style'   , 'pushbutton'   , ...
    'Units'   , 'normalized'   , ...
    'Tag'     , 'PB5'          , ...
    'Callback' , @pushMe        , ...
    'String'  , 'Remove filtered'            , ...
    'Position', [x y 20*DSC(3) 1.2*DSC(4)]);

x =  x + 23*DSC(3);
y = 1 - 3*DSC(4);
PushH{6} = uicontrol (myGuiFig , ...
    'Style'   , 'pushbutton'   , ...
    'Units'   , 'normalized'   , ...
    'Tag'     , 'PB6'          , ...
    'Callback' , @pushMe        , ...
    'String'  , 'Finished!'      , ...
    'Position', [x y 20*DSC(3) 1.2*DSC(4)]);

string4{1} = 'Total number of PIP: ';
string4{2} = 'PIP to be checked  : ';
string4{3} = 'PIP being checked  : ';

x =  pos(1) + DSC(3);
y = 1 - 5*DSC(4);
TextH{1} = uicontrol (myGuiFig , ...
    'Style'   , 'text'         , ...
    'Units'   , 'normalized'   , ...
    'Tag'     , 'Tx1'          , ...
    'String'  , string4{1}     , ...
    'HorizontalAlignment', 'left',...
    'Position', [x y 25*DSC(3) 1.2*DSC(4)]);


x =  pos(1) + DSC(3);
y = 1 - 6.5*DSC(4);
TextH{2} = uicontrol (myGuiFig , ...
    'Style'   , 'text'   , ...
    'Units'   , 'normalized'   , ...
    'Tag'     , 'Tx2'          , ...
    'String'  , string4{2}     , ...
    'HorizontalAlignment', 'left',...
    'Position', [x y 25*DSC(3) 1.2*DSC(4)]);

x =  pos(1) + DSC(3);
y = 1 - 8*DSC(4);
TextH{3} = uicontrol (myGuiFig , ...
    'Style'   , 'text'   , ...
    'Units'   , 'normalized'   , ...
    'Tag'     , 'Tx3'          , ...
    'String'  , string4{3}     , ...
    'HorizontalAlignment', 'left',...
    'Position', [x y 25*DSC(3) 1.2*DSC(4)]);

IdC = 1;
[X, Y, Id2Rem] = filter(myPeakList);
tgi = find(Id2Rem);

fillAxis4(X(tgi(IdC)), Y(tgi(IdC)), myPeakList.LstPIP{1}{tgi(IdC)});

if ishghandle(myGuiFig)
    uiwait(myGuiFig);
end

if ~ishghandle(myGuiFig)
    disp('done'); % if no handles and answer 0 either esc or stop the fuck
end


    function [X, Y, Id2Rem] = filter(myPeakList, Id2Rem)
        
        
        rpl    = length(myPeakList.FOM);
        dt4PLR = [];
        for ii = 1:rpl
            dt4PLR(:,:,ii) = myPeakList.FOM{rpl}.Data;
        end
        X      = mean(dt4PLR(:,5,:), 3);
        Y      = mean(dt4PLR(:,6,:), 3);
        p      = polyfitweighted(X.^2,Y,3);
        Y_     = polyval(p, X.^2);
        
        if nargin == 1
            %[X, Ix] = sort(X); Y = Y(Ix);
            Id2Rem  = abs(Y-Y_) > CI*std(Y-Y_);
        end
        
        scatter(AxisH{3}, X(~Id2Rem), Y(~Id2Rem), 'k')
        hold(AxisH{3}, 'on')
        scatter(AxisH{3}, X(Id2Rem), Y(Id2Rem), 'r')
        [~, Ix]  = sort(X);
        plot(AxisH{3}, X(Ix), Y_(Ix), 'g');
        plot(AxisH{3}, X(Ix), Y_(Ix)+CI*std(Y-Y_), 'r');
        plot(AxisH{3}, X(Ix), Y_(Ix)-CI*std(Y-Y_), 'r');
        hold(AxisH{3}, 'off')
        
        for jj = 1:rpl
            BPP{jj}(:,1) = myPeakList.BPP{jj}.Data(:,1);
            BPP{jj}(:,4) = 0;
            TIP{jj}(:,1) = myPeakList.BPP{jj}.Data(:,1);
            TIP{jj}(:,4) = 0;
            for ii = 1:length(myPeakList.LstPIP{jj})
                cPIP = myPeakList.LstPIP{jj}{ii};
                IdS  = cPIP.IdS;
                dt    = cPIP.y;
                IdE  = IdS + length(dt) - 1;
                try
                    BPP{jj}(IdS:IdE, 2) = max(BPP{jj}(IdS:IdE, 2), dt);
                catch
                    disp('WTF')
                end
                TIP{jj}(IdS:IdE, 2) = TIP{jj}(IdS:IdE, 2) + dt;
                
                if Id2Rem(ii)
                    BPP{jj}(IdS:IdE, 4) = max(BPP{jj}(IdS:IdE, 4), dt);
                    TIP{jj}(IdS:IdE, 4) = TIP{jj}(IdS:IdE, 4) + dt;
                else
                    BPP{jj}(IdS:IdE, 3) = max(BPP{jj}(IdS:IdE, 3), dt);
                    TIP{jj}(IdS:IdE, 3) = TIP{jj}(IdS:IdE, 3) + dt;
                end
            end
        end
        UD.BPP = BPP;
        UD.TIP = TIP;
        AxisH{3}.UserData = UD;
        
        plot(AxisH{2}, BPP{1}(:,1), BPP{1}(:,2), 'k')
        hold(AxisH{2}, 'on')
        plot(AxisH{2}, BPP{1}(:,1), BPP{1}(:,3), 'r')
        hold(AxisH{2}, 'off')
        
        plot(AxisH{1}, BPP{1}(:,1), BPP{1}(:,4), 'k')
        
        TextH{1}.String = [string4{1}, num2str(size(myPeakList.FOM{1}.Data, 1))];
        TextH{2}.String = [string4{2}, num2str(sum(Id2Rem))];
        TextH{3}.String = [string4{3}, num2str(IdC)];
        
    end

    function fillAxis4(x, y, PIP)
        hold(AxisH{3}, 'on')
        scatter(AxisH{3}, x, y, 'g')
        hold(AxisH{3}, 'off')
        
        sizePrm = (PIP.Data(:,2)/max(PIP.Data(:,2))+0.05)*100;
        yyaxis(AxisH{4}, 'left')
        plot(AxisH{4}, PIP.x, PIP.y);
        yyaxis(AxisH{4}, 'right')
        scatter(AxisH{4}, ...
            PIP.x(PIP.Data(:,3)- min(PIP.Data(:,3)) + 1), PIP.Data(:,1), sizePrm);
        
        title('Intensities and accurate masses')
        xlabel([PIP.AxisX.Label, ' / ', PIP.AxisX.Unit]);
        AxisH{4}.YAxis(1).Label.String = ...
            [PIP.AxisZ.Label, ' / ', PIP.AxisZ.Unit];
        AxisH{4}.YAxis(2).Label.String = ...
            [PIP.AxisY.Label, ' / ', PIP.AxisY.Unit];
    end

    function filtPeakList
        rpl    = length(myPeakList.FOM);
        newPeakList = myPeakList;
        UD = AxisH{3}.UserData;
        for ii = 1:rpl
            InfoBPP{ii} = myPeakList.BPP{ii}.InfoTrc;
            newPeakList.BPP{ii} = Trace(InfoBPP{ii},  UD.BPP{ii}(:, [1 3]));
            InfoTIP{ii} = myPeakList.TIP{ii}.InfoTrc;
            newPeakList.TIP{ii} = Trace(InfoTIP{ii},  UD.TIP{ii}(:, [1 3]));
            newPeakList.FOM{ii}.Data(Id2Rem, :) = [];
            newPeakList.FOM{ii}.Data(:,1) = 1:size(newPeakList.FOM{ii}.Data(:, 1), 1);
            newPeakList.LstPIP{ii}(Id2Rem) = [];
        end
        
        myPeakList = newPeakList;
    end


    function pushMe(source, obj)
        switch source.Tag
            case 'PB1'
                
                if IdC == 1
                    IdC = length(tgi);
                else
                    IdC = IdC - 1;
                end
                
                hold(AxisH{3}, 'on')
                scatter(AxisH{3}, X(tgi(IdC)), Y(tgi(IdC)), 'r')
                hold(AxisH{3}, 'off')
                
                fillAxis4(X(tgi(IdC)), Y(tgi(IdC)), myPeakList.LstPIP{1}{tgi(IdC)});
                TextH{3}.String = [string4{3}, num2str(IdC)];
                
            case 'PB2'
                hold(AxisH{3}, 'on')
                if IdC > length(tgi)
                    IdC = 1;
                end
                
                scatter(AxisH{3}, X(tgi(IdC)), Y(tgi(IdC)), 'r')
                hold(AxisH{3}, 'off')
                
                if IdC == length(tgi)
                    IdC = 1;
                else
                    IdC = IdC + 1;
                end
                
                fillAxis4(X(tgi(IdC)), Y(tgi(IdC)), myPeakList.LstPIP{1}{tgi(IdC)});
                TextH{3}.String = [string4{3}, num2str(IdC)];
                
            case 'PB3'
                
                hold(AxisH{3}, 'on')
                scatter(AxisH{3}, X(tgi(IdC)), Y(tgi(IdC)), 'k')
                hold(AxisH{3}, 'off')
                
                Id2Rem(tgi(IdC)) = false;
                tgi = find(Id2Rem);
                
                if IdC == length(tgi)
                    IdC = 1;
                else
                    IdC = IdC + 1;
                end
                
                fillAxis4(X(tgi(IdC)), Y(tgi(IdC)), myPeakList.LstPIP{1}{tgi(IdC)});
                filter(myPeakList, Id2Rem)
                
            case 'PB4'
                
                rpl    = length(myPeakList.FOM);
                hold(AxisH{3}, 'on')
                scatter(AxisH{3}, X(tgi(IdC)), Y(tgi(IdC)), 'w')
                hold(AxisH{3}, 'off')
                for ii = 1:rpl
                    cPIP = myPeakList.LstPIP{ii}{tgi(IdC)};
                    PIPs = cPIP.manualBaseline(['replicate nbr: ', num2str(ii)]);
                    myPeakList.LstPIP{ii}{tgi(IdC)} = cPIP;
                    myPeakList.FOM{ii}.Data(tgi(IdC), :) = ...
                        [myPeakList.FOM{ii}.Data(tgi(IdC), 1), cPIP.FOM];
                end
                
                Id2Rem(tgi(IdC)) = [];
                myPeakList.LstPIP{1}(tgi(IdC)) = [];
                myPeakList.FOM{1}.Data(tgi(IdC), :) = [];
                
                for ii = 1:length(PIPs)
                    k = length(Id2Rem) + 1;
                    Id2Rem(k) = false;
                    myPeakList.LstPIP{1}{k} = PIPs{ii};
                    myPeakList.FOM{1}.Data(k, :) = [k, PIPs{ii}.FOM];
                end
                
                tgi = find(Id2Rem);
                
                if IdC > length(tgi)
                    IdC = 1;
                else
                    IdC = IdC + 1;
                end
                
                fillAxis4(X(tgi(IdC)), Y(tgi(IdC)), myPeakList.LstPIP{1}{tgi(IdC)});
                filter(myPeakList, Id2Rem)
                
                
            case 'PB5'
                filtPeakList
                [X, Y, Id2Rem] = filter(myPeakList);
                tgi = find(Id2Rem);
                if IdC > length(tgi)
                    IdC = 1;
                end
                if isempty(tgi)
                else
                    fillAxis4(X(tgi(IdC)), Y(tgi(IdC)), myPeakList.LstPIP{1}{tgi(IdC)});
                end
                
                 case 'PB6'
                     close(myGuiFig)
        end
        
    end

end