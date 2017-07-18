%% DESCRIPTION
% GUI for @PeakList
%% Copyright
% BSD 3-Clause License
% Copyright 2016-2017 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gui4ClustPlot(aPL, CluLvl, HACA, FOMPIP, options, BPP)

BPEH = subplot(5, 1, 1);
h1 = plot(BPP{1}.Data(:,1), BPP{1}.Data(:,2), 'r');
h1.Tag = 'BPP';
title(BPP{1}.Title);
xlabel([BPP{1}.AxisX.Label, ' / ', BPP{1}.AxisX.Unit]);
ylabel([BPP{1}.AxisY.Label, ' / ', BPP{1}.AxisY.Unit]);

ClPH   = subplot(5,1, 2:5); % CLUSTERS PLOT
ThI    = options.IntThreshold;
ThP    = options.minPIPS;
Id2Plt = FOMPIP(:,2) >= ThP & FOMPIP(:,7) >= ThI;

codeSize = int32(sqrt(FOMPIP(Id2Plt,7))/max(sqrt(FOMPIP(Id2Plt,7)))*(500) + 5);
h2 = scatter(FOMPIP(Id2Plt,4), FOMPIP(Id2Plt,5), codeSize, 'k');
h2.Tag = 'ClstPt';
title('CLUSTERS PLOT - uses the Data cursor to gain more information');
xlabel([aPL.AxisX{1}.Label, ' / ', aPL.AxisX{1}.Unit]);
ylabel([aPL.AxisZ{1}.Label, ' / ', aPL.AxisZ{1}.Unit]);

linkaxes([BPEH,ClPH],'x')

FigInput = gcf;
FigInput.Name = sprintf('Clusters plot at CL = %.2f', CluLvl);
dcm_obj = datacursormode(FigInput);
set(dcm_obj,'UpdateFcn',@myupdatefcn);
figure(FigInput)

    function txt = myupdatefcn(~,event_obj)
        % New data cursor update function
        
        tgt = get(event_obj,'Target');
        pos = get(event_obj,'Position');
        switch tgt.Tag;
            case 'BPP'
                xString = [aPL.AxisX{1}.Label, ' = ', num2str(pos(1), aPL.AxisX{1}.fo), ...
                    ' ', aPL.AxisX{1}.Unit];
                yString = [aPL.AxisZ{1}.Label, ' = ', num2str(pos(2), aPL.AxisZ{1}.fo), ...
                    ' ', aPL.AxisZ{1}.Unit];
                txt = {xString, yString};
            case 'ClstPt'
                Id   = find(FOMPIP(:,4) == pos(1) & FOMPIP(:,5) == pos(2));
                cFOM = FOMPIP(Id, :); Id2HACA = cFOM(1);
                zrtL = ['Cluster # ', num2str(Id2HACA)];
                fstl = ['Mean(M1): ',  num2str(cFOM(4), aPL.AxisX{1}.fo), ...
                    ' ', aPL.AxisX{1}.Unit];
                sndl = ['max(area): ',  num2str(cFOM(6), '%.3e'), ...
                    ' a.u.'];
                thrl = ['Acc. Mass @ BP: ',  num2str(cFOM(5), aPL.AxisY{1}.fo)];
                txt = {zrtL, fstl, sndl, thrl};
                plotCluster(Id2HACA)
        end
        
    end

    function plotCluster(IdIn)
        
        InputFig = figure(                          ...
            'Visible'          , 'on'             , ...
            'Name'             , ['Cluster #', num2str(IdIn)], ...
            'MenuBar'          , 'none'           , ...
            'Units'            , 'normalized'     , ...
            'WindowStyle'      , 'normal'         , ...
            'Resize'           , 'on');
        
        AxisH1 = axes(InputFig ,...
            'Units'            , 'normalized'     , ...
            'OuterPosition'    , [0 0 2/3 1]);
        
        
        title('Superposition of PIPs');
        xlabel([aPL.AxisX{1}.Label, ' / ', aPL.AxisX{1}.Unit]);
        ylabel([aPL.AxisZ{1}.Label, ' / ', aPL.AxisZ{1}.Unit]);

        hold on
        AllFOM = [];
        for ii = 1:length(HACA{IdIn}.LstPIP)
            cPIP = aPL.LstPIP{1}{HACA{IdIn}.LstPIP(ii)};
            plot(AxisH1, cPIP.x, cPIP.y, 'k');
            AllFOM = [AllFOM; cPIP.FOM];
        end
        hold off
        
        AllFOM = sortrows(AllFOM, -1);
        AllFOM(:,10) = AllFOM(:,1)/max(AllFOM(:,1))*100;
        
        String4Edit{1} = '  LIST OF PURE ION PROFILES';
        String4Edit{2} = '____________________________';
        String4Edit{3} = '';
        String4Edit{4} = sprintf('m/z \t Area(a.u.) \t Int(%%)');
        String4Edit{5} = '';
        
        
        for ii =1:length(AllFOM(:,1))
            String4Edit{end+1} = sprintf('%0.5f\t%0.f\t%0.1f', ...
                AllFOM(ii,9), AllFOM(ii, 3), AllFOM(ii, 10));
        end
        
        
        EditH1 = uicontrol(InputFig ,...
            'Style'              , 'Edit'      , ...
            'Visible'            , 'on'        , ...
            'Units'              , 'normalized', ...
            'Max'                , 2           , ...
            'String'             , String4Edit ,...
            'HorizontalAlignment', 'left', ...
            'Position'           , [2/3 0 1/3 1]);
        
    end
end

