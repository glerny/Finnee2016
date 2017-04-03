%% DESCRIPTION
% GUI for @PeakList
%% Copyright
% BSD 3-Clause License
% Copyright 2016-2017 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gui4ClustPlot(aPL, CluLvl, HACA, FOMPIP, options, BPP, LstPIP)

[~, name, ~] = fileparts(aPL.Path2PkL);
strName = sprintf('Clusters plot of %s at CL = %.2f', name, CluLvl);
assignin('base', 'FOMPIP', FOMPIP)
assignin('base', 'HACA', HACA)
BPEH = subplot(5, 1, 1);
h1 = plot(BPP.Data(:,1), BPP.Data(:,2), 'r');
h1.Tag = 'BPP';
title(BPP.Title);
xlabel([BPP.AxisX.Label, ' / ', BPP.AxisX.Unit]);
ylabel([BPP.AxisY.Label, ' / ', BPP.AxisY.Unit]);

ClPH   = subplot(5,1, 2:5); % CLUSTERS PLOT
ThI    = options.IntThreshold;
ThP    = options.minPIPS;
Id2Plt = FOMPIP(:,2) >= ThP & FOMPIP(:,7) >= ThI;

codeSize = int32(FOMPIP(Id2Plt,7)/max(FOMPIP(Id2Plt,7))*(500) + 5);
h2 = scatter(FOMPIP(Id2Plt,4), FOMPIP(Id2Plt,5), codeSize, 'k');
h2.Tag = 'ClstPt';
title('CLUSTERS PLOT - uses the Data cursor to gain more information');
xlabel([aPL.AxisX.Label, ' / ', aPL.AxisX.Unit]);
ylabel([aPL.AxisZ.Label, ' / ', aPL.AxisZ.Unit]);

linkaxes([BPEH,ClPH],'x')

FigInput = gcf;
FigInput.Name = sprintf('Clusters plot of %s at CL = %.2f', name, CluLvl);
dcm_obj = datacursormode(FigInput);
set(dcm_obj,'UpdateFcn',@myupdatefcn);
figure(FigInput)

    function txt = myupdatefcn(~,event_obj)
        % New data cursor update function
        
        tgt = get(event_obj,'Target');
        pos = get(event_obj,'Position');
        tgt.Tag
        switch tgt.Tag;
            case 'BPP'
                xString = [aPL.AxisX.Label, ' = ', num2str(pos(1), aPL.AxisX.fo), ...
                    ' ', aPL.AxisX.Unit];
                yString = [aPL.AxisZ.Label, ' = ', num2str(pos(2), aPL.AxisZ.fo), ...
                    ' ', aPL.AxisZ.Unit];
                txt = {xString, yString};
            case 'ClstPt'
                Id   = find(FOMPIP(:,4) == pos(1) & FOMPIP(:,5) == pos(2));
                cFOM = FOMPIP(Id, :); Id2HACA = cFOM(1);
                zrtL = ['Cluster # ', num2str(Id2HACA)];
                fstl = ['Mean(M1): ',  num2str(cFOM(4), aPL.AxisX.fo), ...
                    ' ', aPL.AxisX.Unit];
                sndl = ['sum(area): ',  num2str(cFOM(6), '%.3e'), ...
                    ' a.u.'];
                thrl = ['Acc. Mass @ BP: ',  num2str(cFOM(5), aPL.AxisY.fo)];
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
        xlabel([aPL.AxisX.Label, ' / ', aPL.AxisX.Unit]);
        ylabel([aPL.AxisZ.Label, ' / ', aPL.AxisZ.Unit]);

        hold on
        AllFOM = [];
        for ii = 1:length(HACA{IdIn}.LstPIP)
            cPIP = LstPIP{HACA{IdIn}.LstPIP(ii)};
            plot(AxisH1, cPIP.x, cPIP.y, 'k');
            AllFOM = [AllFOM; cPIP.FOM];
        end
        hold off
        assignin('base', 'AllFOM', AllFOM)
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

