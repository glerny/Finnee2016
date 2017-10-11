%% DESCRIPTION
%
%% Copyright
% BSD 3-Clause License
% Copyright 2016-2017 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dataOut = peakMatching(obj, parameters, varargin)

options.verified  = false;
options.adjusted  = false;

options.CI4PM     = 2.807;
options.meth4Ans  = 'ttest';
options.prm1      = 0.05; % alpha
options.prm2      = 'Control'; % Tag for reference
options.prm3      = 4; %minimum non nul values
options.prm4      = 1;
options.prm5      = 'E:\iBET\EBC11072017\Athma\Figures';

infoX = obj.Replicates{1}.PeakList.AxisX{1}.InfoAxis;
foX   = obj.Replicates{1}.PeakList.AxisX{1}.fo;
infoZ = obj.Replicates{1}.PeakList.AxisZ{1}.InfoAxis;
foZ   = obj.Replicates{1}.PeakList.AxisZ{1}.fo;

map4tags = false(length(obj.Replicates), 1);
for ii = 1:length(obj.Replicates)
    if strcmpi(obj.Replicates{ii}.Tag, options.prm2)
        map4tags(ii) = true;
    end
end


% STEP 1. Peak matching
% 1.1. Find peaks common to all dataset using original parameters
% 1.1.1. Load default parameters

CI   = options.CI4PM;
P    = strsplit(parameters, ':');
Dmz  = str2double(P{1});
RSDT = str2double(P{2});
ThIt = str2double(P{3});

% (not needed now but done to gain time)
% ----------------
listOfTags = obj.ListOfTags;
tmap       = [];
% ---------------

% 1.1.2. Build a matrix will all FOM
nrpli = length(obj.Replicates);
allFOM = [];
for ii = 1:nrpli
    FOM = obj.Replicates{ii}.PeakList.FOM{1}.Data;
    FOM(:, end+1) = ii;
    allFOM =  [allFOM; FOM];
    BPP{ii}(:,1) = obj.Replicates{ii}.PeakList.AxisX{1}.Data;
    BPP{ii}(:,4) = 0;
    TIP{ii}(:,1) = obj.Replicates{ii}.PeakList.AxisX{1}.Data;
    TIP{ii}(:,4) = 0;
    tmap(ii) = find(strcmp(listOfTags, obj.Replicates{ii}.Tag));
end

% 1.1.2. Remove FOM of max Intnesity lower than ThIT Important to improve
% accuracy (SHOULD ADD CONTROL)

allFOM = sortrows(allFOM, 10);
selFOM = allFOM(allFOM(:,11) >= ThIt, :);
Id = [0; find(diff(selFOM(:,10)) > Dmz); size(selFOM, 1)];
data4norm = [];
for ii = 1:length(Id)-1
    Cut = selFOM (Id(ii)+1:Id(ii+1), :);
    if size(unique(Cut(:,12)), 1) == nrpli
        Cut = sortrows(Cut, 5);
        Id2 = [0; find(diff(Cut(:,5))./Cut(1:end-1, 5)*100 > RSDT); size(Cut, 1)];
        for jj = 1:length(Id2)-1
            Cut2 = Cut (Id2(jj)+1:Id2(jj+1), :);
            if size(unique(Cut2(:,12)), 1) == nrpli
                if size(Cut2, 1) ==  nrpli
                    data4norm(:,:, end+1) = Cut2;
                end
            end
        end
    end
end

data4norm(:,:, 1) = [];

% mz normalisation
dt4mz = squeeze(data4norm(:, 10,:));
mstMz = mean(dt4mz);
for ii = 1:nrpli
    p{ii,1} = polyfitweighted(dt4mz(ii, :), mstMz-dt4mz(ii,:), 1);
    dt4mzc(ii,:) = dt4mz(ii,:) + polyval(p{ii,1}, dt4mz(ii,:));
end
dt4std  = max(dt4mzc) - min(dt4mzc);
DelMz = mean(dt4std) + CI*std(dt4std);
disp(DelMz)

% Time normalisation
dt4time = squeeze(data4norm(:, 5,:));
mstTime = mean(dt4time);
for ii = 1:nrpli
    p{ii,2} = polyfitweighted(dt4time(ii, :), mstTime-dt4time(ii,:), 1);
    dt4timec(ii,:) = dt4time(ii,:) + polyval(p{ii,2}, dt4time(ii,:));
end
dt4std  =  max(dt4timec) - min(dt4timec);
DelTm = mean(dt4std) + CI*std(dt4std);
disp(DelTm)

% 1.2. Peak matching
% 1.2.1 Normalisation
for ii = 1:nrpli
    Id2norm = allFOM(:,12) == ii;
    allFOM(Id2norm, 10) = allFOM(Id2norm, 10) + polyval(p{ii, 1},  allFOM(Id2norm, 10));
    allFOM(Id2norm, 5)  = allFOM(Id2norm, 5)  + polyval(p{ii, 2},  allFOM(Id2norm, 5));
end

% 1.2.2 Peak matching
allFOM = sortrows(allFOM, 10);
Id = [0; find(diff(allFOM(:,10)) > DelMz); size(allFOM, 1)];
IdX = {};

for ii = 1:length(Id)-1
    Cut = allFOM (Id(ii)+1:Id(ii+1), :);
    Cut = sortrows(Cut, 5);
    
    if size(Cut, 1) > 1
        Id2 = [0; find(diff(Cut(:,5)) > DelTm); size(Cut, 1)];
        
        for jj = 1:length(Id2)-1
            Cut2 = Cut (Id2(jj)+1:Id2(jj+1), :);
            
            if size(Cut2, 1) ==  size(unique(Cut2(:,11)), 1)
                Id2PIP = NaN(1, nrpli);
                Id2Int = zeros(1, nrpli);
                Id2PIP(Cut2(:,12)) = Cut2(:,1);
                Id2Int(Cut2(:,12)) = Cut2(:,2);
                IdX(end+1).Id = Id2PIP;
                IdX(end).Int = Id2Int;
                
            else
                
                while ~isempty(Cut2)
                    Cut3 = [];
                    
                    for ll = 1:nrpli
                        if  sum(Cut2(:,12) == ll) == 1
                            Cut3(end+1, :) = Cut2(Cut2(:,12) == ll,:);
                            Cut2(Cut2(:,12) == ll,:) = [];
                        end
                    end
                    
                    if ~isempty(Cut2)
                        Id2Ctrl = unique(Cut2(:,12));
                        
                        for ll = 1:length(Id2Ctrl)
                            Idd = find(Cut2(:,12) == Id2Ctrl(ll));
                            if isempty(Cut3)
                                try
                                    Cut3(end+1, :) = Cut2(Idd(1),:);
                                catch
                                    disp('whattheheck')
                                end
                                Cut2(Idd(1),:) = [];
                            else
                                
                                dEuc = (mean(Cut3(:,5))/DelTm).^2 + (mean(Cut3(:,10))/DelMz).^2;
                                [~, IdF] = sort(abs((Cut2(Idd,5)/DelTm).^2 + (Cut2(Idd,10)/DelMz).^2 - dEuc));
                                if isempty(IdF)
                                    disp('WTF')
                                end
                                Cut3(end+1, :) = Cut2(Idd(IdF(1)),:);
                                Cut2(Idd(IdF(1)),:) = [];
                            end
                        end
                    end
                    Id2PIP = NaN(1, nrpli);
                    Id2Int = zeros(1, nrpli);
                    Id2PIP(Cut3(:,12)) = Cut3(:,1);
                    Id2Int(Cut3(:,12)) = Cut3(:,2);
                    IdX(end+1).Id = Id2PIP;
                    IdX(end).Int = Id2Int;
                end
                
            end
        end
    else
        Id2PIP = NaN(1, nrpli);
        Id2Int = zeros(1, nrpli);
        Id2PIP(Cut(:,12)) = Cut(:,1);
        Id2Int(Cut(:,12)) = Cut(:,2);
        IdX(end+1).Id = Id2PIP;
        IdX(end).Int = Id2Int;
    end
end

% STEP 2. Data analysis

switch options.meth4Ans
    case 'ttest'
        
        AxisX = obj.Replicates{1}.PeakList.AxisX{1};
        AxisY = obj.Replicates{1}.PeakList.AxisZ{1};
        
        % 1. Organise the data
        matrixOfIntensity = zeros(size(IdX, 2), max(size(obj.Replicates)));
        allFiltData{size(IdX, 2)}       = {};
        for ii = 1:size(IdX, 2)
            cLi = IdX(ii);
            IdR = find(~isnan( cLi.Id));
            sFOM = [];
            for jj = 1:length(IdR)
                IdPi = cLi.Id(IdR(jj));
                cPIP = obj.Replicates{IdR(jj)}.PeakList.LstPIP{1}{IdPi};
                sFOM = [sFOM; [IdR(jj), IdPi, cPIP.FOM]];
            end
            matrixOfIntensity(ii, :) = cLi.Int;
            allFiltData{ii} = sFOM;
        end
        
        assignin('base', 'Markers', allFiltData);
        assignin('base', 'IntMarkers', matrixOfIntensity)
        return
        disp('done')
        
         % 2. Find data that will be used for references
         %2.1. Find data with no zeros
         IdFull = find(sum(matrixOfIntensity == 0, 2) == 0);
         
        % 2.2 Normalised and do ttest2
        x = matrixOfIntensity(IdFull, map4tags);
        for ii = 1:size(x,2)
            xnorm (:,ii) =  x(:,ii)/sum(x(:,ii));
        end
        y = matrixOfIntensity(IdFull, ~map4tags);
        for ii = 1:size(y,2)
            ynorm (:,ii) =  y(:,ii)/sum(y(:,ii));
        end
        h = ttest2(xnorm,ynorm,'Vartype','unequal', 'Dim', 2,  'Alpha',options.prm1);
        
        % 2.3. Calculate data for normalisation
        Int4Refs = sum(matrixOfIntensity(IdFull( h == 0),:), 1);
        assignin('base', 'Data4Reference', allFiltData(IdFull( h == 0)));
        
        % 3. Analysis
        % 3.1. Find entries with more or equal options.prm3 non nuls
        % values
        Id2Anal = find(sum(matrixOfIntensity ~= 0, 2) >= options.prm3);
        
         % 3.2 Normalised and do ttest2
        XY = matrixOfIntensity(Id2Anal, :);
        for ii = 1:size(XY,2)
            XYnorm (:,ii) =  XY(:,ii)/Int4Refs(:,ii);
        end
        assignin('base', 'XYnorm', XYnorm)
        h = ttest2(XYnorm(:, map4tags), XYnorm(:, ~map4tags),'Vartype','unequal', 'Dim', 2,  'Alpha',options.prm1);
        assignin('base', 'Markers', allFiltData(Id2Anal( h == 1)));
        assignin('base', 'IntMarkers', XYnorm(Id2Anal( h == 1), :))
        
        return

        % 4. Plot results
        Id2Markers = Id2Anal(h == 1);
        dt2plot    = [];
        
        for ii = 1:length(Id2Markers)
            cL = allFiltData{Id2Markers(ii)};
            dt2plot(ii, 1) = min(cL(:,4) - 6*sqrt(cL(:,7)));
            dt2plot(ii, 2) = max(cL(:,4) + 6*sqrt(cL(:,7)));
            dt2plot(ii, 3) = min(cL(:,9)- 1.96*sqrt(cL(:,10)));
            dt2plot(ii, 4) = max(cL(:,9)+ 1.96*sqrt(cL(:,10)));
            dt2plot(ii, 5) = mean(cL(:,4));
            dt2plot(ii, 6) = mean(cL(:,11));
        end
        
        ListProf = {};
        ListFile      = {};
        for ii = 1:length(obj.Replicates)
            ii
            for jj = 1:1%obj.Replicates{ii}.nbrReplicates
                F2L = fullfile(obj.Replicates{ii}.path2fin{jj}, 'myFinnee.mat');
                ListFile{ii, jj} = F2L{1};
                load(F2L{1});
                for kk = 1:size(dt2plot, 1)
                    ListProf{kk, ii, jj} = ...
                        myFinnee.Datasets{options.prm4}.getProfile([dt2plot(kk, 3) dt2plot(kk, 4)],...
                        'XLim', [dt2plot(kk, 1) dt2plot(kk, 2)], 'SetSize', 1);
                end
            end
        end
        
        for ii = 1:size(ListProf, 1)
            InputFig = figure(                       ...
                'Visible'          , 'on'           , ...
                'Units'            , 'normalized'   ,...
                'Position'         , [0.25 0.25 0.5 0.5]);
            
            AxisH1 = axes(InputFig ,...
                'Units'            , 'normalized'     , ...
                'OuterPosition'    , [0 0 2/3 1]);
            
            for jj = 1:size(ListProf, 2)
                for kk = 1:size(ListProf, 3)
                    hold on
                    if any(jj == [1, 2, 3, 4, 5])
                        plot(AxisH1, ListProf{ii, jj, kk}.Data(:,1), ListProf{ii, jj, kk}.Data(:,2), 'k')
                    else
                        plot(AxisH1, ListProf{ii, jj, kk}.Data(:,1), ListProf{ii, jj, kk}.Data(:,2), 'r')
                    end
                    hold off
                    title('Red: Control; Black: Other');
                    xlabel([AxisX.Label, ' / ', AxisX.Unit]);
                    ylabel([AxisY.Label, ' / ', AxisY.Unit]);
                end
            end
            
            
            cL    = allFiltData{Id2Markers(ii)};
            
            String4Edit{1} = '  FIGURES OF MERITS';
            String4Edit{2} = '_____________________';
            String4Edit{3} = '';
            String4Edit{4} = sprintf(['Peak center   : '...
                , AxisX.fo, ' %s'], mean(cL(:,6)),  AxisX.Unit);
            String4Edit{5} = sprintf(['Peak variance: '...
                , '%.3e', ' %s^2'], mean(cL(:,7)), AxisX.Unit);
            String4Edit{6} = sprintf(['m/z = ', obj.Replicates{1}.PeakList.AxisY{1}.fo, ...
                ' +/- ' obj.Replicates{1}.PeakList.AxisY{1}.fo], mean(cL(:,11)), 1.96*std(cL(:,11)));
            
            EditH1 = uicontrol(InputFig ,...
                'Style'              , 'Edit'      , ...
                'Visible'            , 'on'        , ...
                'Units'              , 'normalized', ...
                'Max'                , 2           , ...
                'String'             , String4Edit ,...
                'HorizontalAlignment', 'left', ...
                'Position'           , [2/3 0 1/3 1]);
            
            saveas(InputFig, fullfile(options.prm5, ['figure', num2str(ii), '.jpg']));
        end
        
end

% FIn = subplot(2, 1, 1);
% hold on
% Fot = subplot(2, 1, 2);
% hold on
% assignin('base', 'allFiltData', allFiltData)
%
% for ii = 1:nrpli
%     cId = find(strcmp(listOfTags, obj.Replicates{ii}.Tag));
%
%     p1{ii} = plot(FIn, BPP{ii}(:,1), BPP{ii}(:,3), 'Color', 'k');
%     p1{ii}.Tag = 'plot';
%
%     p2{ii} = plot(Fot, BPP{ii}(:,1), BPP{ii}(:,4), 'Color','k');
%     p2{ii}.Tag = 'plot';
%
% end
%
% linkaxes([FIn, Fot],'x')
% s = {};
%
% ind2sct = data4ctgm(:, 4) == 0;
% if any(ind2sct)
%     s{1} = stem(FIn, data4ctgm(ind2sct, 2), data4ctgm(ind2sct, 3), ...
%         'Color', 'k',...
%         'MarkerFaceColor','red');
%     s{1}.Tag = 'stem';
% end
%
% axis(FIn);
% FIn.Title.String = 'Data remaining (Base Peak Profile)';
% FIn.XAxis.Label.String = [infoX.Label, ' / ', infoX.Unit];
% FIn.YAxis.Label.String = [infoZ.Label, ' / ', infoZ.Unit];
% hold off;
%
% axis(Fot);
% title('Data removed (Base Peak Profile)')
% xlabel([infoX.Label, ' / ', infoX.Unit]);
% ylabel([infoZ.Label, ' / ', infoZ.Unit]);
% hold off;
%
% FigInput = gcf;
% FigInput.Name = sprintf('Study centergram');
% dcm_obj = datacursormode(FigInput);
% set(dcm_obj,'UpdateFcn',@myupdatefcn);
% figure(FigInput)
%
% ToglH1 = uicontrol(FigInput ,...
%     'Style'   , 'togglebutton', ...
%     'Visible' , 'on'        , ...
%     'Units'   , 'normalized', ...
%     'String'  , 'popup PIP on',...
%     'Tag'     , 'ToglH1    ',...
%     'Callback', @pushMe     , ...
%     'Position', [2/3 0 1/3 1]);
% wPB1            = ToglH1.Extent(3);
% hPB1            = ToglH1.Extent(4);
% ToglH1.Position = [1-1.3*wPB1 1-1.3*hPB1  1.2*wPB1  1.2*hPB1];
%
%     function txt = myupdatefcn(~,event_obj)
%         % New data cursor update function
%
%         tgt = get(event_obj,'Target');
%         pos = get(event_obj,'Position');
%         switch tgt.Tag;
%             case 'plot'
%                 xString = [infoX.Label, ' = ', num2str(pos(1), foX), ...
%                     ' ', infoX.Unit];
%                 yString = [infoZ.Label, ' = ', num2str(pos(2), foZ), ...
%                     ' ', infoZ.Unit];
%                 txt = {xString, yString};
%             case 'stem'
%                 M1   = pos(1);
%                 IMax = pos(2);
%                 Id   = find(data4ctgm(:,2) == M1 & data4ctgm(:,3) == IMax);
%                 cFOM = [];
%                 IdXc = IdX(data4ctgm(Id, 1)).Id;
%                 Id4rc = find(~isnan(IdXc));
%                 for sf1 = 1:length(Id4rc)
%                     IdPi = IdXc(Id4rc(sf1));
%                     cPIP = obj.Replicates{Id4rc(sf1)}.PeakList.LstPIP{1}{IdPi};
%                     cFOM = [cFOM; cPIP.FOM];
%                 end
%
%                 fo   = ['%.', num2str(signFig(std(cFOM(:,4)))), 'f'];
%                 fstl = ['Centers: ',  num2str(mean(cFOM(:,4)), fo), ...
%                     ' (', num2str(std(cFOM(:,4)), fo), ...
%                     ' n=', num2str(sf1),') ', infoX.Unit];
%                 fo   = ['%.', num2str(signFig(std(cFOM(:,9)))), 'f'];
%                 sndl = ['Acc. Mass: ',  num2str(mean(cFOM(:,9)), fo), ...
%                     ' (', num2str(std(cFOM(:,9)), fo), ...
%                     ' n=', num2str(sf1),') '];
%                 thrl = ['Intensities: '];
%                 for sf1 =  1:length(IdX(data4ctgm(Id, 1)).Int)
%                     thrl = [thrl, num2str(IdX(data4ctgm(Id, 1)).Int(sf1)), '|'];
%                 end
%                 txt = {fstl, sndl, thrl};
%
%                 if ToglH1.Value == 0
%                     aPIP = {};
%                     for sf1 = 1:length(Id4rc)
%                         IdPi = IdXc(Id4rc(sf1));
%                         aPIP{sf1} = obj.Replicates{Id4rc(sf1)}.PeakList.LstPIP{1}{IdPi};
%                     end
%
%                     mergePIP(aPIP)
%                 end
%         end
%     end
%
%     function pushMe(source, ~)
%         if source.Value == 1
%             source.String = 'popup PIP off';
%         else
%             source.String = 'popup PIP on';
%         end
%     end

end

