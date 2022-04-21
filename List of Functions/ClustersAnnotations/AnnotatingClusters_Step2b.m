function Formula = AnnotatingClusters_Step2(MSSpectra, Adducts, Data, NeutralLoss, p, Delta, BestIons, ClustersIOns)

Formula.Clust    = [];
Formula.mf       = {};
Formula.Possible = {};

Rp = 60000;
ClustersIOns = sortrows(ClustersIOns, 'ID');
% flatten the data
allFOMs = table();
for ii = 1:size(ClustersIOns, 1)
    ClusterID = ii*ones(size(ClustersIOns.FOM{ii}, 1), 1);
    allFOMs = [allFOMs; [table(ClusterID) ClustersIOns.FOM{ii}]];
end
FolderOut = uigetdir();
Options.plot = 'best'; % 'all'/'best'/'none'

%Correcting the Data
allFOMs.Corrmz = allFOMs.mean_AccMass - polyval(p, allFOMs.mean_AccMass).*allFOMs.mean_AccMass/1000000;

% Creation of the Master mz axis
mzLim = [min(allFOMs.Corrmz) max(allFOMs.Corrmz)];
aMZ   = mzLim(1)-4*mzLim(1)/(Rp*2.354);
step  = (mzLim(1)/(Rp*2.354*5));

while 1
    aMZ(end+1) = aMZ(end)+step;
    
    if aMZ(end) > mzLim(2)+4*mzLim(2)/(Rp*2.354)
        break
    end
end

ix = MSSpectra.BPmz >= aMZ(1) & MSSpectra.BPmz <= aMZ(end);
MSSpectra = MSSpectra(ix, :);
clear ix

LstClusters = unique(allFOMs.ClusterID);
for cl =  1:length(LstClusters)
    
    if any(BestIons.ClusterID == cl) && strcmp(Options.plot, 'all')
        FOM = sortrows(allFOMs(allFOMs.ClusterID == LstClusters(cl), :), 'mean_Area','descend');
        cTrace = ClustersIOns.Trace{cl};
        h = figure('units','normalized','outerposition',[0 0 1 1]);
        subplot(2, 2, 1)
        plot(cTrace.Data(:,1), cTrace.Data(:,2))
        xlabel([cTrace.AxisX.Label ' (' cTrace.AxisX.Unit ')']);
        ylabel('Absorbance');
        title(cTrace.Title);
        
        subplot(2, 2, 2)
        bar(FOM.mean_AccMass, FOM.mean_Area, 0)
        xlabel('m/z');
        ylabel('Peak Area');
        xtickformat('%.3f')
        ytickformat('%.0f')
        title('Bar plot of ions in cluster');
        
        close(h)
        saveas(h, fullfile(FolderOut, ['cluster #', num2str(cl), '.png']))
        
    elseif any(BestIons.ClusterID == cl) && ~strcmp(Options.plot, 'None')
        I2C = BestIons(BestIons.ClusterID == cl, :);
        FOM = allFOMs(allFOMs.ClusterID == cl, :);
        I2C = sortrows(I2C, 'Coeff');
        
        for ii = 1:size(I2C, 1)
            FOM = sortrows(allFOMs(allFOMs.ClusterID == LstClusters(cl), :), 'mean_Area','descend');
            cTrace = ClustersIOns.Trace{cl};
            
            h = figure('units','normalized','outerposition',[0 0 1 1]);
            subplot(4, 2, [1 3])
            plot(cTrace.Data(:,1), cTrace.Data(:,2))
            xlabel([cTrace.AxisX.Label ' (' cTrace.AxisX.Unit ')']);
            ylabel('Absorbance');
            ax = gca;
            ax.FontSize = 12;
            title(cTrace.Title, 'FontSize', 14);
            xtickformat('%.2f')
            ytickformat('%.0f')
            
            subplot(4, 2, [2 4])
            bar(FOM.mean_AccMass, FOM.mean_Area, 0)
            xlabel('m/z');
            ylabel('Peak Area');
            xtickformat('%.3f')
            ytickformat('%.0f')
            ax = gca;
            ax.FontSize = 12;
            title('Bar plot of ions in cluster', 'FontSize', 14);
            
            MS2algn = I2C.MS{ii};
            MSLim = [min(MS2algn(:,1) - (3*Delta)*MS2algn(:,1)/1000000) ...
                max(MS2algn(:,1) + (3*Delta)*MS2algn(:,1)/1000000)];
            MSLim(1) = MSLim(1) - 100*MSLim(1)/(Rp*2.354);
            MSLim(2) = MSLim(2) + 4*MSLim(1)/(Rp*2.354);
            XY4Ref  =  centr2prof(MS2algn, Rp,  aMZ(aMZ >= MSLim(1) & aMZ <= MSLim(2)));
            XY4Ref(2, :) = XY4Ref(2, :)/max(XY4Ref(2, :))*100;
            subplot(4, 2, 5)
            plot(XY4Ref(1, :), XY4Ref(2, :), 'r','Linewidth',1)
            ylabel('Rel. int. (%)')
            ytickformat('%.0f')
            MS2match = [I2C.Ions{ii}(:,1)+I2C.Shift(ii) I2C.Ions{ii}(:,2)];
            hold on
            bar(MS2match(:,1), MS2match(:,2)/max(MS2match(:,2))*100, 0);
            ax = gca;
            hold off
            xlim(MSLim)
            XY4Alf =  centr2prof(MS2match, Rp,  aMZ(aMZ >= MSLim(1) & aMZ <= MSLim(2)));
            XY4Alf(2, :) = XY4Alf(2, :)/(XY4Alf(2, :) /XY4Ref(2, :));
            ax = gca;
            ax.FontSize = 12;
            set(ax,'xticklabel',[])
            legend([I2C.mf{ii} ' (' num2str(I2C.Coeff(ii)) ')'], 'FontSize', 12);
            
            subplot(4, 2, 7)
            plot(XY4Ref(1, :), XY4Ref(2, :) - XY4Alf(2, :), 'b','Linewidth',1)
            xlabel('m/z')
            ylabel('Rel. int. (%)')
            xtickformat('%.3f')
            ytickformat('%.0f')
            ax = gca;
            ax.FontSize = 12;
            xlim(MSLim)
            
            %%%% FIND ADDUCTS %%%%
            IM = fndAllMotif(I2C.mf{ii}, Data, Adducts);
            
            PosMetabo = table();
            TotAdd = {};
            
            for jj = 1:size(IM)
                AdductsPattern = mkAddPtrn(IM.mf{jj}, Data, Adducts, NeutralLoss);
                [~, iAP] = unique(AdductsPattern.mf);
                AdductsPattern = AdductsPattern(iAP, :);
                AdductsPattern.mi = AdductsPattern.mi./AdductsPattern.Z;
                AdductsPattern.IDFeature = nan(size(AdductsPattern, 1), 1);
                AdductsPattern.AccuMass = nan(size(AdductsPattern, 1), 1);
                AdductsPattern.Area = nan(size(AdductsPattern, 1), 1);
                
                for kk = 1:size(FOM, 1)
                    i2FOM = find(abs(FOM.Corrmz(kk) - AdductsPattern.mi)./AdductsPattern.mi*1000000 <= Delta);
                    
                    if length(i2FOM) == 1
                         AdductsPattern.IDFeature(i2FOM) = FOM.ClusterID(kk);
                         AdductsPattern.AccuMass(i2FOM) = FOM.mean_AccMass(kk);
                         AdductsPattern.Area(i2FOM) = FOM.mean_Area(kk);
                    elseif length(i2FOM) > 1
                        error('pp')
                    end
                end
                
                disp('pp');
                AdductsPattern(isnan(AdductsPattern.AccuMass), :) = [];
                [~, i2max] = max(AdductsPattern.Area);
                TotalAdducts = size(AdductsPattern, 1);
                PosMetabo = [PosMetabo; [AdductsPattern(i2max, :), table(TotalAdducts)]];
                TotAdd{jj} = AdductsPattern;
            end
            
            Formula.Clust(end+1, 1)    = cl;
            Formula.mf{end+1, 1}       =  I2C.mf{ii};
            Formula.Possible{end+1, 1} = PosMetabo;
            
            subplot(4, 2, [6 8]);
            plot(3);
            pos = get( subplot(4, 2, [6 8]),'position');
            delete( subplot(4, 2, [6 8]))
            
            TString = mkString(PosMetabo);
            annotation(gcf,'Textbox','String',TString, 'FontSize', 9, 'Position', pos);
            
            saveas(h, fullfile(FolderOut, ['cluster #' num2str(cl) '_' num2str(ii) '.png']))
            close(h)
        end
        
    elseif any(BestIons.ClusterID == cl) && strcmp(Options.plot, 'None')
        FOM = allFOMs(allFOMs.ClusterID == cl, :);
        I2C = sortrows(I2C, 'Coeff');
        
        for ii = 1:size(I2C, 1)
            FOM = sortrows(allFOMs(allFOMs.ClusterID == LstClusters(cl), :), 'mean_Area','descend');
            cTrace = ClustersIOns.Trace{cl};
            
            %%%% FIND ADDUCTS %%%%
            IM = find(strcmp(IonsInESI.mf, I2C.mf{ii}));
            
            PosMetabo = table;
            for jj = 1:size(IM)
                Im = find(strcmp(IonsInESI.Motif, IonsInESI.Motif{IM(jj)}));
                pCut = IonsInESI(Im ,:);
                pCut.AccMass = NaN(size(pCut,1), 1);
                pCut.Area    = NaN(size(pCut,1), 1);
                
                for kk = 1:size(pCut,1)
                    I1 = find((abs(FOM.Corrmz - pCut.BPmz(kk))./...
                        pCut.BPmz(kk)*1000000 <= 3*Delta));
                    if ~isempty(I1)
                        if I1 == 1
                            if ~strcmp(pCut.mf{kk}, I2C.mf{ii}), continue; end
                        end
                        pCut.AccMass(kk) = FOM.mean_AccMass(I1);
                        pCut.Area(kk) = FOM.mean_Area(I1);
                        
                    end
                    
                end
                pCut(isnan(pCut.Area), :) = [];
                PosMetabo = [PosMetabo; pCut];
            end
            
            Formula.Clust(end+1, 1)    = cl;
            Formula.mf{end+1, 1}       =  I2C.mf{ii};
            Formula.Possible{end+1, 1} = PosMetabo;
        end
        
    end  
end
Formula = struct2table(Formula);
    function myText = mkString(myArray)
        
        myText{1} = sprintf('POSSIBLE MOLECULES');
        myText{2} = sprintf('##Motif##<nbrIons>: [Adduct (n=?,Z=?)](mz error(ppm)) [Adduct (n=?,Z=?)](mz error(ppm)) ...');
        
        % REMOVE MULTIPLON
        LstMotifs = unique(myArray.Motif);
        for mT1 = 1:length(LstMotifs)
            ILC = find(strcmp(myArray.Motif, LstMotifs{mT1}));
            if length(ILC) > 1
                [mf, ia, ic] = unique(myArray.mf(ILC));
                
                if length(mf) < length(ILC)
                    if length(mf) == 1
                        myArray.AccMass(ILC(2:end)) = NaN;
                    else
                        disp('hya')
                    end
                end
            end
        end
        myArray(isnan(myArray.AccMass), :) = [];
        %%% END REMOVE MULTIPLON
        
        [mf, ia, ic] = unique(myArray.Motif);
        for mT1 = 1:length(mf)
            frq(mT1) = length(unique(myArray.BPmz(strcmp(myArray.Motif, mf{mT1}))));
        end
        
        [~, iS] = sort(frq, 'descend');
        for mT1 = 1:length(iS)
            myText{mT1+3} = sprintf('##%s##<%.0f>: ', ...
                mf{iS(mT1)}, frq(iS(mT1)));
            Ix = find(strcmp(myArray.Motif, mf{iS(mT1)}));
            [~, ia] = sort((myArray.Area(Ix)), 'descend');
            Ix = Ix(ia);
            
            for mT2 = 1:length(Ix)
                myText{mT1+3} = [myText{mT1+3} ...
                    sprintf(' [%s (n=%.0f,Z=%.0f)](%.4f %.1f)', ...
                    myArray.Adducts{Ix(mT2)}, myArray.Repetition(Ix(mT2)),...
                    myArray.Charge(Ix(mT2)), myArray.BPmz(Ix(mT2)), ...
                    (myArray.BPmz(Ix(mT2)) - myArray.AccMass(Ix(mT2)))/myArray.BPmz(Ix(mT2))*1000000)];
            end
        end
    end

    function XYout = centr2prof(XYin, Rp, axis)
        XYout = zeros(size(axis));
        for lp1 = 1: size(XYin, 1)
            XYout = XYout + XYin(lp1, 2)*...
                pdf('normal', axis, XYin(lp1, 1), (XYin(lp1, 1)/Rp)/2.354);
        end
        XYout = [axis; XYout];
    end
end