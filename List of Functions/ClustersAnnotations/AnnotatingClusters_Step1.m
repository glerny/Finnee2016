function [BestIOns, ClustersIOns] = AnnotatingClusters_Step1(MSSpectra, ClustersIOns)

Rp = 60000;% Estimated Resolving power of the MS instrument
Thrs_mz = 5;% ppm
ThresCoeff = 0.003;
MzInt = [0.3 3];

BestIOns = table();
ClustersIOns = sortrows(ClustersIOns, 'ID');
FolderOut = uigetdir();
ClustersIOns.CanDo = false(size(ClustersIOns, 1), 1);

if ischar(FolderOut)
    SavePrint = true;
else
    SavePrint = false;
end

% flatten the database
allFOMs = table();
for ii = 1:size(ClustersIOns, 1)
    ClusterID = ii*ones(size(ClustersIOns.FOM{ii}, 1), 1);
    allFOMs = [allFOMs; [table(ClusterID) ClustersIOns.FOM{ii}]];
end
mzLim = [min(allFOMs.mean_AccMass) max(allFOMs.mean_AccMass)];
step  = (mzLim(1)/(Rp*2.354*5));

% Creation of the Master mz axis
aMZ   = mzLim(1)-4*mzLim(1)/(Rp*2.354);
while 1
    aMZ(end+1) = aMZ(end)+step;
    
    if aMZ(end) > mzLim(2)+4*mzLim(2)/(Rp*2.354)
        break
    end
    step  = (aMZ(end)/(Rp*2.354*5));
end


ix = MSSpectra.BPmz >= aMZ(1) & MSSpectra.BPmz <= aMZ(end);
MSSpectra = MSSpectra(ix, :);
clear ix

LstClusters = unique(allFOMs.ClusterID);
for cl = 1:size(ClustersIOns, 1)
    FOM = sortrows(allFOMs(allFOMs.ClusterID == LstClusters(cl), :), ...
        'mean_Area','descend');
    
    ix = FOM.mean_AccMass >= FOM.mean_AccMass(1) + MzInt(1) &...
        FOM.mean_AccMass <= FOM.mean_AccMass(1) + MzInt(2);
    if sum(ix) == 0
        continue
    end
    ClustersIOns.CanDo(cl) = true;
    MSexp = [FOM.mean_AccMass, FOM.mean_Area];
    
    IX = find(abs(MSSpectra.BPmz - FOM.mean_AccMass(1))./MSSpectra.BPmz*1000000 < Thrs_mz);
    if isempty(IX)
        continue
    end
    
    Match = MSSpectra(IX, :);
    Match.Shift = NaN(length(IX), 1);
    Match.Coeff = NaN(length(IX), 1);
    
    for ii = 1: length(IX)
        
        MS2algn = Match.MS{ii};
        MSLim = [min(MS2algn(:,1) - Thrs_mz*MS2algn(:,1)/1000000) ...
            max(MS2algn(:,1) + Thrs_mz*MS2algn(:,1)/1000000)];
        MSLim(1) = MSLim(1) - 1.1;
        MSLim(2) = MSLim(2) + 2.2;
        
        XY4Ref  =  centr2prof(MS2algn, Rp,  aMZ(aMZ >= MSLim(1) & aMZ <= MSLim(2)));
        XY2algn = centr2prof(MSexp, Rp,  aMZ(aMZ >= MSLim(1) & aMZ <= MSLim(2)));
        FOMext  = sortrows(MSexp(MSexp(:,1) >= MSLim(1) &...
            MSexp(:,1) <= MSLim(2), :), 1);
        [ Match.Shift(ii), Match.Coeff(ii)] = ...
            doAlignMSProfile(XY4Ref, XY2algn, Thrs_mz);
        Match.Ions{ii} = FOMext;
    end
    Match = Match(Match.Coeff <= ThresCoeff, :);
    if ~isempty(Match)
        ind = ones(size(Match, 1), 1);
        BestIOns = [BestIOns;  [FOM(ind, [1 2 4 6 8 10]), Match]];
    end
    
    if ~isempty(Match) && SavePrint
        h = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off');
        
        Match = sortrows(Match, 'Coeff');
        for jj = 1:min(2, size(Match, 1))
            
            MS2algn = Match.MS{jj};
            MSLim = [min(MS2algn(:,1) - Thrs_mz*MS2algn(:,1)/1000000) ...
                max(MS2algn(:,1) + Thrs_mz*MS2algn(:,1)/1000000)];
            MSLim(1) = MSLim(1) - 100*MSLim(1)/(Rp*2.354);
            MSLim(2) = MSLim(2) + 4*MSLim(1)/(Rp*2.354);
            
            XY4Ref  =  centr2prof(MS2algn, Rp,  aMZ(aMZ >= MSLim(1) & aMZ <= MSLim(2)));
            XY4Ref(2, :) = XY4Ref(2, :)/max(XY4Ref(2, :))*100;
            subplot(2, 2, jj)
            plot(XY4Ref(1, :), XY4Ref(2, :), 'r','Linewidth',1)
            ax = gca;
            ax.FontSize = 16;
            xlabel('m/z')
            ylabel('Relative intensity')
            xtickformat('%.3f')
            ytickformat('%.1f')
            ax.LineWidth = 1.5;
            
            hold on
            bar(Match.Ions{jj}(:,1)+Match.Shift(jj), Match.Ions{jj}(:,2)/max(Match.Ions{jj}(:,2))*100, 0)
            ax = gca;
            hold off
            title([char(Match.mf{jj}) ' (' num2str(Match.Coeff(jj)) ')'])
            xlim(MSLim)
            
            XY4Alf =  centr2prof([Match.Ions{jj}(:,1)+Match.Shift(jj) Match.Ions{jj}(:,2)], Rp,  aMZ(aMZ >= MSLim(1) & aMZ <= MSLim(2)));
            XY4Alf(2, :) = XY4Alf(2, :)/(XY4Alf(2, :) /XY4Ref(2, :));
            subplot(2, 2, jj+2)
            plot(XY4Ref(1, :), XY4Ref(2, :) - XY4Alf(2, :), 'b','Linewidth',1)
            ax = gca;
            ax.FontSize = 16;
            xlabel('m/z')
            ylabel('Relative intensity')
            xtickformat('%.3f')
            ytickformat('%.1f')
            ax.LineWidth = 1.5;
            title('% residuals')
            xlim(MSLim)
            
        end
        sgtitle(['Cluster #' num2str(cl)])
        saveas(h, fullfile(FolderOut, ['cluster #', num2str(cl), '.png']))
        close(h)
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