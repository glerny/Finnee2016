function obj = getEnvelope(obj, tgtQC, CorelThres)

%% 1- INITIALISATION AND OPTIONS;

switch tgtQC
    case 1
        HACA     = obj.QC.Method1.HACA;
        Profiles = obj.QC.Method1.Profiles;
        
    case 2
        HACA     = obj.QC.Method2.HACA;
        Profiles = obj.QC.Method2.Profiles;
        
    case 3
        HACA     = obj.QC.Method3.HACA;
        Profiles = obj.QC.Method3.Profiles;
        
    otherwise
        error
end

%% 2- Group profiles @ CC = CorelThres
Clusters.FOM = {};
Clusters.Profiles = {};

for ii = 1:length(HACA)
    FOM = HACA{ii}.FOM;
    Z   = HACA{ii}.Z;
    maxI = size(Z, 1) + 1;
    allI = 1:maxI;
    
    if isempty(Z)
        Clusters.FOM{end + 1} = FOM;
        continue
    end
    
    while any(Z(:,4) > CorelThres)
        IX = find(Z(:,4) > CorelThres, 1, 'last');
        [Z, IX2] = unfold(Z, IX);
        if ~isempty(IX2)
            Clusters.FOM{end + 1} = FOM(IX2, :);
            [~, pos] = intersect(allI, IX2);
            allI(pos) = [];
        end
    end
    
    if ~isempty(allI)
        for jj = 1:length(allI)
            Clusters.FOM{end + 1} = FOM(allI(jj), :);
        end
    end
    
end

%% 3- Make and check averaged profiles
AxisTime = obj.QC.Axis.AxisX;
AxisY    = obj.QC.Axis.AxisY;
allProf = zeros(length(AxisTime.Data(:, 1)), length(Clusters.FOM));
for ii = 1:length(Clusters.FOM)
    Vector = AxisTime.Data;
    for jj = 1:size(Clusters.FOM{ii}, 1)
        XY = Profiles{Clusters.FOM{ii}.IDFeature(jj)};
        [~, pos] = intersect(Vector(:,1), XY(:,1));
        Vector(pos, jj+1) = XY(:,2);
    end
    allProf(:, ii) =  sum(Vector(:, 2:end), 2);
end

%% 4- Final verifications
CC = corrcoef(allProf);
while  max(max(triu(CC, 1))) >= CorelThres
    max(max(triu(CC, 1)))
    [col, lin] = find(CC == max(max(triu(CC, 1))));
    allProf(:, col(1)) = sum([allProf(:, col(1)), allProf(:, lin(1))], 2);
    allProf(:, lin(1)) = [];
    Clusters.FOM{col(1)} = [Clusters.FOM{lin(1)}; Clusters.FOM{col(1)}];
    Clusters.FOM(lin(1)) = [];
    CC = corrcoef(allProf);
end

%% 5- Get Figures of Merits

for ii = 1:size(Clusters.FOM, 2)
    FOM = sortrows(Clusters.FOM{ii}, 'mean_Area', 'descend');
    Summary.ID(ii, 1)       = ii;
    Summary.nbrIons(ii, 1)  = size(FOM, 1);
    Summary.BI_Area(ii, 1)  = FOM.mean_Area(1);
    Summary.BI_mz(ii, 1)    = FOM.mean_AccMass(1);
    Summary.Tm_ave(ii, 1)   = mean(FOM.mean_M1);
    Summary.Tm_std(ii, 1)   = std(FOM.mean_M1);
    
    XY          = trailRem([AxisTime.Data, allProf(:,ii)], 2);
    info.Title  = ['clustered ions #' num2str(ii) ' (ions: ' num2str(size(FOM, 1)) ')'];
    info.FT     = '';
    info.TT     = 'SEP';
    info.AxisX  =  Axis(AxisTime.InfoAxis, []);
    info.AxisY  =  Axis(AxisY.InfoAxis, []);
    info.P2Fin  = '';
    info.Loc    = 'intrace';
    info.AdiPrm = {};
    Summary.Trace{ii, 1}   = Trace(info, XY);
    
    Vector = AxisTime.Data;
    for jj = 1:size(FOM, 1)
        XY = Profiles{FOM.IDFeature(jj)};
        [~, pos] = intersect(Vector(:,1), XY(:,1));
        Vector(pos, jj+1) = XY(:,2);
    end
    Corr2MeanP = corrcoef([allProf(:,ii), Vector(:, 2:end)]);
    FOM = addvars(FOM, Corr2MeanP(2:end, 1), 'NewVariableNames', 'Correl2mean');
    Summary.FOM{ii, 1} = FOM;
    
end
    
%% 6- Save and close

switch tgtQC
    case 1
        obj.QC.Method1.IonsClust = struct2table(Summary);
        
    case 2
        obj.QC.Method2.IonsClust = struct2table(Summary);
        
    case 3
        obj.QC.Method3.IonsClust = struct2table(Summary);
        
    otherwise
        error
end

myMaster = obj;
save(fullfile(obj.Path, obj.Name), 'myMaster')


    function [Y, Id] = unfold(Z, IX)
        Y = Z;
        Id1 =  Z(IX, 1:2);
        Y(IX, :) = [];
        LoopMe   = true;
        
        while LoopMe
            
            LoopMe = false;
            I2 = [];
            I3 = [];
            for lp1 = 1:length(Id1)
                I4 = find(Y(:,3) == Id1(lp1));
                if isempty(I4)
                    I2 = [I2, Id1(lp1)];
                else
                    I2 = [I2, Y(I4, 1:2)];
                    I3 = [I3, I4];
                    LoopMe   = true;
                end
            end
            Id1 = I2;
            Y(I3, :) = [];
        end
        Id = Id1;
    end
end

