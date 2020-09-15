function [Isotopomers, clusterFOM]  = doGroupIons(myMaster, CC)

cHA  = myMaster.QC.Method3.HACA;
cFOM = myMaster.QC.Method3.FOM;
cPro = myMaster.QC.Method3.Profiles;

Isotopomers = {};
for ii = 1:length(cHA)
    
    clear Thread
    xZ = cHA{ii}.Z;
    xFOM = cHA{ii}.FOM;
    
    for jj = 1:size(xFOM, 1)
        Thread{jj} = xFOM(jj, :);
        
    end
    
    if isempty(xZ)
        Ix = 0;
    else
        Ix = find(xZ(:,4) >= CC, 1, 'last');
        if isempty(Ix), Ix = 0; end
    end
  
    for jj = 1:Ix
            Thread{xZ(jj, 3)} = [Thread{xZ(jj, 1)}; Thread{xZ(jj, 2)}];
            Thread{xZ(jj, 1)} = [];
            Thread{xZ(jj, 2)} = [];
    end
    
    Thread = Thread(find(~cellfun(@isempty, Thread)));
    Isotopomers = [Isotopomers, Thread];
end

% Checking if some clusters corrolate at more than 0.95.
AxisTime = myMaster.QC.Axis.AxisX.Data;
while 1
    
    allProf  = zeros(length(AxisTime), size(Isotopomers, 1));
    for ii = 1:length(Isotopomers)
        cThread = Isotopomers{ii};
        cProf = zeros(length(AxisTime), size(cThread, 1));
        
        for jj = 1:size(cThread, 1)
            XY = myMaster.QC.Method3.Profiles{cThread.ID(jj)};
            [C,ia,ib] = intersect(XY(:, 1), AxisTime);
            cProf(ib, jj) = XY(:, 2);
        end
        
        allProf(:, ii) = sum(cProf, 2);
    end
    
    R = corrcoef(allProf);
    max_CC = max(max(triu(R,1)));
    
    if max_CC < CC, break; end
    
    [linMin, colMin] = find(triu(R,1) == max_CC);
    Isotopomers{end+1} = [Isotopomers{linMin(1)}; Isotopomers{colMin(1)}];
    Isotopomers([linMin(1), colMin(1)]) = [];
end
% Getting CLusters' FOM
clusterFOM = {};

for ii = 1:length(Isotopomers)
    clusterFOM{ii, 1}         = ii;
    clusterFOM{ii, 2}         = size(Isotopomers{ii}, 1);
    [clusterFOM{ii, 3}, Ix2M] = max(Isotopomers{ii}.Area);
    clusterFOM{ii, 4}         = Isotopomers{ii}.AccurateMass(Ix2M);
    clusterFOM{ii, 5}         = mean(Isotopomers{ii}.CtrTime);
    clusterFOM{ii, 6}         = std(Isotopomers{ii}.CtrTime);
    
end

clusterFOM = cell2table(clusterFOM,...
    'VariableNames',{'ID2allThread' 'nbrIon' 'BIArea' 'BIam' 'time' 'std_time'});

end

