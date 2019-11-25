function Master = Markers

% I. INITIALISATION
% I.1. Constants
WR      = 170; % weight ratio mz/tm. 0.5 min as the same weight as 0.001 ppm
cMz     = 0.005;    % constrain condition for mz
cTm     = 0.5;      % constrain condition for tm
cArea   = 2;
nthre   = 0.75;
thrSN   = 10;
thrCC   = 0.7;
TLim    = [0 inf];
Master.QC.PeakLists = {};

% I.2. Load QC samples
Master.path      = pwd;
fprintf('\n\n STARTING \n')
dirs = uigetdirs(Master.path, 'select Finnee QC folders');
Legend = {};

figure('Name', 'check me!'),
hold on

for ii = 1:length(dirs)
    myPLQ{ii} = load(fullfile(dirs{ii}, 'myPeakList_bak.mat'));
    BPP = myPLQ{ii}.myPeakList.BPP{1, 1}.Data;
    IdB = strfind(dirs{ii}, '\');
    Name = dirs{ii}(IdB(end)+1:end-4);
    Master.QC.PeakLists{end+1}.ID =  Name;
    Master.QC.PeakLists{end}.path = dirs{ii};
    Master.QC.PeakLists{end}.BPP = BPP;
    data =  myPLQ{ii}.myPeakList.FOM{1}.Data;
    idx = data(:,5) >= TLim(1) & data(:,5) <= TLim(2);
    Master.QC.PeakLists{end}.FOM  = data(idx, :);
    fprintf('%s included\n', Master.QC.PeakLists{end}.ID)
    Legend{end+1} = Name;
    plot(BPP(:,1), BPP(:,2))
end
repli   = ii;
QCnbr   = ii;
fprintf('\n%.0f QC files used\n', repli)
xlabel('time / min');
ylabel('Intensity')
legend(Legend)
title('superimposed BPP from QC samples')
hold off


% II. Finding markers without tm correction
%II.1. CONSTRAINED CLUSTERING USING NEAREST NEIGHBOURG
Cluster = ConsClust(Master);

%III.2. Check profiles in each cluster
ii = 1;
Mcc = [];
allPrf = {};
while 1
    if ii > length(Cluster)
        break
    end
    Profiles = BPP(:,1);
    cData = Cluster{ii};
    if size(cData, 1) > 1
        for jj = 1:size(cData, 1)
            I2PeL = cData(jj, 9);
            I2PIP = cData(jj, 1);
            
            if I2PeL <= QCnbr
                cPIP = myPLQ{I2PeL}.myPeakList.LstPIP{1}{I2PIP};
            else
                cPIP = myPLS{I2PeL-QCnbr}.myPeakList.LstPIP{1}{I2PIP};
            end
            
            x = cPIP.x;
            vq = interp1(x, cPIP.y, Profiles(:,1));
            vq(isnan(vq)) = 0;
            Profiles(:, end+1) = vq;
        end
        Profiles(sum(Profiles(:, 2:end) == 0, 2) == jj, :) = [];
        CC = corrcoef([mean(Profiles(:, 2:end), 2), Profiles(:, 2:end)]);
        
        if length(unique(cData(:,9))) ~= length(cData(:,9)) || min(CC(2:end, 1)) < thrCC
            CC(:,1) = []; CC(1, :) = [];
            [I1, I2] = find(CC == min(min(CC)));
            CC1 = corrcoef([Profiles(:, I1(1)+1), Profiles(:, 2:end)]);
            CC2 = corrcoef([Profiles(:, I2(1)+1), Profiles(:, 2:end)]);
            [~, iz] = max([CC1(2:end, 1), CC2(2:end, 1)], [], 2);
            Cluster{ii} = cData(iz == 1, :);
            Cluster{end+1} = cData(iz == 2, :);
        else
            ii = ii +1;
            Mcc(end+1, 1) = mean(CC(2:end, 1));
            Mcc(end  , 2) = min(CC(2:end, 1));
            allPrf{end+1} = [Profiles(:,1), mean(Profiles(:, 2:end), 2)];
        end
    else
        ii = ii +1;
        Mcc(end+1, 1) = NaN;
        Mcc(end  , 2) = NaN;
        allPrf{end+1} = [];
    end
end

FOM = zeros(length(Cluster), 10);
Areas = zeros(length(Cluster), repli - QCnbr);
for ii = 1:length(Cluster)
    vect = zeros(1, repli - QCnbr);
    cData = Cluster{ii};
    FOM(ii, 1) = ii;
    FOM(ii, 2) = sum(cData(:,9) <= QCnbr);
    FOM(ii, 3) = sum(cData(:,9) > QCnbr);
    FOM(ii, 4) = mean(cData(:,3));
    FOM(ii, 5) = std(cData(:,3));
    FOM(ii, 6) = mean(cData(:,6));
    FOM(ii, 7) = std(cData(:,6));
    idc = cData(:,9) <= QCnbr;
    FOM(ii, 8) = mean(cData(idc,2));
    FOM(ii, 9) = std(cData(idc,2));
    FOM(ii, 10) = FOM(ii, 9)/FOM(ii, 8)*100;
    FOM(ii, 11) = Mcc(ii,1);
    FOM(ii, 12) = Mcc(ii,2);
    idc = cData(:,9) > QCnbr;
    vect(cData(idc,9)-QCnbr) = cData(idc,2);
    Areas(ii, :) = vect;
end


minRplt = ceil(nthre* QCnbr);
Ix = find(FOM(:,2) < minRplt);
FOM(Ix, :) = [];
Areas(Ix, :) = [];
Cluster(Ix) = [];
allPrf(Ix)  = [];
FOM(:,1) = 1:length(Cluster);
Im = findSupClusters(FOM);
FOM(:, end+1)= Im';

Master.paramters.WR      = WR;
Master.paramters.cTm     = cTm;
Master.paramters.cMz     = cMz;
Master.paramters.FOM     = FOM;
Master.paramters.cArea   = cArea;
Master.paramters.nthre   = nthre;
Master.paramters.thrSN   = thrSN;
Master.paramters.thrCC   = thrCC;
Master.paramters.TLim    = TLim;

Master.Markers.Cluster = Cluster;
Master.Markers.FOM     = FOM;
Master.Markers.SAreas  = Areas;
Master.Markers.Profi   = allPrf;
% fprintf('Total markers found (at least common to %.0f QC): \t%.0f\n',  [minRplt length(Cluster)]);
% fprintf('Markers common to the %.0f QC samples: \t%.0f\n',  [repli sum(FOM(:,2) == repli)]);

    function cluster = ConsClust(Mst)
        
        PLs = Mst.QC.PeakLists;
        %1. randomise the ordre of selection of replicates
        IX = randperm(length(PLs));
        cluster = {};
        cFOM    = [];
        
        %2. All features in the first datasets start a cluster
        
        data = PLs{IX(1)}.FOM...
            (PLs{IX(1)}.FOM(:, 12) >= thrSN, [1, 4:7, 10, 12, 11]);
        data(:, end+1) = IX(1);
        for ccii = 1:size(data, 1)
            cluster{ccii} = data(ccii, :);
            cFOM(ccii, :) = [1, data(ccii, [2, 3, 6])];
        end
        
        %3. Then do the remaining datasets
        for ccii = 2:length(IX)
            data =  PLs{IX(ccii)}.FOM...
                (PLs{IX(ccii)}.FOM(:, 12) >= thrSN, [1, 4:7, 10, 12, 11]);
            data(:, end+1) = IX(ccii);
            
            target = zeros(size(cluster));
            for jj = 1:size(data, 1)
                
                %3.a find the closest neighbourgh
                test = cFOM(:, [3, 4]) - data(jj, [3 6]);
                [~, Id4min] = min(sqrt((test(:,1)/WR).^2 +(test(:,2)).^2));
                
                %3.b. check constrains
                if abs(cFOM(Id4min,3) - data(jj, 3)) <= cTm &&...
                        abs(cFOM(Id4min,4) - data(jj, 6)) <= cMz &&...
                        (data(jj,2) >= cFOM(Id4min,2)/cArea) && ...
                        (data(jj,2) <= cFOM(Id4min,2)*cArea)
                    ccData    = [cluster{Id4min};  data(jj, :)];
                    cluster{Id4min}   = ccData;
                    cFOM(Id4min, 2:4) = mean(ccData(:, [2 3 6]), 1);
                    cFOM(Id4min, 1)   = size(ccData,1);
                    target(Id4min)    = target(Id4min)+1;
                    
                else
                    cluster{end+1} = data(jj, :);
                    cFOM(end+1, :) = [1, data(jj, [2, 3, 6])];
                    target(end+1)  = 1;
                    
                end
            end
        end
        nbrQC = ccii;
        
        FOM = zeros(length(cluster), 8);
        for ccii = 1:length(cluster)
            ccData = cluster{ccii};
            FOM(ccii, 1) = ccii;
            FOM(ccii, 2) = size(ccData,1);
            FOM(ccii, 3) = mean(ccData(:,3));
            FOM(ccii, 4) = std(ccData(:,3));
            FOM(ccii, 5) = mean(ccData(:,6));
            FOM(ccii, 6) = std(ccData(:,6));
            FOM(ccii, 7) = mean(ccData(:,2));
            FOM(ccii, 8) = std(ccData(:,2));
            
        end
        
        [FOM, cluster] = Concatenate(FOM, cluster);
        for ccii = 1:length(cluster)
            ccData = cluster{ccii};
            if length(unique(ccData(:,9))) == length(ccData(:,9))
                FOM(ccii, 9) = 1;
            end
            
        end
    end




end

function cluster = Clust_2(Data, cluster, WR)

% find multiplons
Xx      = unique(Data(:,9));
[~, tg] = max(histc(Data(:,9), Xx));
lst = find(Data(:,9) == Xx(tg(1)));

for ll = 1:length(lst)
    ccluster{ll} = Data(lst(ll), :);
    cFOM(ll, :) = [1, Data(lst(ll), [3, 6])];
end
Data(lst, :) = [];

for ll = 1:size(Data, 1)
    test = cFOM(:, [2, 3]) - Data(ll, [3 6]);
    [~, Id4min] = min(sqrt((test(:,1)/WR).^2 +(test(:,2)).^2));
    cD    = [ccluster{Id4min};  Data(ll, :)];
    ccluster{Id4min}   = cD;
    cFOM(Id4min, 2:3) = mean(cD(:, [3 6]), 1);
    cFOM(Id4min, 1)   = size(cD,1);
end

for ll = 1:size(ccluster, 2)
    cD = ccluster{ll};
    if length(unique(cD(:,9))) ~= length(cD(:,9))
        cluster = Clust_2(cD, cluster, WR);
    else
        cluster{end+1} = cD;
    end
end

end

function [FOM, Cluster] = Concatenate(FOM, Cluster)

%III. CONCATENATE CLUSTER THAT OVERLAP AT 4xSIGMA
%!!! CORRECT FORM FOR WRONG SIGMA (IF n < 4);
while 1
    
    redo = false;
    cii = 1;
    while 1
        mStm = median(FOM(FOM(:,2) >= 4,4));
        mSmz = median(FOM(FOM(:,2) >= 4,6));
        FOM(:,9) = FOM(:,4);
        FOM(:,10) = FOM(:,6);
        FOM(FOM(:,9) < mStm , 9) = mStm;
        FOM(FOM(:,10) < mSmz , 10) = mSmz;
        IdCct =...
            find(abs(FOM(cii, 3) - FOM(:,3)) <= 2*(FOM(cii,9) + FOM(:,9)) & ...
            abs(FOM(cii, 5) - FOM(:,5)) <= 2*(FOM(cii, 10) + FOM(:, 10)));
        
        if length(IdCct) > 1
            ccData = [Cluster{IdCct(1)}; Cluster{IdCct(2)}];
            Cluster{IdCct(1)} = ccData;
            FOM(IdCct(1), :) = [IdCct(1), size(ccData, 1), mean(ccData(:,3)),...
                std(ccData(:,3)), mean(ccData(:,6)), std(ccData(:,6)),...
                mean(ccData(:,2)), std(ccData(:,2)), 0, 0];
            FOM(IdCct(1), 9) = max(mStm, FOM(IdCct(1), 4));
            FOM(IdCct(1), 10) = max(mSmz, FOM(IdCct(1), 6));
            FOM(IdCct(2), :)  = [];
            Cluster(IdCct(2)) = [];
            redo = true;
            %disp('merging')
        end
        cii = cii+1;
        if cii > size(FOM, 1), break, end
    end
    
    if ~redo
        break
    else
        %disp('one more cycle')
    end
end
FOM = FOM(:, 1:8);
end

function Ix = findSupClusters(FOM)
mStm = median(FOM(FOM(:,2) >= 4,4));
mSmz = median(FOM(FOM(:,2) >= 4,6));
FOM(FOM(:,4) < mStm , 4) = mStm;
FOM(FOM(:,6) < mSmz , 6) = mSmz;
for ii = 1:size(FOM, 1)
    IdX = find(abs(FOM(ii, 3) - FOM(:,3)) <= 2*(FOM(ii,4) + FOM(:,4)) & ...
        abs(FOM(ii, 5) - FOM(:,5)) <= 2*(FOM(ii,6) + FOM(:,6)));
    
    if length(IdX) > 1
        Ix(ii) = 1;
    else
        Ix(ii) = 0;
    end
end


end
