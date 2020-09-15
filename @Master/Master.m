%% DESCRIPTION
% ROI is the class that deals with three-dimensional representations.
% Those can a cut from the dataset after mzAxis alignment to the Matser mz Axis
% others.
%
%% LIST OF THE CLASS'S PROPERTIES
%
%% LIST OF THE CLASS'S METHODS
%
%% Copyright
% BSD 3-Clause License
% Copyright 2016-2017 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Master
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Path
        Name
        QC      = {};
        Samples = {};
        NPL     = 'myPeakList.mat'
    end
    
    methods
        function obj = Master(cMz, cTm, WR, nthre)
            %% 1- INITIALISATION
            % 1.1. Constants
            
            [obj.Name, obj.Path] = uiputfile('*.mat','Location of master file');
            cArea   = 2;
            thrSN   = 10;
            thrCC   = 0.7;
            TLim    = [0 inf];
            Mw      = 50000;
            
            % 1.2. Load QC samples
            fprintf('\n\n STARTING \n')
            dirs = uigetdirs(pwd, 'select Finnee QC folders');
            
            figure
            hold on
            Lim = [inf, 0];
            for ii = 1:length(dirs)
                PL{ii} = load(fullfile(dirs{ii}, obj.NPL));
                plot(PL{ii}.myPeakList.BPP{1}.Data(:,1), PL{ii}.myPeakList.BPP{1}.Data(:,2))
                PL{ii}.myPeakList.FOM{1}.Data(:,1) = ...
                    1:length(PL{ii}.myPeakList.FOM{1}.Data(:,1));
                myPkLsts{ii} = PL{ii}.myPeakList.FOM{1}.Data;
                obj.QC.Files{ii} = dirs{ii};
                axis{ii} = PL{ii}.myPeakList.AxisX{1, 1}.Data;
                Lim = [min(Lim(1), min(axis{1})) max(Lim(2), max(axis{1}))];
            end
            
            p = polyfull(axis{1}(2:end), diff(axis{1}), 2, []);
            AxisTm = Lim(1);
            while 1
                AxisTm(end+1) = AxisTm(end) + polyval(p, AxisTm(end));
                if AxisTm(end) > Lim(2), break; end
                
            end
            hold off
            repli   = ii;
            QCnbr   = ii;
            
            %% 2-. Finding markers without tm correction
            % 2.1. CONSTRAINED CLUSTERING USING NEAREST NEIGHBOURG
            Cluster = ConsClust(myPkLsts);  
            
            %% 3- Check profiles in each cluster
            ii     = 1;
            Mcc    = [];
            allPrf = {};
            id2del =  size(Cluster);
            
            while 1
                if ii > length(Cluster)
                    break
                end
                
                Profiles = AxisTm';
                % 3.1 DO test
                cData = Cluster{ii};
                if size(cData, 1) > 1
                    for jj = 1:size(cData, 1)
                        I2PeL = cData(jj, 9);
                        I2PIP = cData(jj, 1);
                        cPIP = PL{I2PeL}.myPeakList.LstPIP{1}{I2PIP};
                        
                        x = cPIP.x;
                        vq = interp1(x, cPIP.y, Profiles(:,1));
                        vq(isnan(vq)) = 0;
                        Profiles(:, end+1) = vq;
                    end
                    PrfMinus = circshift(Profiles(:, 2:end), -1);
                    PrfPlus  = circshift(Profiles(:, 2:end), 1);
                    Profiles(sum([PrfMinus, Profiles(:, 2:end), PrfPlus] == 0, 2) == 3*jj, :) = [];
                    
                    if any(sum(Profiles) == 0)
                        disp('HISE')
                    end
                    
                    CC = corrcoef([mean(Profiles(:, 2:end), 2), Profiles(:, 2:end)]);
                    if isempty(CC)
                        do('wtfh')
                    end
                    
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
            
            % 3.2. Calculate FOM
            FOM = zeros(length(Cluster), 11);
            for ii = 1:length(Cluster)
                vect = zeros(1, repli - QCnbr);
                cData = Cluster{ii};
                FOM(ii, 1) = ii;
                FOM(ii, 2) = sum(cData(:,9) <= QCnbr);
                FOM(ii, 3) = mean(cData(:,3));
                FOM(ii, 4) = std(cData(:,3));
                FOM(ii, 5) = mean(cData(:,6));
                FOM(ii, 6) = std(cData(:,6));
                idc = cData(:,9) <= QCnbr;
                FOM(ii, 7) = mean(cData(idc,2));
                FOM(ii, 8) = std(cData(idc,2));
                FOM(ii, 9) = FOM(ii, 8)/FOM(ii, 7)*100;
                FOM(ii, 10) = Mcc(ii,1);
                FOM(ii, 11) = Mcc(ii,2);
            end
            
            % 3.3. Filter bad clusters
            minRplt = ceil(nthre* QCnbr/100);
            Ix = find(FOM(:,2) < minRplt);
            FOM(Ix, :) = [];
            Cluster(Ix) = [];
            allPrf(Ix)  = [];
            FOM(:,1) = 1:length(Cluster);
            
            %% 4- Save and quit
            obj.QC.parameters.WR      = WR;
            obj.QC.parameters.cTm     = cTm;
            obj.QC.parameters.cMz     = cMz;
            obj.QC.parameters.cArea   = cArea;
            obj.QC.parameters.nthre   = nthre;
            obj.QC.parameters.thrSN   = thrSN;
            obj.QC.parameters.thrCC   = thrCC;
            obj.QC.parameters.TLim    = TLim;
            obj.QC.parameters.nbrQC   = QCnbr;
            obj.QC.parameters.minRplt = minRplt;
            
            
            obj.QC.Method1.Cluster = Cluster;
            obj.QC.Method1.FOM    = array2table(FOM,...
                'VariableNames',{'ID', 'Replicates', ...
                'CtrTime', 'std_CtrTime', ...
                'AccurateMass', 'std_AccurateMass', ...
                'Area', 'std_ARea', 'RSD'...
                'Mean_PearsonCoef', 'Min_PearsonCoef'});
            obj.QC.Method1.Profiles   = allPrf;
            
            obj.QC.Axis.AxisX = Axis(PL{1}.myPeakList.AxisX{1}.InfoAxis, AxisTm');
            obj.QC.Axis.AxisY = PL{1}.myPeakList.AxisY{1};
            
            myMaster = obj;
            save(fullfile(obj.Path, obj.Name), 'myMaster')
            
            %% NESTED FUNCTIONS
            
            % NF1- COnstrained clustering
            function cluster = ConsClust(PLs)
                
                %1. randomise the ordre of selection of replicates
                IX = randperm(length(PLs));
                cluster = {};
                cFOM    = [];
                
                %2. All features in the first datasets start a cluster
                
                data = PLs{IX(1)}...
                    (PLs{IX(1)}(:, 12) >= thrSN, [1, 4:7, 10, 12, 11]);
                data(:, end+1) = IX(1);
                for ccii = 1:size(data, 1)
                    cluster{ccii} = data(ccii, :);
                    cFOM(ccii, :) = [1, data(ccii, [2, 3, 6])];
                end
                
                %3. Then do the remaining datasets
                for ccii = 2:length(IX)
                    data =  PLs{IX(ccii)}...
                        (PLs{IX(ccii)}(:, 12) >= thrSN, [1, 4:7, 10, 12, 11]);
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
            
            % NF2- CONCATENATE CLUSTER THAT OVERLAP AT 4xSIGMA
            function [FOM, Cluster] = Concatenate(FOM, Cluster)
                
                while 1
                    
                    redo = false;
                    cii = 1;
                    while 1
                        mStm = median(FOM(FOM(:,2) >= 4,4));
                        mSmz = median(FOM(FOM(:,2) >= 4,6));
                        FOM(:,9) = FOM(:,4);
                        FOM(:,10) = FOM(:,6);
                        %FOM(FOM(:,9) < mStm , 9) = mStm;
                        %FOM(FOM(:,10) < mSmz , 10) = mSmz;
                        IdCct =...
                            find(abs(FOM(cii, 3) - FOM(:,3)) <= 3*(FOM(cii,9) + FOM(:,9)) & ...
                            abs(FOM(cii, 5) - FOM(:,5)) <= 3*(FOM(cii, 10) + FOM(:, 10)));
                        
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
                        end
                        cii = cii+1;
                        if cii > size(FOM, 1), break, end
                    end
                    
                    if ~redo
                        break
                    else
                        disp('one more cycle')
                    end
                end
                FOM = FOM(:, 1:8);
                
            end
        end
    end
end

