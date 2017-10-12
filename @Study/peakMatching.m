%% DESCRIPTION
%
%% Copyright
% BSD 3-Clause License
% Copyright 2016-2017 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [output, Tag] = peakMatching(obj, parameters, varargin)

options.CI4PM     = 2.58;
% STEP 1. Peak matching
% 1.1. Find peaks common to all dataset using original parameters
% 1.1.1. Load default parameters

CI   = options.CI4PM;
P    = strsplit(parameters, ':');
Dmz  = str2double(P{1});
RSDT = str2double(P{2});

% 1.1.2. Build a matrix will all FOM
nrpli = length(obj.Replicates);
allFOM = [];
for ii = 1:nrpli
    FOM = obj.Replicates{ii}.PeakList.FOM{1}.Data;
    FOM(:, end+1) = ii;
    allFOM = [allFOM; FOM];
    Tag{ii} = [obj.Replicates{ii}.Tag, '/', ...
        obj.Replicates{ii}.name, '/', ...
        num2str(obj.Replicates{ii}.nbrReplicates)];
end

% 1.1.2. Remove FOM of max Intnesity lower than ThIT Important to improve
% accuracy (SHOULD ADD CONTROL)

allFOM = sortrows(allFOM, 10);
Id = [0; find(diff(allFOM(:,10)) > Dmz); size(allFOM, 1)];
data4norm = [];
for ii = 1:length(Id)-1
    Cut = allFOM (Id(ii)+1:Id(ii+1), :);
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
dt4mz   = squeeze(data4norm(:, 10,:));
dt4std  = max(dt4mz) - min(dt4mz);
DelMz   = mean(dt4std) + 2*CI*std(dt4std);
fprintf('\nOptimized Dmz = %.4f\n', DelMz)

% Time normalisation
dt4time = squeeze(data4norm(:, 5,:));
dt4std  =  max(dt4time) - min(dt4time);
DelTm   = mean(dt4std) + 2*CI*std(dt4std);
fprintf('Optimized DTm = %.2f\n', DelTm)

% 1.2.2 Peak matching
allFOM = sortrows(allFOM, 10);
Id = [0; find(diff(allFOM(:,10)) > DelMz); size(allFOM, 1)];
IdX      = {};
output   = [];

for ii = 1:length(Id)-1
    Cut = allFOM (Id(ii)+1:Id(ii+1), :);
    Cut = sortrows(Cut, 5);
    
    if size(Cut, 1) > 1
        Id2 = [0; find(diff(Cut(:,5)) > DelTm); size(Cut, 1)];
        
        for jj = 1:length(Id2)-1
            Cut2 = Cut (Id2(jj)+1:Id2(jj+1), :);
            
            if size(Cut2, 1) ==  size(unique(Cut2(:,11)), 1)
                data2add             = NaN(nrpli, 12);
                data2add(Cut2(:,12),:) = Cut2;
                output(:, :, end+1)  = data2add;
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
                    data2add             = NaN(nrpli, 12);
                    data2add(Cut3(:,12), :) = Cut3;
                    output(:, :, end+1)  = data2add;
                end
            end
        end
    else
        data2add             = NaN(nrpli, 12);
        data2add(Cut(:,12), :) = Cut;
        output(:, :, end+1)  = data2add;
    end
end

output(:, :, 1) = [];

end

