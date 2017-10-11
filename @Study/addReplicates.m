%% DESCRIPTION
%
%% Copyright
% BSD 3-Clause License
% Copyright 2016-2017 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function obj = addReplicates(obj, pklIn, tag, name, parameters, varargin)

% AT THE MOMENT DO THE CALIBRATION BUT DON'T USE IT FOR OPTIMISATION
% PARAMETERS

options.verified = true;
options.adjusted = true;
options.CI       = 2.58;
options.merged   = true;
options.percRpts = 0.75;

CI   = options.CI;
m    = length(obj.Replicates) + 1;
rpts = length(pklIn);
obj.Replicates{m}.nbrReplicates = rpts;
for ii = 1:rpts
    obj.Replicates{m}.path2fin{ii} = pklIn{ii}.Path2Fin;
    obj.Replicates{m}.Log2crea{ii} = pklIn{ii}.Log2crea;
end

obj.Replicates{m}.Tag = tag;
if ~any(strcmp(tag, obj.ListOfTags))
    Id = length(obj.ListOfTags) + 1;
    obj.ListOfTags{Id} = tag;
end

obj.Replicates{m}.name = name;

if rpts > 1
    
    P    = strsplit(parameters, ':');
    Dmz  = str2double(P{1});
    RSDT = str2double(P{2});
    
    PIP2match = [];
    dt4norm   = [];
    
    for ii = 1:rpts
        cPIP = pklIn{ii}.FOM{1}.Data(:, [1, 5, 10, 2]);
        cPIP(:, end+1) = ii;
        PIP2match = [cPIP; PIP2match];
    end
    
    PIP2match = sortrows(PIP2match, 3);
    Id = [0; find(diff(PIP2match(:,3)) > Dmz); size(PIP2match, 1)];
    for ii = 1:length(Id)-1
        Cut = PIP2match (Id(ii)+1:Id(ii+1), :);
        if size(unique(Cut(:,5)), 1) == rpts
            Cut = sortrows(Cut, 2);
            Id2 = [0; find(diff(Cut(:,2))./Cut(1:end-1, 2)*100 > RSDT); size(Cut, 1)];
            for jj = 1:length(Id2)-1
                Cut2 = Cut (Id2(jj)+1:Id2(jj+1), :);
                if size(unique(Cut2(:,5)), 1) == rpts && size(Cut2(:,2), 1) == rpts
                    dt4norm(:,:, end+1) = Cut2;
                end
            end
        end
    end
    dt4norm(:,:,1) = [];
    
    % mz normalisation
    dt4mz = squeeze(dt4norm(:, 3,:));
    mstMz = mean(dt4mz);
    for ii = 1:rpts
        p{ii,1} = polyfitweighted(dt4mz(ii, :), mstMz-dt4mz(ii,:), 3);
    end
    dt4std = max(dt4mz) - min(dt4mz);
    DelMz = mean(dt4std) + 2*CI*std(dt4std);
    fprintf('\nOptimized Dmz = %.4f\n', DelMz) 
    
    % Time normalisation
    dt4time = squeeze(dt4norm(:, 2,:));
    mstTime = mean(dt4time);
    for ii = 1:rpts
        p{ii,2} = polyfitweighted(dt4time(ii, :), mstTime-dt4time(ii,:), 3);
    end
    dt4std =  max(dt4time) - min(dt4time);
    DelTm = mean(dt4std) + 2*CI*std(dt4std);
    fprintf('Optimized DTm = %.2f\n', DelTm) 
    
    %  Intensity normalisation
    
    dt4Int = squeeze(dt4norm(:, 4,:));
    dt4std = std(dt4Int)./mean(dt4Int);
    RsdInt =  (mean(dt4std) + 2*CI*std(dt4std))*100;
    
    % NEW LOOP WITH OPTIMIZED PARAMETERS
    PIP2match = [];
    IxPIP     = {};
    PIPIn     = {};
    
    for ii = 1:rpts
        cPIP = pklIn{ii}.FOM{1}.Data(:, [1, 5, 10, 2, 6]);
        cPIP(:, end+1) = ii;
        PIP2match = [cPIP; PIP2match];
        PIPIn{ii} = false(size(cPIP,1), 1);
    end
    
    dt2keep   = [];
    PIP2match = sortrows(PIP2match, 3);
    Id = [0; find(diff(PIP2match(:,3)) > DelMz); size(PIP2match, 1)];
    
    for ii = 1:length(Id)-1
        Cut = PIP2match (Id(ii)+1:Id(ii+1), :);
        if size(unique(Cut(:,6)), 1) >= int32(rpts*options.percRpts)
            Cut = sortrows(Cut, 2);
            Id2 = [0; find(diff(Cut(:,2)) > DelTm); size(Cut, 1)];
            for jj = 1:length(Id2)-1
                Cut2 = Cut (Id2(jj)+1:Id2(jj+1), :);
                if size(unique(Cut2(:,6)), 1)>= int32(rpts*options.percRpts)
                    if size(unique(Cut2(:,6)), 1) == size(Cut2(:,6), 1)
                        dt2add = nan(rpts, 6);
                        dt2add(Cut2(:,6),:) = Cut2;
                        dt2keep(:,:, end+1) = dt2add;
                    else
                        
                        nbrClst = 0;
                        for ll = 1:rpts
                            nbrClst = max(nbrClst, sum(Cut2(:,6) == ll));
                        end
                        
                        % Normalisaing data 4 k-mean
                        D4KM = Cut2(:, 2:4);
                        D4KM(:,1) = D4KM(:,1)/DelTm;
                        D4KM(:,2) = D4KM(:,2)/DelMz;
                        D4KM(:,3) = D4KM(:,3)./(RsdInt/100*D4KM(:,3));
                        [idx,~] = kmeans(D4KM,nbrClst);
                        
                        IdC = unique(idx);
                        for mm = 1:length(IdC)
                            Cut3 = Cut2 (idx == IdC(mm), :);
                            if size(unique(Cut3(:,6)), 1)  >= int32(rpts*options.percRpts)
                                if size(unique(Cut3(:,6)), 1) == size(Cut3(:,6), 1)
                                    dt2add = nan(rpts, 6);
                                    dt2add(Cut3(:,6),:) = Cut3;
                                    dt2keep(:,:, end+1) = dt2add;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    dt2keep(:,:, 1) = [];
    for ii = 1:rpts
        PIPinCOM(:,ii) = squeeze(dt2keep(ii,1,:));
    end
    
    for ii = 1:rpts
        BPP{ii}(:,1) = pklIn{ii}.AxisX{1}.Data;
        BPP{ii}(:,4) = 0;
        TIP{ii}(:,1) = pklIn{ii}.AxisX{1}.Data;
        TIP{ii}(:,4) = 0;
        
        for jj = 1:length(pklIn{ii}.LstPIP{1})
            IdS = pklIn{ii}.LstPIP{1}{jj}.IdS;
            y   = pklIn{ii}.LstPIP{1}{jj}.y;
            IdE = IdS +  size(y, 1) - 1;
            BPP{ii}(IdS:IdE, 2) = max(BPP{ii}(IdS:IdE, 2), y);
            TIP{ii}(IdS:IdE, 2) = TIP{ii}(IdS:IdE, 2)+ y;
            if ~isempty(find(PIPinCOM(:,ii) == jj))
                BPP{ii}(IdS:IdE, 3) = max(BPP{ii}(IdS:IdE, 3), y);
                TIP{ii}(IdS:IdE, 3) = TIP{ii}(IdS:IdE, 3)+ y;
            else
                BPP{ii}(IdS:IdE, 4) = max(BPP{ii}(IdS:IdE, 4), y);
                TIP{ii}(IdS:IdE, 4) = TIP{ii}(IdS:IdE, 4)+ y;
            end
        end
    end
    
    if options.verified
        
        figure; hold on
        for ii = 1:rpts
            [~, name] = fileparts(pklIn{ii}.Path2Fin{1});
            Legend{ii} = name;
            plot(TIP{ii}(:,1), TIP{ii}(:,3))
        end
        hold off
        legend(Legend)
        title('Remaining data')
        
        figure; hold on
        for ii = 1:rpts
            plot(TIP{ii}(:,1), TIP{ii}(:,4))
        end
        hold off
        legend(Legend)
        title('Filtered out data')
    end
    
    % sorting dt2keep by replicates
    for ii = 1:size(dt2keep, 3)
        dt2keep(:,:, ii) = sortrows(squeeze(dt2keep(:,:, ii)), 6);
    end
    obj.Replicates{m}.Summary = dt2keep;
    
    % Create the peak List
    myPeakList = pklIn{1};
    myPeakList.Type = 'merged';
    myPeakList.p4norm = p;
    
    for ii = 1:rpts
        AxisX{ii} = pklIn{ii}.AxisX{1}.Data + polyval(p{ii,2}, pklIn{ii}.AxisX{1}.Data);
        myPeakList.Path2Fin{ii} = pklIn{1}.Path2Fin{1};
        myPeakList.Log2crea{ii} = pklIn{1}.Log2crea{1};
    end
    myPeakList.AxisX{1} = Axis( myPeakList.AxisX{1}.InfoAxis, AxisX{1});
    
    BPPm(:,1) = AxisX{1};
    BPPm(:,2) = 0;
    TIPm(:,1) = AxisX{1};
    TIPm(:,2) = 0;
    mergedPIP{size(dt2keep, 3)} = {};
    myPeakList.FOM{1}.Headings = {'Id', 'IntMax', 'Tm@IM', 'M0', 'M1',...
        'M2', 'M3', 'mean(m/z)', 'std(m/z)', 'Acc. Mass', 'DeltaInt'};
    myPeakList.FOM{1}.Data = zeros(size(dt2keep, 3), 11);
    
    for jj = 1:size(dt2keep, 3)
        data2merged.Int  = zeros(size(AxisX{1}));
        data2merged.mass = zeros(size(AxisX{1}));
        data2merged.ntz  = zeros(size(AxisX{1}));
        tFOM = [];
        
        for ii = 1:rpts
            IdT = dt2keep(ii, 1,jj);
            if isnan(IdT)
                disp('ttt')
            else
            cData = pklIn{ii}.LstPIP{1}{IdT}.Data;
            
            tFOM = [tFOM; pklIn{ii}.LstPIP{1}{IdT}.FOM];
            cData(:,1) = cData(:,1) + polyval(p{ii,1}, cData(:,1));
            cData(:,3) = AxisX{ii}(cData(:,3));
            
            tm = unique(cData(:,3));
            uData = zeros(length(tm), 3);
            for kk = 1:length(tm)
                Id2B = cData(:,3) == tm(kk);
                uData(kk,3) = tm(kk);
                uData(kk,2) = sum(cData(Id2B, 2));
                uData(kk,1) = sum(cData(Id2B, 1).*cData(Id2B, 2))/sum(cData(Id2B, 2));
            end
            
            data2add = interp1(uData(:,3), uData(:,1), AxisX{1}, 'linear');
            data2add(isnan(data2add)) = 0;
            data2merged.mass = data2merged.mass + data2add;
            data2add = interp1(uData(:,3), uData(:,2), AxisX{1}, 'linear');
            data2add(isnan(data2add)) = 0;
            data2merged.Int  = data2merged.Int + data2add;
            nz               = find(data2add~=0);
            data2merged.ntz(nz)  = data2merged.ntz(nz)+1;
            end
        end
        
        data2merged.mass = data2merged.mass./ data2merged.ntz;
        data2merged.Int  = round(data2merged.Int/rpts);
        data2merged.Int(data2merged.Int < 0) = 0;
        BPPm(:,2) = max(BPPm(:,2), data2merged.Int);
        TIPm(:,2) = TIPm(:,2) + data2merged.Int;
        IdNZ = find(data2merged.Int ~= 0);
        try
            mergedPIP{jj} = pklIn{1}.LstPIP{1}{dt2keep(1, 1,jj)};
            mergedPIP{jj}.IdS = IdNZ(1);
        catch
            disp('WTF')
        end
        mergedPIP{jj}.Data = ...
            [data2merged.mass(IdNZ), data2merged.Int(IdNZ), IdNZ];
        mergedPIP{jj}.x =  AxisX{1}(min(IdNZ):max(IdNZ));
        
        myPeakList.FOM{1}.Data(jj, :) = [jj, mean(tFOM,1,'omitnan')];
        myPeakList.sFOM               = [nan, std(tFOM,1,'omitnan')];
        myPeakList.nFOM               = [nan, sum(~isnan(t))];
    end
    
    InfoBPP = pklIn{1}.BPP{1}.InfoTrc;
    InfoBPP.Loc = 'inTrace';
    myPeakList.BPP{1} = Trace(InfoBPP, BPPm);
    InfoTIP = pklIn{ii}.TIP{1}.InfoTrc;
    InfoTIP.Loc = 'inTrace';
    myPeakList.TIP{1} = Trace(InfoBPP, TIPm);
    myPeakList.LstPIP{1} = mergedPIP;
    
    
    obj.Replicates{m}.PeakList = myPeakList; %.PeakListReview;
else
    obj.Replicates{m}.PeakList = pklIn{1};
end

if isempty(obj.Path2Std)
    obj.Path2Std = uigetdir(pwd, 'Select the folder of destination');
    if ~ischar(obj.Path2Std)
        error('myApp:argChk', 'Cancel by user');
    end
end


save(obj);

end

