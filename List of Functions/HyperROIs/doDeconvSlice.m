function FOM_2 = doDeconvSlice(FOM)

count = 0;
ThreCmpt = 95;

for ii = 3:size(FOM.DeconvMS.deconMS, 2)
    
    for jj = 1:size(FOM.DeconvMS.ModelMZ{ii}, 2)
        cMS = [FOM.DeconvMS.axisMS{ii}(:, 1), FOM.DeconvMS.ModelMZ{ii}(:, jj)];
        
        % check if center far from edges
        [~, IdMax] = max(cMS(:,2));
        D2Edges = min(IdMax, size(cMS, 1) - IdMax);
        if D2Edges > 3
            % CHECK SIGNAL/TO NOISE
            all_Prof = [FOM.DeconvMS.axisTm{ii}(:, 1) ...
                squeeze(FOM.DeconvMS.deconMS{ii}(:, jj, :))];
            sum_Prof = sum(all_Prof(:, 2:end), 2);
            Thres = prctile(sum_Prof, 25);
            Signal = max(sum_Prof) - Thres;
            Noise  = 3*std(sum_Prof(sum_Prof <= Thres));
            
            if Signal/Noise >= 10
                sum_Prof(sum_Prof <= 0.001*max(sum_Prof)) = 0;
                trace = [];
                bool_trace = false;
                start = 0;
                for kk = 1:length(sum_Prof)
                    
                    if sum_Prof(kk) > 0
                        if ~bool_trace
                            bool_trace = true;
                            start = kk;
                        end
                        
                    elseif bool_trace
                        bool_trace = false;
                        if kk - start > 10
                            trace(end+1, :) = [start kk];
                        end
                    end
                end
                
                if ~isempty(trace)
                    for kk = 1:size(trace, 1)
                        count = count + 1;
                        FOM_2.Slice.Id_HyperROIs(count) = ii;
                        XY = all_Prof(trace(kk, 1):trace(kk, 2), :);
                        lm = LocalMaxima([XY(:,1) sum(XY(:, 2:end), 2)], 1, max(sum(XY(:, 2:end), 2))*0.05);
                        FOM_2.Slice.lm{count} = lm;
                        
                        pca_Data = XY(:, 2:end);
                        pca_Data(isnan(pca_Data)) = 0;
                        warning off
                        [~,~,~, ~, explained] = pca(pca_Data, 'algorithm','als');
                        warning on
                        NbrSpct = find(cumsum(explained) >= ThreCmpt, 1, 'first');
                        if isempty(NbrSpct), NbrSpct = 1; end
                        FOM_2.Slice.NbrProf_pca{count} = NbrSpct;
                        FOM_2.Slice.axisMS{count} = cMS(:, 1);
                        FOM_2.Slice.ModelMZ{count} = cMS(:, 2:end);
                        FOM_2.Slice.axisTm{count} = XY(:, 1);
                        FOM_2.Slice.Profiles{count} = XY(:, 2:end);
                    end
                end
            end
            
        end
    end
end