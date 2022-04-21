function [Aligned_Data, Shifts] = doAlignment3D_MinPearson_WithRef(Data3D, Ref, leeway, dir)

options = optimset('Display','off');
for ii = 1:size(Data3D, 3)
    XY(:, 1) = Ref(:, 1);
    XY(:, 2) = sum(squeeze(Data3D(:, :, ii)), dir, 'omitnan');
    Ori(:, ii) = XY(:, 2);
    
    if dir == 1
        Shifts(ii) = fminsearch(@align, 0, options);
        [~, Cor(:, ii)] = align(Shifts(ii));
        cData = interp1(XY(:, 1) + Shifts(ii), squeeze(Data3D(:, :, ii))', XY(:, 1));
        cData(isnan(cData)) = 0;
        Aligned_Data(:, :, ii) = cData';
        
    elseif dir == 2
        Shifts(ii) = fminsearch(@align, 0, options);
        [~, Cor(:, ii)] = align(Shifts(ii));
        cData = interp1(XY(:, 1) + Shifts(ii), squeeze(Data3D(:, :, ii)), XY(:, 1));
        cData(isnan(cData)) = 0;
        Aligned_Data(:, :, ii) = cData;
    end
end

    function [f, Aligned_data] = align(x)
        if nnz(XY(:, 2)) < 4
            x = 0;
            f = 1;
            Aligned_data = XY(:,2);
            return
        end
        
        Aligned_data = interp1(XY(:, 1) +x , XY(:, 2), XY(:, 1));
        if abs(x) > leeway, Aligned_data = XY(:, 2)*0; end
        Aligned_data(isnan(Aligned_data)) = 0;
        C = corrcoef(Aligned_data, Ref(:, 2));
        f = 1 - C(1, 2);
        
    end

end

