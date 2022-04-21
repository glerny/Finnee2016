function OptiMe = OptimisedMe_mtd2(HyperROIs)

% definition
TBL = load(HyperROIs.links2dat);
links2dat = TBL.links2dat; clear TBL;
links2dat.S2N = nan(size(links2dat, 1), 1);
links2dat.alignement = zeros(size(links2dat, 1), 3);
minCntPts = 10;

%% 2D SG filter
wx_filter = 5;
wy_filter = 3;
nx_filter = 2;
ny_filter = 2;
h_filter = sgsdf_2d...
    (-wx_filter:wx_filter, -wy_filter:wy_filter, nx_filter, ny_filter);

% TODO: parfor me
for ii = 1:size(links2dat, 1)
    file2add = links2dat.Datafiles{ii};
    while 1
        [fidRead, errmsg]  = fopen(file2add, 'r');
        if fidRead ~= -1
            break
        end
    end
    
    fseek(fidRead,  0, 'bof');
    data = fread(fidRead,inf, links2dat.format{ii});
    Size = links2dat.Size;
    X_Tm = data(1:Size(ii,2));
    Y_mz = data(Size(ii,2)+1:Size(ii,2)+Size(ii,1));
    id1 = Size(ii,2)+Size(ii,1)+1;
    cHROI = reshape(data(id1:end), Size(ii,1), Size(ii,2), []);
    fclose(fidRead);
    
    data = [X_Tm, squeeze(mean(cHROI, 1))];
    isNotValid = false(size(data, 2)-1, 1);
    for jj = 2:size(data, 2)
        Vector = data(:, jj);
        Vector(:,2) = circshift(Vector(:,1) , 1);
        Vector(:,3) = circshift(Vector(:,1) , 2);
        
        if max(diff([0; find(sum(Vector, 2, 'omitnan') == 0); size(Vector, 1)+1]))...
                <= minCntPts
            isNotValid(jj-1) = true;
        end
        
        [IntMax(jj-1), ID] = max(data(:, jj), [], 'omitnan');
        TmMax(jj-1) = data(ID, 1);
    end
    
    [Aligned_Data, Opt] = doAlignment_MinPearson(data, true(size(data, 2)-1, 1));
    Opt(isNotValid, :) = NaN;
    INan = isnan(Aligned_Data);
    Aligned_Data(INan) = 0;
    XY = [Aligned_Data(:,1), mean(Aligned_Data(:, 2:end), 2)];
    
    clear Noi
    for jj = 1:size(Aligned_Data, 2)-1
        Noi(jj) = ...
            std(Aligned_Data(:,jj+1)' - Aligned_Data(:,jj+1)'/XY(:,2)'*XY(:,2)');
    end
    OptiMe.ID(ii, 1) = ii;
    OptiMe.isValid(ii, 1) = sum(~isNotValid);
    OptiMe.meanS2N(ii, 1) = mean(max(Aligned_Data(:, 2:end), [], 'omitnan')./(3*Noi), 'omitnan');
    OptiMe.minS2N(ii, 1) = min(max(Aligned_Data(:, 2:end), [], 'omitnan')./(3*Noi), [], 'omitnan');
    OptiMe.meanCC(ii, 1) = mean(Opt(:,2), 'omitnan');
    OptiMe.minCC(ii, 1) = min(Opt(:,2),[], 'omitnan');
    OptiMe.meanMaxInt(ii, 1) =  mean(IntMax, 'omitnan');
    OptiMe.stdMaxInt(ii, 1) =  std(IntMax, [], 'omitnan');
    OptiMe.meanTmATMI(ii, 1) =  mean(TmMax, 'omitnan');
    OptiMe.stdTmATMI(ii, 1) =  std(TmMax, 'omitnan');
    OptiMe.Opt{ii, 1} = Opt;
    OptiMe.INV{ii, 1} = isNotValid;
    OptiMe.ITN{ii, 1} = [IntMax' TmMax'];
end
OptiMe = struct2table(OptiMe);
