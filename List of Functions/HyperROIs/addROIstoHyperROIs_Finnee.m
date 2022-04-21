function HyperROIs = addROIstoHyperROIs_Finnee(HyperROIs, TAG)

%%%%!!!!!!!
%%%% TO BE CHANGED BY MEASURINH mzstart AND mzend DURING CENTROID
%%%% TRANSFORMATIOM
%%%% !!!!!
% TODO: Make an object of HyperROIs.

% 1. Introduction
dts    = HyperROIs.parameters.dts;
Limits = HyperROIs.parameters.limits;

if ~strcmp(TAG, unique(HyperROIs.ListOfTags))
    HyperROIs.ListOfTags{end+1} = TAG;
end

%% 2D SG filter
wtm_filter = 3;
wmz_filter = 2;
ntm_filter = 1;
nmz_filter = 1;
h_filter = sgsdf_2d...
    (-wtm_filter:wtm_filter, -wmz_filter:wmz_filter, ntm_filter, nmz_filter);

dirs = uigetdirs(pwd, 'select Finnee Samples folders');
tb = load(HyperROIs.links2dat);
links2dat = tb.links2dat; clear tb

% Separation of elements for parfor loop
ID            = links2dat.ID;
mz_Interval   = links2dat.mz_Interval;
time_Interval = links2dat.time_Interval;
Size          = links2dat.Size;
format        = links2dat.format;
Datafiles     = links2dat.Datafiles;
%ChrMoment     = links2dat.ChrMoment;

for cF = 1:length(dirs)
    [~, N, ~] = fileparts(dirs{cF});
    fprintf('\t adding %s\n', N)
    MF = load(fullfile(dirs{cF}, 'myFinnee.mat'));
    myFinnee = MF.myFinnee; clear MF;
    
    warning off
    %TODO: Check if the file is already there
    HyperROIs.Files.ID(end+1) = HyperROIs.Files.ID(end) + 1;
    HyperROIs.Files.Name{end} = myFinnee.FileID;
    HyperROIs.Files.Tags{end} = TAG;
    HyperROIs.Files.Path{end} = fullfile(dirs{cF}, 'myFinnee.mat');
    warning on
    
    [ROI, X, Y] = myFinnee.Datasets{dts}.mkMnROI_2...
        (Limits.mz_min, Limits.mz_max, Limits.Time_min, Limits.Time_max);
    
    parfor ii = 1:size(Limits, 1)
        file2add = Datafiles{ii};
        [fidWriteDat, errmsg]  = fopen(file2add, 'a+');
        while fidWriteDat == -1
            disp('pp')
            [fidWriteDat, errmsg]  = fopen(file2add, 'a+');
        end
        fseek(fidWriteDat,  0, 'bof');
        data = fread(fidWriteDat,inf, format{ii});
        if isempty(data)
            disp([file2add ' is empty'])
            continue
        end
        X_Tm = data(1:Size(ii,2));
        Y_mz = data(Size(ii,2)+1:Size(ii,2)+ Size(ii,1));
        
        [Xq, Yq] = meshgrid(X_Tm, Y_mz);
        [Xn, Yn] = meshgrid(Y{ii}', X{ii});
        newROI = interp2(Xn, Yn, ROI{ii}, Xq, Yq);
        newROI(isnan(newROI)) = 0;
        fROI = filter2(h_filter, newROI, 'same');
        fROI(isnan(fROI)) = 0;
               
        % TODO: optimise writting size
        fwrite(fidWriteDat, fROI, 'double');
        fclose(fidWriteDat);
    end
end

links2dat = table(ID, mz_Interval, time_Interval, Size, format, Datafiles);
save(HyperROIs.links2dat, 'links2dat')
save(HyperROIs.myPath, 'HyperROIs')




