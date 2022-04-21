%% DESCRIPTION
%
%% Copyright
% BSD 3-Clause License
% Copyright 2021 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function myROIs = mkMnROI_dat(obj, target, FullFolder)
% TODO: help and description

narginchk(2, 3)
myROIs = {};
target = sortrows(target,'ID','ascend');

if nargin == 2
    [file, path] = uiputfile('*.roi', 'Create ROI''s folder');
    if isequal(file,0)
        disp('User selected Cancel');
        return
        
    else
        % TODO: check indx and extension .h5 present
        % TODO: Check if file already exist, if yes ask to delete and
        % replace
    end
end

FullFolder = fullfile(path, file);
if ~exist(FullFolder, 'dir')
    mkdir(FullFolder);
else
%TODO: if full folfer exist ask if delete is require
end

myROIs.AxisX      = obj.AxisX; AxisX = obj.AxisX.Data;
myROIs.AxisY      = obj.AxisY; AxisY = obj.AxisY.Data;
myROIs.myPath     = fullfile(FullFolder, 'myROIs.mat');
myROIs.links2dat  = fullfile(FullFolder, 'links2dat.mat');
myROIs.Path2Dat   = fullfile(FullFolder, 'datafile.dat');
myROIs.FinneeFile = obj.Path2Fin;
myROIs.Dataset    = obj.Log;

nrow                    = size(target, 1);
links2dat.ID            = target.ID;
links2dat.mz_Interval   = [target.mz_min, target.mz_max];
links2dat.time_Interval = [target.Time_min, target.Time_max];
links2dat.size          = zeros(nrow, 2);
links2dat.format        = repmat({''}, nrow, 1);
links2dat.position      = zeros(nrow, 2);
links2dat.ChrMoment     = {};

indMZ = zeros(size(target, 1), 2);
Data4ROI{size(target, 1)} = [];
X{size(target, 1)} = [];
Y{size(target, 1)} = [];

for ii = 1:size(target, 1)
    ix = findCloser(target.mz_min(ii), AxisY); indMZ(ii, 1) = max(1, ix);
    ix = findCloser(target.mz_max(ii), AxisY);
    indMZ(ii, 2) = min(ix, length(AxisY));
    X{ii} = AxisY(indMZ(ii, 1):indMZ(ii, 2));
end


TLim      = [min(target.Time_min), max(target.Time_max)];
IdTLim1 = find(AxisX < TLim(1), 1, 'last');
if isempty(IdTLim1), IdTLim1 = 1; end

IdTLim2 = find(AxisX >  TLim(2), 1, 'first');
if isempty(IdTLim2), IdTLim2 = length(AxisX); end

IdTLim = [IdTLim1 IdTLim2];
IsOpenRois = true(size(target, 1), 1);
[fidWriteDat, errmsg]  = fopen(myROIs.Path2Dat, 'wb');

switch obj.Format
    case 'profile'
        if ~isempty(obj.AxisY.Data)
            
            h = waitbar(0,'Making ROI, please wait');
            for ii = IdTLim(1):IdTLim(2)
                waitbar(ii/length(AxisX(:,1)))
                
                XMS = xpend(obj, obj.ListOfScans{ii});
                id = find(target.Time_min <= AxisX(ii) & target.Time_max >= AxisX(ii));
                for jj = 1:length(id)
                    Data4ROI{id(jj)}(:, end+1) = XMS.Data(indMZ(id(jj),1):indMZ(id(jj),2), 2);
                    Y{id(jj)}(end+1) =  AxisX(ii);
                end
                
                id = find(target.Time_max < AxisX(ii) & IsOpenRois);
                for jj = 1:length(id)
                    fullData = [[0, Y{id(jj)}]...
                        ; [X{id(jj)}, Data4ROI{id(jj)}]];
                    M = ChrMoment3D...
                        (fullData(1, 2:end)', fullData(2:end, 1),fullData(2:end, 2:end));
                    
                    links2dat.size(id(jj), :)     = size(fullData);
                    links2dat.format{id(jj)}      = 'double';
                    links2dat.position(id(jj), 1) = ftell(fidWriteDat);
                    fwrite(fidWriteDat, fullData(:), 'double');
                    links2dat.position(id(jj), 2) = ftell(fidWriteDat);
                    links2dat.ChrMoment{id(jj), 1} = M;
                    
                    IsOpenRois(id(jj)) = false;
                    Data4ROI{id(jj)} = [];
                    X{id(jj)} = [];
                    Y{id(jj)} = [];
                end
            end
            
            id = find(IsOpenRois);
            for jj = 1:length(id)
                fullData = [[0, Y{id(jj)}]...
                    ; [X{id(jj)}, Data4ROI{id(jj)}]];
                M = ChrMoment3D...
                    (fullData(1, 2:end)', fullData(2:end, 1),fullData(2:end, 2:end));
                
                links2dat.size(id(jj), :)     = size(fullData);
                links2dat.format{id(jj)}      = 'double';
                links2dat.position(id(jj), 1) = ftell(fidWriteDat);
                fwrite(fidWriteDat, fullData(:), 'double');
                links2dat.position(id(jj), 2) = ftell(fidWriteDat);
                links2dat.ChrMoment{id(jj), 1} = M;
                IsOpenRois(id(jj)) = false;
                Data4ROI{id(jj)} = [];
                X{id(jj)} = [];
                Y{id(jj)} = [];
            end
            
            
        end
        
    case 'centroid'
        error('L37')
end
if ishandle(h), close(h); end
fclose(fidWriteDat)
links2dat = struct2table(links2dat);
save(myROIs.links2dat, 'links2dat')
save(myROIs.myPath, 'myROIs')

end
