%% DESCRIPTION
%
%% Copyright
% BSD 3-Clause License
% Copyright 2021 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [info] = mkMnROI_h5(obj, mz_Int, tm_Int, FullFile)
% TODO: help and description

narginchk(3, 4)
info = {};

if nargin == 3
    [file, path, ~] = uiputfile('*.h5','Create HDF5 fike');
    if isequal(file,0)
        disp('User selected Cancel');
        return
        
    else
        % TODO: check indx and extension .h5 present
        % TODO: Check if file already exist, if yes ask to delete and
        % replace
        disp(['User selected ', fullfile(path, file)]);
        FullFile =  fullfile(path, file);
    end
end

AxisX     = obj.AxisX.Data;
AxisY     = obj.AxisY.Data;

% NOTE: %HDF5 create with High-Level Functions
h5create(FullFile, '/Axis_Time', size(AxisX))
h5write(FullFile, '/Axis_Time', AxisX)
h5writeatt(FullFile, '/Axis_Time', 'Label', obj.AxisX.Label)
h5writeatt(FullFile, '/Axis_Time', 'Unit', obj.AxisX.Unit)

h5create(FullFile, '/Axis_mz', size(AxisY))
h5write(FullFile, '/Axis_mz', AxisY)
h5writeatt(FullFile, '/Axis_mz', 'Label', obj.AxisY.Label)
h5writeatt(FullFile, '/Axis_mz', 'Unit', obj.AxisY.Unit)

h5writeatt(file, '/', 'DateOfCreation', datestr(datetime))
h5writeatt(FullFile, '/', 'File_Of_Origin',  obj.Path2Fin)
h5writeatt(FullFile, '/', 'Log_Of_Transformation',  obj.Log)

h5create(FullFile, '/ROIs_Limits', [size( mz_Int, 1) 5])
h5write(FullFile, '/ROIs_Limits', [(1:size( mz_Int, 1))', mz_Int, tm_Int])
h5writeatt(FullFile, '/ROIs_Limits', 'Label', 'ID, mz_start, mz_end, tm_start, tm_end') 


indMZ = zeros(size(mz_Int));
Data4ROI{size(mz_Int, 1)} = [];
X{size(mz_Int, 1)} = [];
Y{size(mz_Int, 1)} = [];

for ii = 1:size(mz_Int, 1)
    ix = findCloser(mz_Int(ii, 1), AxisY); indMZ(ii, 1) = max(1, ix);
    ix = findCloser(mz_Int(ii, 2), AxisY);
    indMZ(ii, 2) = min(ix, length(AxisY));
    X{ii} = AxisY(indMZ(ii, 1):indMZ(ii, 2));
end


TLim      = [min(tm_Int(:,1)), max(tm_Int(:,2))];
IdTLim1 = find(AxisX < TLim(1), 1, 'last');
if isempty(IdTLim1), IdTLim1 = 1; end

IdTLim2 = find(AxisX >  TLim(2), 1, 'first');
if isempty(IdTLim2), IdTLim2 = length(AxisX); end

IdTLim = [IdTLim1 IdTLim2];

IsOpenRois = true(size(mz_Int, 1), 1);
switch obj.Format
    case 'profile'
        if ~isempty(obj.AxisY.Data)
            
            h = waitbar(0,'Making ROI, please wait');
            for ii = IdTLim(1):IdTLim(2)
                waitbar(ii/length(AxisX(:,1)))
                
                XMS = xpend(obj, obj.ListOfScans{ii});
                
                id = find(tm_Int(:,1) <= AxisX(ii) & tm_Int(:, 2) >= AxisX(ii));
                for jj = 1:length(id)
                    Data4ROI{id(jj)}(:, end+1) = XMS.Data(indMZ(id(jj),1):indMZ(id(jj),2), 2);
                    Y{id(jj)}(end+1) =  AxisX(ii);
                end
                
                id = find(tm_Int(:, 2) < AxisX(ii) & IsOpenRois);
                for jj = 1:length(id)
                    loc = sprintf('/ROIs/ROI4ID#%i', id(jj));
                    fullData = [[0, Y{id(jj)}]...
                        ; [X{id(jj)}, Data4ROI{id(jj)}]];
                    h5create(FullFile, loc, size(fullData))
                    h5write(FullFile, loc, fullData)
                    h5writeatt(FullFile, loc, 'X_Label', obj.AxisX.Label)
                    h5writeatt(FullFile, loc, 'X_Unit', obj.AxisX.Unit)
                    h5writeatt(FullFile, loc, 'Y_Label', obj.AxisY.Label)
                    h5writeatt(FullFile, loc, 'Y_Unit', obj.AxisY.Unit)
                    h5writeatt(FullFile, loc, 'Time_Interval', tm_Int(id(jj), :))
                    h5writeatt(FullFile, loc, 'mz_Interval', mz_Int(id(jj), :))
                    
                    M = ChrMoment3D...
                        (fullData(1, 2:end)', fullData(2:end, 1),fullData(2:end, 2:end));
                    h5writeatt(FullFile, loc, 'Statistical_Moments_Values', table2array(M))
                    h5writeatt(FullFile, loc, 'Statistical_Moments_Labels'...
                        , 'M0\time, M1\time, M2\time, M3\time, M0\mz, M1\mz, M2\mz, M3\mz');
                    
                    IsOpenRois(id(jj)) = false;
                    Data4ROI{id(jj)} = [];
                    X{id(jj)} = [];
                    Y{id(jj)} = [];
                end
                
            end
            
            
        end
        
    case 'centroid'
        error('L37')
end
if ishandle(h), close(h); end

end
