%% DESCRIPTION
%
%% Copyright 
% BSD 3-Clause License
% Copyright 2016-2017 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


classdef Options4Finnee
    
    properties
        %% Options for fileInfo
        FileIn      = ''
        FolderOut   = ''
        FileID      = ''
        Path2Fin    = ''
        ID4Saving   = 'myFinnee'
        Rounding    = 3
        Overwrite   = false
        FileFormat  = 'mzML'
        MaxFileSize = 5000000000
        
        % Limits
        MZlim       = [0 inf]
        TMLim       = [0 inf]
        
        % Remove spikes
        RemSpks     = false
        SpksSz      = 2
    end
    
    methods
        function obj = Options4Finnee(VAR)
            input = @(x) find(strcmpi(VAR,x),1);
            
            tgtIx               = input('tLim');
            if ~isempty(tgtIx)
                tLim = VAR{tgtIx +1};
                obj.TMLim(1) = min(tLim);
                obj.TMLim(2) = max(tLim);
            end
            tgtIx               = input('mzLim');
            if ~isempty(tgtIx)
                mzLim = VAR{tgtIx +1};
                obj.MZlim(1) = min(mzLim);
                obj.MZlim(2) = max(mzLim);
            end
            tgtIx               = input('rounding');
            if ~isempty(tgtIx)
                obj.Rounding = VAR{tgtIx +1};
            end
            tgtIx               = input('remSpikes');
            if ~isempty(tgtIx)
                obj.RemSpks     = true;
                obj.SpksSz      = VAR{tgtIx +1};
            end
            tgtIx               = input('fileIn');
            if ~isempty(tgtIx);
                obj.FileIn = VAR{tgtIx +1};
            end
            tgtIx               = input('folderOut');
            if ~isempty(tgtIx);
                obj.FolderOut = VAR{tgtIx +1};
            end
            tgtIx               = input('fileID');
            if ~isempty(tgtIx);
                obj.FileID = VAR{tgtIx +1};
            end
            tgtIx               = input('overwrite');
            if ~isempty(tgtIx);
                obj.Overwrite = true;
            end
            
            % fileIn
            if isempty(obj.FileIn)
                ext = 'pwd/*.mzML';
                txtStg = 'Select the mzML file to load';
                [fileName, pathName] = uigetfile(ext, txtStg);
                if ~ischar(fileName) && ~ischar(pathName)
                    error('myApp:argChk', 'User cancel file selection');
                end
                obj.FileIn = fullfile(pathName, fileName);
            end
            
            % folderOut
            if isempty(obj.FolderOut)
                obj.FolderOut = uigetdir(pwd, 'Select the folder of destination');
                if ~ischar(obj.FolderOut)
                    error('myApp:argChk', 'Cancel by user');
                end
            end
            
            % fileID
            if isempty(obj.FileID)
                [~, dfltName, ~] = fileparts(obj.FileIn);
                answer = inputdlg('Enter a name for Finnee', ...
                    'fileID', 1, {dfltName});
                if isempty(answer)
                    error('myApp:argChk', 'Cancel by user')
                elseif isempty(answer{1})
                    error('myApp:argChk', 'Cancel by user')
                end
                obj.FileID = answer{1};
            end
            
            if exist(fullfile(obj.FolderOut, [obj.FileID '.fin']), 'dir') == 7
                if obj.Overwrite
                    rmdir(fullfile(obj.FolderOut, [obj.FileID '.fin']), 's')
                else
                    error('error \nThe directory %s already exist. \nDelete it, change the name or use ''Overwrite''', ...
                        fullfile(obj.FolderOut, [obj.FileID '.fin']));
                end
            end
            obj.Path2Fin = fullfile(obj.FolderOut, [obj.FileID '.fin']);
            mkdir(obj.Path2Fin);
        end
    end
    
end

