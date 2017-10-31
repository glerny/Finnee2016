%% DESCRIPTION
% FINNEE is the entry object of the Finnee2016 toolbox. Finnee is used to
% record all the information about a single experiment and to link with
% the various binary files where all the data are stored. A Finnee object
% may contain various *Datasets* where each dataset will contain all the
% MS scans that described  a full experiment and a series
% of methods allowing to modify, correct and create a particular dataset. Finnee2016 was
% engineered for data recorded using separation techniques hyphenated with
% high-resolution mass spectrometry (X-HRMS
% <https://github.com/glerny/Finnee2016/wiki/Definition-of-terms> ). Finnee 
% only recognised mzML files. For
% more information or questions, visite the (Finnee Blog<>) or the (wiki
% <>)
%
%% LIST OF THE CLASS'S PROPERTIES
% *FileId*        : The generic name for the Finnee object
% DateOfCreation  : The date of creation of the Finne object
% *FileIn*        : The original mzML file
% *Datasets*      : The list of all *datasets* that are linked with this
% object.
% *Options*       : (Hidden) Options used for the creation of the Finnee
%   object.
% *Path2Fin*      : (Hidden) The path to the parent directory
%   (directory with the .fin extension)
% *MZMLDump*      : (Hidden) A copy of the information that was present in
%   the original mzML file.
%
%% LIST OF THE CLASS'S METHODS
% *Finnee*            : The constructor method.
% *queryMZMLDump*     : Not implemented yet
% *Save*              : Record changes
% *BaselineCorrection*: (Profile dataset only). Allow correcting profiles
%   in a dataset from baseline drift.
% *DoCentroid*        : (Profile dataset only). Will transform each profile
%   MS scan in a dataset to centroid MS scan.
% *FilterDataset*     : Filter a dataset to remove unnecessary information.
% *Align2newMZ*       : (Profile dataset only). Adjust every MS scans to a common MZ axis
%
%% REFERENCES
% * Erny, GL, Acunga, T, Simó, C, Alves, A, "Background correction in
%   separation techniques hyphenated to high-resolution mass spectrometry
%   –Thorough correction with MS scans recorded as profile spectra"
%   Journal of Chromatography A
%   doi: <http://dx.doi.org/10.1016/j.chroma.2017.02.052>
% * Erny, GL, Acunga, T, Simó, C, Alves, A, "Finnee—A Matlab toolbox for
%   separation techniques hyphenated high resolution mass spectrometry
%   dataset"
%   Chemometrics and Intelligent Laboratory Systems 155, 138-144
%   doi: <http://dx.doi.org/10.1016/j.chemolab.2016.04.013>
% * Erny, GL, Acunga, T, Simó, C, Cifuentes, A. Alves, A, "Algorithm for
%   comprehensiva analysis of datasets from hyphenated high resolution mass
%   spectrometric techniques using single ion profiles and cluster
%   analysis"
%   Journal of Chromatography A 1429, 134-141
%   doi: <http://dx.doi.org/10.1016/j.chroma.2015.12.005>
% * Erny, GL, Simó, C, Cifuentes, A, Esteves, VI, "Introducing the concept
%   of centergram. A new tool to squeeze data from separation techniques–
%   mass spectrometry couplings"
%   Journal of Chromatography A 1330, 89-96.
%   doi: <http://dx.doi.org/10.1016/j.chroma.2014.01.014>
% * <https://github.com/glerny/Finnee2016/wiki/>
% * <https://finneeblog.wordpress.com/>
%
%% COPYRIGHT
% Copyright BSD 3-Clause License Copyright 2016-2017 G. Erny
% (guillaume@fe.up.pt), FEUP, Porto, Portugal
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Finnee
    
    properties
        FileID         % The generic name
        DateOfCreation % The date of creation
        FileIn         % Name and path of the original mzML file
        Datasets       % The array of datasets linked to this object
        Options        % The options that were used by the constructor
                       % method
        Path2Fin       % The paths to the .fin folder
        MZMLDump       % Some general information from the mzML file
    end
    
    methods
        function obj = Finnee(varargin)
            %% DESCRIPTION
            % FINNEE is the constructor method that decifer the mzML file 
            % and parse the information into the Finnee objet and oter 
            % assoiated objets. The constructor method will create the 
            % first dataset that can either be MS centroid or MS profile 
            % type depending on the format of the mzML file. !FINNEE only 
            % kept MS1 informations, MSn scans will be disarded!. Finnee 
            % can be used without any input parameters or using the 
            % following options:
            %
            % *fileIn    Followed by a string that is the full path of the 
            %            mzML file.
            % *folderOut Followed by a tring that is the path to the 
            %            destination folder for this FINNEE object.
            % *fileID    Followed by a string that is the generic name for 
            %            this object. the Finnee objects and all associated
            %            data files will be recorded in the folder
            %            'folderOut\fileID.fin'.
            % *overwrite delete the folder 'folderOut\fileID.fin' and ALL 
            %            files in this folder if it already exits.
            % *tLim      Followed by a 2x1 array of numbers
            %            (default [0 inf]). Only records scans between 
            %            tLim(1) and tLim(2)
            % *mzLim     Followed by a 2x1 array of numbers
            %            (default [0 inf]). For each scan only kept mz 
            %            values between mzLim(1) and mzLim(2). 
            % *spikes    Followed by an integer between 0 and 3 (default 1)
            %            . (see the method @Finnee\FilterDataset) for
            %            additional information. Remove spikes in every MS 
            %            scans If used, where spikes are any peaks in each
            %            MS of length equal or lower that the integer.
            %            'spikes' followed by 0 allows to to turn off 
            %            spikes removal.
            %
            %% EXAMPLES
            % * myFinnee = Finnee;
            %   Will create a Finnee object without any options. The target
            %   mzML file, destination folder and generic name will be
            %   asked when running the script.
            % *	myFinnee = Finnee('fileIn'   , 'K:\Data\test.mzML', ...
            %                     'folderOut', 'K:Finnee',          ...
            %                     'Overwrite',                      ...
            %                     'fileID'   , 'test1',             ...
            %                     'Spikes'   , 0);
            %   Here all parameters are defined by the optional parameters
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %% CORE OF THE FUNCTION
            
            % 1- Initialisation and options
   
            % Check the options and create the Finnee object
            options = checkVarargin(varargin{:});
            
            % Initialization of Finnee
            obj.DateOfCreation = datetime;
            obj.FileID         = options.FileID;
            obj.FileIn         = options.FileIn;
            obj.Datasets       = {};
            obj.MZMLDump       = {};
            obj.Options        = options;
            obj.Path2Fin       = '';
            
            % Check Options, ask for mzML file/folder and name if the
            % corresponding options are empty
            if isempty(options.FileIn)
                ext = 'pwd/*.mzML';
                txtStg = 'Select the mzML file to load';
                [fileName, pathName] = uigetfile(ext, txtStg);
                if ~ischar(fileName) && ~ischar(pathName)
                    error('myApp:argChk', 'User cancel file selection');
                end
                obj.FileIn = fullfile(pathName, fileName);
            end
            
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
            
            if isempty(options.FolderOut)
                options.FolderOut = uigetdir(pwd, 'Select the folder of destination');
                if ~ischar(options.FolderOut)
                    error('myApp:argChk', 'Cancel by user');
                end
            end
            
            if exist(fullfile(options.FolderOut, [obj.FileID '.fin']), 'dir') == 7
                if options.Overwrite
                    rmdir(fullfile(options.FolderOut, [obj.FileID '.fin']), 's')
                else
                    error('error \nThe directory %s already exist. \nDelete it, change the name or use ''Overwrite''', ...
                        fullfile(options.FolderOut, [obj.FileID '.fin']));
                end
            end
            
            % Create the .fin folder
            obj.Path2Fin = fullfile(options.FolderOut, [obj.FileID '.fin']);
            mkdir(obj.Path2Fin);
            
            % 2- Run the domzML2Finnee 
            % Scripts for other format could be implemented here if needed
            [obj, MSn]      = domzML2Finnee(obj);
            
            % MSn are the MSn scans that were present in the mzML file,
            % they are not used yet; but are accessible via the following
            % variable that will be created in the woekspace.
            assignin('base', 'MSn4myFinnee', MSn);
            
            % 3- Save and done
            myFinnee = obj; %#ok<*NASGU>
            save(fullfile(obj.Path2Fin, 'myFinnee.mat'), 'myFinnee')
            
            %% SUB FUNCTIONS
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %% CHECKVARARGIN
            function options = checkVarargin(varargin)
                % CHECKVARARGIN is used to check the optional paramters
                % and create the options parameter.
                
                % 1- Defaults optional parameters
                options.FileIn      = '';
                options.FolderOut   = '';
                options.FileID      = '';
                options.Path2Fin    = '';
                options.Overwrite   = false;
                options.FileFormat  = 'mzML';
                options.MaxFileSize = 10000000000;
                options.XLim        = [0 inf];
                options.YLim        = [0 inf];
                options.RemSpks     = true;
                options.SpkSz       = 1;
                
                % 2- Decipher varargin and update options when relevamt
                input = @(x) find(strcmpi(varargin,x),1);
              
                tgtIx = input('fileIn');
                if ~isempty(tgtIx);
                    options.FileIn = varargin{tgtIx +1};
                end
                
                tgtIx = input('folderOut');
                if ~isempty(tgtIx);
                    options.FolderOut = varargin{tgtIx +1};
                end
                
                tgtIx = input('fileID');
                if ~isempty(tgtIx);
                    options.FileID = varargin{tgtIx +1};
                end
                
                tgtIx = input('overwrite');
                if ~isempty(tgtIx);
                    options.Overwrite = true;
                end  
                
                tgtIx = input('tLim');
                if ~isempty(tgtIx)
                    tLim         = varargin{tgtIx +1};
                    options.XLim = [min(tLim) max(tLim)];
                end
                
                tgtIx = input('mzLim');
                if ~isempty(tgtIx)
                    mzLim        = varargin{tgtIx +1};
                    options.YLim = [min(mzLim) max(mzLim)];
                end
                
                tgtIx = input('spikes');
                if ~isempty(tgtIx)
                    spks = varargin{tgtIx +1};
                    if spks == 0
                        options.RemSpks = false;
                    else
                        options.RemSpks = true;
                        options.SpkSz   =  spks;
                    end
                end
            end
            
            
            %% DOMZML2FINNEE
            function [obj, MSn] = domzML2Finnee(obj)
                % DOMZML2FINNEE is THE function that read the mzML file, 
                % organize and store the information.
                
                % 1- Initialisations
                % Open and check the mzML file
                
                fidRead = fopen(obj.FileIn, 'r'); % original mzML file
                if fidRead == -1 && fidWriteDat == -1
                    error('myApp:argChk', ...
                        'Error while opening the files. Type help hyphMSdata2struct for more information')
                end
                MSn = {};
                
                % Read the mzML file until the <run> tag is fund. All 
                % information before this tag are copied to MZMLDump
                obj.MZMLDump = readMZML( fidRead, 'run');
                
                % 2- Find key informations.
                % The first MS scan is read to obtain key informations
                % (i.e. centroid or profile scan, axes labels and units,
                % migration/retantion times. Will create the Axis objects
                keptPl = ftell(fidRead); % position of the first scan
                curLine = fgetl(fidRead);
                [LDR, ~] = getMZMLCamp(curLine, fidRead);
                [~, rndStr] = fileparts(tempname);
                
                isMS1 = true;
                while ~strcmp(LDR.label, '/spectrum')
                    if strcmp(LDR.label, 'spectrumList')
                        [bool, ~, provField] = fieldfind( LDR.attributes, 'count');
                        if bool
                            scanCount = str2double(provField{1});
                        else
                            error('attribute count not present in tag spectrumList');
                        end
                        
                    elseif strcmp(LDR.label, 'spectrum')
                        [bool, ~, ~] = fieldfind( LDR.attributes, 'defaultArrayLength');
                        if bool
                        else
                            error('attribute defaultArrayLength not present in tag spectrum');
                        end
                        
                    elseif strcmp(LDR.label, 'cvParam')
                        [bool, ~, field] = fieldfind( LDR.attributes, 'accession');
                        if bool
                            switch field{1}
                                case 'MS:1000127'
                                    infoDts.Format = 'centroid';
                                    infoDts.Log    = 'CTR=1';
                                    
                                case 'MS:1000128'
                                    infoDts.Format = 'profile';
                                    infoDts.Log    = 'PRF=1';
                                    
                                case 'MS:1000016'
                                    [~, ~, provField] = fieldfind( LDR.attributes, 'unitName');
                                    infoX.Label   = 'Time';
                                    infoX.Unit    = provField{1};
                                    infoX.dp      = 2;
                                    infoX.Loc     = 'None';
                                    AxisX         = Axis(infoX);
                                    infoDts.AxisX = AxisX;
                                    infoPrf.AxisX = AxisX;
                                    
                                case 'MS:1000574'
                                    compression = field{1};
                                    
                                case 'MS:1000514'
                                    [~, ~, provField] = fieldfind( LDR.attributes, 'unitName');
                                    infoY.Label       = 'Mass';
                                    infoY.Unit        = provField{1};
                                    infoY.dp          = 4;
                                    infoY.Loc         = 'None';
                                    AxisY             = Axis(infoY);
                                    infoDts.AxisY     = AxisY;
                                    infoScn.AxisX     = AxisY;
                                    
                                case 'MS:1000523'
                                    dataFormat = field{1};
                                    
                                case 'MS:1000515'
                                    [~, ~, provField] = fieldfind( LDR.attributes, 'unitName');
                                    infoZ.Label       = 'Intensity';
                                    infoZ.Unit        = provField{1};
                                    infoZ.dp          = 0;
                                    infoZ.Loc         = 'None';
                                    AxisZ             = Axis(infoZ);
                                    infoDts.AxisZ     = AxisZ;
                                    infoPrf.AxisY     = AxisZ;
                                    infoScn.AxisY     = AxisZ;
                            end
                        end
                    end
                    curLine = fgetl(fidRead);
                    
                    if ~ischar(curLine)
                        error('myApp:endOfFile', ...
                            'The mzML file is not complete')
                    end
                    [LDR, ~] = getMZMLCamp(curLine, fidRead);
                end
                
                % 3- Process each scan
                fseek(fidRead, keptPl, 'bof'); % Go back at the start of the spectrumList
                curLine  = fgetl(fidRead);
                [LDR, ~] = getMZMLCamp(curLine, fidRead);
                allProfiles         = [];
                infoDts.Title       = 'Original dataset';
                infoDts.ListOfScans = {};
                AxisMZ              = [];
                count               = 1;
                fln                 = 1;
                MZlim               = [inf 0];
                
                h = waitbar(0,'processing scans');
                while count <= scanCount
                    waitbar(count/scanCount)
                    
                    % Find each scan
                    if strcmp(LDR.label, 'spectrum'), boolArray = 1; end
                    
                    if strcmp(LDR.label, '/mzML'), break; end
                    
                    if strcmp(LDR.label, 'cvParam') && ~isempty(LDR.attributes)
                        [~, ~, field] = fieldfind( LDR.attributes, 'accession');
                        if strcmp(field{1}, 'MS:1000016')
                            [~, ~, provField] = fieldfind( LDR.attributes, 'value');
                            tm = str2double(provField{1});
                            if tm >= options.XLim(1) && tm <= options.XLim(2)
                                    doRecord = true;
                                    if isMS1
                                        allProfiles(end+1, 1) = str2double(provField{1}); 
                                    else
                                        MSn{end+1}.Time = str2double(provField{1}); 
                                    end
                            elseif tm > options.XLim(2)
                                break
                            else
                                doRecord = false;
                            end
                            
                        end
                        
                        if strcmp(field{1}, 'MS:1000580') 
                             isMS1 = false;
                        end
                        
                        if strcmp(field{1}, 'MS:1000827') && doRecord
                            [~, ~, provField] = fieldfind( LDR.attributes, 'value');
                            MSn{end}.Precursor(1) =  str2double(provField{1}); 
                        end
                        
                        if strcmp(field{1}, 'MS:1000828') && doRecord
                            [~, ~, provField] = fieldfind( LDR.attributes, 'value');
                            MSn{end}.Precursor(2) =  str2double(provField{1}); %#ok<*AGROW>
                        end
                        
                        if strcmp(field{1}, 'MS:1000829') && doRecord
                            [~, ~, provField] = fieldfind( LDR.attributes, 'value');
                            MSn{end}.Precursor(3) =  str2double(provField{1});
                        end
                    end
                    
                    if strcmp(LDR.label, 'binary')
                        input = LDR.text;
                        if ~isempty(input) && doRecord
                            switch dataFormat
                                case 'MS:1000523'
                                    output = base64decode(input);
                                otherwise
                                    error('precision not recognized')
                            end
                            switch compression
                                case 'MS:1000574'
                                    output = zlibdecode(output);
                                otherwise
                                    error('compression not recognized')
                            end
                            output= typecast(uint8(output),'double');
                        else
                            output = 0;
                        end
                        if boolArray
                            boolArray = 0;
                            mzValue = output';
                        else
                            boolArray = 1;
                            intValue = output';
                            MS =  [mzValue intValue];
                            count = count + 1;
                            
                            if ~isMS1 && doRecord
                                MSn{end}.Data = MS;
                                isMS1 = true;
                            else
                                
                                %% FOR PROFILE MS SCANS
                                if strcmp(infoDts.Format, 'profile') && doRecord
                                    infoDts.Path2Dat{fln} = fullfile(obj.Path2Fin, rndStr);
                                    infoScn.Title     = ['Profile scan #', num2str(count)];
                                    infoScn.FT        = infoDts.Log;
                                    infoScn.TT        = 'PRF';
                                    infoScn.P2Fin     = obj.Path2Fin;
                                    infoScn.Path2Dat  = infoDts.Path2Dat{fln};
                                    infoScn.Precision = 'single';
                                    infoScn.Loc       = 'inFile';
                                    infoScn.AdiPrm    = {};
                                    
                                    % Cut MS to size
                                    if ~isempty(MS)
                                        Id2Rem = MS(:,1) < obj.Options.YLim(1) ...
                                            | MS(:,1) > obj.Options.YLim(2);
                                        MS(Id2Rem, :) = [];
                                    end
                                        
                                    % Filter spikes if needed
                                    if ~isempty(MS) && obj.Options.RemSpks
                                        spkSz = obj.Options.SpkSz;
                                        MS    = spikesRemoval(MS, spkSz );
                                    end
                                    
                                    % reduced trailing zero in excess
                                    if ~isempty(MS)
                                        provMat      = [MS(2:end, 2); 0];
                                        provMat(:,2) = MS(:, 2);
                                        provMat(:,3) = [0; MS(1:end-1, 2)];
                                        MS           = MS(sum(provMat, 2) > 0, :);
                                    end
                                    
                                    % recorded each scans
                                    if isempty(MS)
                                        infoDts.ListOfScans{end+1} = Trace(infoScn);
                                    else
                                        infoDts.ListOfScans{end+1} = Trace(infoScn, MS);
                                        if min(MS(:,1)) == 0
                                            disp('wtf')
                                        end 
                                        MZlim(1) = min(MZlim(1), min(MS(:,1)));
                                        MZlim(2) = max(MZlim(2), max(MS(:,1)));
                                        
                                    end
                                    
                                    s = dir(infoDts.Path2Dat{fln});
                                    if s.bytes > obj.Options.MaxFileSize;
                                        [~, rndStr]           = fileparts(tempname);
                                        fln                   = fln + 1;
                                        infoDts.Path2Dat{fln} = fullfile(obj.Path2Fin, rndStr);
                                        infoScn.Path2Dat      = infoDts.Path2Dat{fln};
                                    end
                                    
                                    %% FOR CENTROID MS SCANS
                                elseif strcmp(infoDts.Format, 'centroid')  && doRecord
                                    
                                    % 2. Do if a scan is centroid mode
                                    infoDts.Path2Dat{fln} = fullfile(obj.Path2Fin, rndStr);
                                    infoScn.Title     = ['Centroid scan #', num2str(count)];
                                    infoScn.FT        = infoDts.Log;
                                    infoScn.TT        = 'CTR';
                                    infoScn.Path2Dat  = infoDts.Path2Dat{fln};
                                    infoScn.Precision = 'single';
                                    infoScn.Loc       = 'inFile';
                                    infoScn.AdiPrm    = {};
                                    infoScn.P2Fin     = obj.Path2Fin;
                                    
                                    % Cut MS to size
                                    if ~isempty(MS)
                                        Id2Rem = MS(:,1) < obj.Options.YLim(1) ...
                                            | MS(:,1) > obj.Options.YLim(2);
                                        MS(Id2Rem, :) = [];
                                    end
                                    
                                    % recorded each scans
                                    if isempty(MS)
                                        infoDts.ListOfScans{end+1} = Trace(infoScn);
                                    else
                                        infoDts.ListOfScans{end+1} = Trace(infoScn, MS);
                                        if min(MS(:,1)) == 0
                                            disp('wtf')
                                        end 
                                        MZlim(1) = min(MZlim(1), min(MS(:,1)));
                                        MZlim(2) = maz(MZlim(2), max(MS(:,1)));
                                        
                                    end
                                    
                                    s = dir(infoDts.Path2Dat{fln});
                                    if s.bytes > obj.Options.MaxFileSize;
                                        [~, rndStr]           = fileparts(tempname);
                                        fln                   = fln + 1;
                                        infoDts.Path2Dat{fln} = fullfile(obj.Path2Fin, rndStr);
                                        infoScn.Path2Dat      = infoDts.Path2Dat{fln};
                                    end
                                end
                                
                                if  doRecord
                                    if  ~isempty(MS)
                                        allProfiles(end, 2) = sum(MS(:,2));
                                        allProfiles(end, 3) = max(MS(:,2));
                                    else
                                        allProfiles(end, 2) = 0;
                                        allProfiles(end, 3) = 0;
                                    end
                                end
                            end
                        end
                        
                    end
                    
                    curLine = fgetl(fidRead);
                    if ~ischar(curLine)
                        warning('myApp:endOfFile', ...
                            'The mzML file is not complete')
                        break
                    end
                    [LDR, ~] = getMZMLCamp(curLine, fidRead);
                end
                if ishandle(h), close(h); end
                fclose(fidRead);
                
                while allProfiles(end, 1) == 0
                    allProfiles(end, :) = [];
                end
                
                % 4- Create all axes and profiles and save data
                infoAxis           = AxisX.InfoAxis; %Axis time
                infoAxis.Loc       = 'inFile';
                infoAxis.Precision = 'single';
                infoAxis.Path2Dat  = infoDts.Path2Dat{fln};
                if isempty(allProfiles)
                    infoDts.AxisX  = Axis(infoAxis);
                else
                    infoDts.AxisX  = Axis(infoAxis, allProfiles(:,1));
                end
                
                infoAxis           = AxisY.InfoAxis; %Axis m/z
                infoAxis.Loc       = 'inFile';
                infoAxis.Precision = 'single';
                infoAxis.Path2Dat  = infoDts.Path2Dat{fln};
                if isempty(AxisMZ)
                    infoDts.AxisY  = Axis(infoAxis);
                else
                    infoDts.AxisY  = Axis(infoAxis, AxisMZ(:,1));
                end
                
                infoPrf.AxisX      = AxisX;
                infoPrf.AxisY      = AxisZ;
                infoPrf.Loc       = 'inFile';
                infoPrf.Precision = 'single';
                infoPrf.Path2Dat  = infoDts.Path2Dat{fln};
                infoPrf.P2Fin     = obj.Path2Fin;
                infoPrf.FT        = infoDts.Log;
                infoPrf.TT        = 'SEP';
                infoPrf.AdiPrm    = {};
                
                % *Base Peak Profile (BPP)*: most intensed ions in each MS
                % scan.
                infoPrf.Title     = 'Base Peak Profile';
                if isempty(allProfiles)
                    infoDts.BPP   = Trace(infoPrf);
                else
                    infoDts.BPP   = Trace(infoPrf, [allProfiles(:,1), allProfiles(:,3)]);
                end
                
                % *Total Ion Profile (TIP)*: sum of all ions in each scan
                infoPrf.Title     = 'Total Ion profile';
                if isempty(allProfiles)
                    infoDts.TIP   = Trace(infoPrf);
                else
                    infoDts.TIP   = Trace(infoPrf, [allProfiles(:,1), allProfiles(:,2)]);
                end
                
                % *Total Ion Spectrum (TIS)*, empty with centroid data. With
                % profile data, the TIS is an MS spectrum that is the sum
                % of all MS spectra
                infoPrf.AxisX      = AxisY;
                infoPrf.Title     = 'Total Ion Spectrum';
                if isempty(AxisMZ)
                    infoDts.TIS   = Trace(infoPrf);
                else
                    infoDts.TIS   = Trace(infoPrf, [AxisMZ(:,1), AxisMZ(:,2)]);
                end
                
                % *Frequency Ion Spectrum (FIS)*, empty with centroid data.
                % With profile data, the FIS is the recording for each m/z
                % increment of the number of non-zeros intensities values
                % along all scans
                infoPrf.Title     = 'Frequency Ion Spectrum';
                if isempty(AxisMZ)
                    infoDts.FIS   = Trace(infoPrf);
                else
                    infoDts.FIS   = Trace(infoPrf, [AxisMZ(:,1), AxisMZ(:,3)]);
                end
                
                % *Base Ion Spectrum (BIS)*, empty with centroid data.
                % With profile data, the BIS is the most intense ions
                % measured in every scans at each m/z increment
                infoPrf.Title     = 'Base Ion Spectrum';
                if isempty(AxisMZ)
                    infoDts.BIS   = Trace(infoPrf);
                else
                    infoDts.BIS   = Trace(infoPrf, [AxisMZ(:,1), AxisMZ(:,4)]);
                end
                
                infoDts.LAST      = Trace();
                
                % save and go
                infoDts.Option4crt.function = 'domzML2Finnee.m';
                infoDts.P2F        = obj.Path2Fin;
                infoDts.MZlim      = MZlim;
                obj.Datasets{end+1}= Dataset(infoDts);

            end
        end
        
        function save(obj)
            %% DESCRIPTION
            % SAVE allows to save the finnee object.
            %
            %% EXAMPLES
            % myFinnee.save will save the last changes. This is normally
            % done automaticall with most methods.
            %
            myFinnee = obj;
            save(fullfile(obj.Path2Fin, 'myFinnee.mat'), 'myFinnee')
        end
    end
end
