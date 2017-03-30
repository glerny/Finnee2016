%% DESCRIPTION
% FINNEE is the main class of the Finnee2016 toolbox. Finnee is used to
% record the information about the sample and separation; to maintain the
% links with the various binary files (often name dat files) necessary to
% store all the experimental data; and will contain multiple datasets to
% follow each transformation. Finnee2016 was engineered for data recorded
% using separation techniques hyphenated with high-resolution mass
% spectrometry (X-HRMS
% <https://github.com/glerny/Finnee2016/wiki/Definition-of-terms> ). The
% Finnee class contains a list of properties to describe and organise the
% information and specific method to explore, display or transform those
% data.
%
%% LIST OF THE CLASS'S PROPERTIES
% *FileId*        : The generic name
% DateOfCreation:
% *FileIn*        : The original mzML file
% *Datasets*      : The list of all datasets that are linked with this Finnee
%   object. Initially, only one dataset will exist, the original one, but
%   new datasets can easily be created.
% *Options*       : (Hidden) Options used for the creation of the Finnee
%   object.
% *Path2Fin*      : (Hidden) The path to the parent directory
%   (directory with the .fin extension)
% *MZMLDump*      : (Hidden) A copy of the information that was present in
%   the original mzML file.
%
%% LIST OF THE CLASS'S METHODS
% *Finnee*            : The constructor method.
% *queryMZMLDump*     : not implemented yet
% *Save*              : record changes
% *BaselineCorrection*: (Profile dataset only). Allow correcting profiles
%   in a dataset from baseline drift. BaselineCorrection will create a new
%   dataset.
% *DoCentroid*        : (Profile dataset only). Will transform each profile
%   MS scan in a dataset to centroid MS scan. DoCentroid will create a new
%   dataset.
% *FilterDataset*     : Allow to filter a dataset using different method.
%   FilterDataset will create a new dataset. The list of filter are
%   * _RemoveSpikes_
%   * _RemoveNoise_
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
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Finnee
    
    properties
        FileID         % The generic name
        DateOfCreation % The date of creation
        FileIn         % Name and path of the original mzML file
        Datasets       % The array of datasets linked to this object
    end
    
    properties (Hidden = true, SetAccess = immutable)
        Options        % The options that were used by the constructor 
                       % method
        Path2Fin       % The path to the .fin folder
    end
    
    properties (Hidden = true)
        MZMLDump       % Some general information from the mzML file
    end
    
    methods
        function obj = Finnee(varargin)
            %% DESCRIPTION
            % Finnee is the constructor method and can be used with
            % various optional parameters
            % * _fileIn_ followed by the full path and name of the mzML 
            %   file.
            % * _folderOut_ followed by the path to the destination folder. 
            % * _fileID_ followed by a generic name.
            % * _overwrite_ delete the folder at folderOut\fileID.fin if it 
            %   exists
            % * _ tLim_ not implemented
            % * _ mzLim_ not implemented
            % * _ Mztolerance_ used to normalised m/z axis of every scans 
            %   to the first axis. The default value is 3.
            % * _ spikes_ (see the method @Finnee\FilterDataset) for 
            %   additional information. Remove spikes in every MS scan. 
            %   If used, the options should be followed by 0,1, 2 ot 3 to 
            %   define the size of spikes to be considered for removal. The
            %   default value is 2; 0 allows to to turn off spikes removal.
            %
            %% EXAMPLES
            % * myFinnee = Finnee; 
            %   Will create a Finnee object without any options. The target 
            %   mzML file, destination folder and generic name will be
            %   asked during the creation of the object
            % *	myFinnee = Finnee('fileIn'   , 'K:\Data\test.mzML', ...
            %                     'folderOut', 'K:Finnee',          ...
            %                     'Overwrite',                      ...
            %                     'fileID'   , 'test1',             ...
            %                     'Spikes'   , 0);
            %   Here all parameters are defined by the optional parameters
            
            %% CORE OF THE FUNCTION
            % 1. Initialisation and options
            
            % Check the options and create the Finnee object
            options = checkVarargin(varargin{:});
            
            % Initialize Finnee
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
            
            % Run the domzML2Finnee to read the mzML file and parse the
            % information to the Finnee object
            obj      = domzML2Finnee(obj);
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
                options.ID4Saving   = 'myFinnee';
                options.MZtol       = 3;
                options.Overwrite   = false;
                options.FileFormat  = 'mzML';
                options.MaxFileSize = 10000000000;
                options.XLim        = [0 inf];
                options.YLim        = [0 inf];
                options.RemSpks     = true;
                options.SpkSz       = 2;
                
                % 2- Decipher varargin and update options when relevamt
                input = @(x) find(strcmpi(varargin,x),1);
                
                tgtIx = input('tLim');
                if ~isempty(tgtIx)
                    tLim         = varargin{tgtIx +1};
                    options.XLim = [min(tLim) max(tLim)];
                end
                
                tgtIx = input('mzLim');
                if ~isempty(tgtIx)
                    mzLim          = varargin{tgtIx +1};
                    options.Ylim(1)= [min(mzLim) max(mzLim)];
                end
                
                tgtIx  = input('MZtolerance');
                if ~isempty(tgtIx)
                    options.MZtol = varargin{tgtIx +1};
                end
                
                tgtIx = input('spikes');
                if ~isempty(tgtIx)
                    spks = varargin{tgtIx +1};
                    if spks == 0
                        options.RemSpks = false;
                    else
                        options.RemSpks = true;
                        options.SpkSz  =  spks;
                    end
                end
                
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
            end
            
                        
            %% DOMZML2FINNEE
            function obj = domzML2Finnee(obj)
                % DOMZML2FINNEE is THE function that will read the mzML
                % file, organize and store the information in the Finnee
                % object and .dat file.
                
                % 1- Initialisations
                % Open and check the mzML file
                fidRead = fopen(obj.FileIn, 'r'); % original mzML file
                if fidRead == -1 && fidWriteDat == -1
                    error('myApp:argChk', ...
                        'Error while opening the files. Type help hyphMSdata2struct for more information')
                end
                
                % Read the mzML file until the <run> tag is fund. This tag
                % indicate that a MS can has been recorded. All information
                % before this tag are copied to MZMLDump
                obj.MZMLDump = readMZML( fidRead, 'run');
                
                % 2- Find key informations.
                % The first MS scan is read to obtain key informations
                % (i.e. centroid or profile scan, axes label and unit,
                % migration/retantion times. Will create the Axis objects
                keptPl = ftell(fidRead); % position of the first scan
                curLine = fgetl(fidRead);
                [LDR, ~] = getMZMLCamp(curLine, fidRead);
                [~, rndStr] = fileparts(tempname);
                
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
                                    infoX.Label  = 'Time';
                                    infoX.Unit   = provField{1};
                                    infoX.dp     = 2;
                                    infoX.Loc    = 'None';
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
                                    AxisY              = Axis(infoY);
                                    infoDts.AxisY      = AxisY;
                                    infoScn.AxisX      = AxisY;
                                    
                                case 'MS:1000523'
                                    dataFormat = field{1};
                                    
                                case 'MS:1000515'
                                    [~, ~, provField] = fieldfind( LDR.attributes, 'unitName');
                                    infoZ.Label       = 'Intensity';
                                    infoZ.Unit        = provField{1};
                                    infoZ.dp          = 0;
                                    infoZ.Loc         = 'None';
                                    AxisZ              = Axis(infoZ);
                                    infoDts.AxisZ      = AxisZ;
                                    infoPrf.AxisY      = AxisZ;
                                    infoScn.AxisY      = AxisZ;
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
                
                allProfiles         = zeros(scanCount, 3);
                infoDts.Title       = 'Original dataset';
                infoDts.ListOfScans = {};
                AxisMZ               = [];
                count               = 1;
                fln                 = 1;
                
                h = waitbar(0,'processing scans');
                while count <= scanCount
                    waitbar(count/scanCount)
                    
                    % Find each scan
                    if strcmp(LDR.label, 'spectrum'), boolArray = 1; end
                    
                    if strcmp(LDR.label, 'cvParam') && ~isempty(LDR.attributes)
                        [~, ~, field] = fieldfind( LDR.attributes, 'accession');
                        if strcmp(field{1}, 'MS:1000016')
                            [~, ~, provField] = fieldfind( LDR.attributes, 'value');
                            allProfiles(count, 1) = str2double(provField{1});
                        end
                    end
                    
                    if strcmp(LDR.label, 'binary')
                        input = LDR.text;
                        if ~isempty(input)
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
                            
                            %% FOR PROFILE MS SCANS
                            if strcmp(infoDts.Format, 'profile')
                                infoDts.Path2Dat{fln} = fullfile(obj.Path2Fin, rndStr);
                                infoDts.tol4MZ    = obj.Options.MZtol;
                                infoScn.Title     = ['Profile scan #', num2str(count)];
                                infoScn.FT        = infoDts.Log;
                                infoScn.TT        = 'PRF';
                                infoScn.Path2Dat  = infoDts.Path2Dat{fln};
                                infoScn.Precision = 'single';
                                infoScn.Loc       = 'inFile';
                                infoScn.AdiPrm{1}.type = 'AxisMZ normalisation';
                                infoScn.AdiPrm{1}.var  = 0;
                                
                                % sort out the whole mz axes and normalised scan from
                                % each scan
                                indNotZeros = intValue ~= 0;
                                if isempty(AxisMZ)
                                    AxisMZ = mzValue;
                                    AxisMZ(:,2) = intValue;
                                    AxisMZ(indNotZeros,3) = 1;
                                    AxisMZ(:,4) = intValue;
                                else
                                    
                                    tol = obj.Options.MZtol;
                                    
                                    % Normalised the new mzAxis to the Axis from the first
                                    % scan
                                    [~,ia,ib] = intersect(round(AxisMZ(:,1)*10^tol)/10^tol,...
                                        round(mzValue*10^tol)/10^tol);
                                    cst = mean((AxisMZ(ia, 1) - mzValue(ib))./AxisMZ(ia));
                                    infoScn.AdiPrm{1}.type = 'AxisMZ normalisation';
                                    infoScn.AdiPrm{1}.var  = cst;
                                    
                                    % Check if the two axes contain different values
                                    axis = mzValue/(1-cst);
                                    [tf, loc] = ismember(round(axis*10^tol)/10^tol, ...
                                        round(AxisMZ(:,1)*10^tol)/10^tol);
                                    indZeros = find(tf == 0);
                                    if ~isempty(indZeros)
                                        [tfc, locc] = ismember(ceil(axis*10^tol)/10^tol, ...
                                            ceil(AxisMZ(:,1)*10^tol)/10^tol);
                                        loc(indZeros) = locc(indZeros);
                                        tf(indZeros) = tfc(indZeros);
                                    end
                                    
                                    % Update AxisMZ
                                    if any(tf)
                                        AxisMZ(loc(tf), 2) = ...
                                            AxisMZ(loc(tf), 2) + intValue(loc ~=0); %#ok<*AGROW>
                                        freqAtt = intValue(loc ~=0) > 0;
                                        AxisMZ(loc(tf), 3) = ...
                                            AxisMZ(loc(tf), 3) + freqAtt;
                                        AxisMZ(loc(tf), 4) = max([AxisMZ(loc(tf), 4), ...
                                            intValue(loc ~=0)], [], 2);
                                    end
                                    
                                    if any(~tf)
                                        MZ2add = axis(~tf);
                                        MZ2add(:,2) = intValue(~tf);
                                        indNotZeros = MZ2add(:,2) ~= 0;
                                        MZ2add(indNotZeros,3) = 1 ;
                                        MZ2add(:,4) = intValue(~tf);
                                        AxisMZ = [AxisMZ; MZ2add];
                                        AxisMZ = sortrows(AxisMZ,1);
                                    end
                                end
                                
                                % Filter spikes if needed
                                if obj.Options.RemSpks
                                    spkSz = obj.Options.SpkSz;
                                    MS    = spikesRemoval(MS, spkSz );
                                end
                                
                                % reduced trailing zero in excess
                                provMat      = [MS(2:end, 2); 0];
                                provMat(:,2) = MS(:, 2);
                                provMat(:,3) = [0; MS(1:end-1, 2)];
                                MS           = MS(sum(provMat, 2) > 0, :);
                                
                                % recorded each scans
                                if isempty(MS)
                                    infoDts.ListOfScans{count} = Trace(infoScn);
                                else
                                    infoDts.ListOfScans{count} = Trace(infoScn, MS);
                                end
                                
                                s = dir(infoDts.Path2Dat{fln});
                                if s.bytes > obj.Options.MaxFileSize;
                                    [~, rndStr]           = fileparts(tempname);
                                    fln                   = fln + 1;
                                    infoDts.Path2Dat{fln} = fullfile(obj.Path2Fin, rndStr);
                                    infoScn.Path2Dat      = infoDts.Path2Dat{fln};
                                end
                                
                            
                            %% FOR CENTROID MS SCANS
                            elseif strcmp(infoDts.Format, 'centroid')
                                % 2. Do if a scan is centroid mode
                                infoDts.Path2Dat{fln} = fullfile(obj.Path2Fin, rndStr);
                                infoDts.tol4MZ    = NaN;
                                
                                infoScn.Title     = ['Centroid scan #', num2str(count)];
                                infoScn.FT        = infoDts.Log;
                                infoScn.TT        = 'CTR';
                                infoScn.Path2Dat  = infoDts.Path2Dat{fln};
                                infoScn.Precision = 'single';
                                infoScn.Loc       = 'inFile';
                                infoScn.AdiPrm    = {};
                                
                                % recorded each scans
                                if isempty(MS)
                                    infoDts.ListOfScans{count} = Trace(infoScn);
                                else
                                    infoDts.ListOfScans{count} = Trace(infoScn, MS);
                                end
                                
                                s = dir(infoDts.Path2Dat{fln});
                                if s.bytes > obj.Options.MaxFileSize;
                                    [~, rndStr]           = fileparts(tempname);
                                    fln                   = fln + 1;
                                    infoDts.Path2Dat{fln} = fullfile(obj.Path2Fin, rndStr);
                                    infoScn.Path2Dat      = infoDts.Path2Dat{fln};
                                end
                                
                            else
                                error('The data type was not recognised')
                            end
                            
                            allProfiles(count, 2) = sum(MS(:,2));
                            allProfiles(count, 3) = max(MS(:,2));
                            count = count + 1;
                        end
                    end
                    
                    curLine = fgetl(fidRead);
                    if ~ischar(curLine)
                        error('myApp:endOfFile', ...
                            'The mzML file is not complete')
                    end
                    [LDR, ~] = getMZMLCamp(curLine, fidRead);
                end
                try close(h), catch, end
                fclose(fidRead);
                
                % 4- Create all azxes and profiles and save data
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
