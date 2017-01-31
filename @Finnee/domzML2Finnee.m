%% Description
% DOMZML2FINNEE takes as entry paramter a options4finnee object and will
% create a Finnee object. DOMZML2FINNEE will load, read and convert an mzML
% file to a Finnee object.
%
%% Copyright 
% BSD 3-Clause License
% Copyright 2016-2017 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function obj = domzML2Finnee(obj)

%% CORE OF THE FUNCTION
% 1. INITIALISATION

% INITFUNCTION is used to test the entries, load the target MS dataset and
% load the default values
obj.Path2Fin = obj.Options.Path2Fin;
obj.FileID   = obj.Options.FileID;

fidRead = fopen(obj.Options.FileIn, 'r'); % original mzML file
if fidRead == -1 && fidWriteDat == -1
    error('myApp:argChk', ...
        'Error while opening the files. Type help hyphMSdata2struct for more information')
end

% Read and store mzML information up to <run ....
obj.MZMLDump = readMZML( fidRead, 'run');

% Get axes Unit and Labels
% i. Initialisation and default values
keptPl = ftell(fidRead);
curLine = fgetl(fidRead);
[LDR, ~] = getMZMLCamp(curLine, fidRead);
obj.Datasets{1} = Dataset;
obj.Datasets{1}.Path2Fin = obj.Options.Path2Fin;
[~, rndStr] = fileparts(tempname);

% ii. get axes labels, units, and MS scan type (centroid or profile)
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
                    obj.Datasets{1}.Format = 'MS centroid';
                case 'MS:1000128'
                    obj.Datasets{1}.Format = 'MS profile';
                case 'MS:1000016'
                    [~, ~, provField] = fieldfind( LDR.attributes, 'unitName');
                    obj.Datasets{1}.XUnit = provField{1};
                case 'MS:1000574'
                    compression = field{1};
                case 'MS:1000514'
                    [~, ~, provField] = fieldfind( LDR.attributes, 'unitName');
                    obj.Datasets{1}.YUnit = provField{1};
                case 'MS:1000523'
                    dataFormat = field{1};
                case 'MS:1000515'
                    [~, ~, provField] = fieldfind( LDR.attributes, 'unitName');
                    obj.Datasets{1}.ZUnit = provField{1};
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

fseek(fidRead, keptPl, 'bof'); % Go back at the start of the spectrumList
curLine = fgetl(fidRead);
[LDR, ~] = getMZMLCamp(curLine, fidRead);

allProfiles = zeros(scanCount, 3);
rdg = obj.Options.Rounding;
obj.Datasets{1}.RndFactor = rdg;
obj.Datasets{1}.Title = 'Original dataset';
obj.Datasets{1}.Path2Fin = obj.Path2Fin;
obj.Datasets{1}.Log = 'MSPROFDTS DTS=1';
axeMZ = [];
count = 1;
fln = 1;


h = waitbar(0,'processing scans');
while count <= scanCount
% iii. Find each scan
    if strcmp(LDR.label, 'spectrum'), boolArray = 1; end
    
    if strcmp(LDR.label, 'cvParam') && ~isempty(LDR.attributes)
        [~, ~, field] = fieldfind( LDR.attributes, 'accession');
        if strcmp(field{1}, 'MS:1000016')
            [~, ~, provField] = fieldfind( LDR.attributes, 'value');
            allProfiles(count, 1) = str2double(provField{1});
            
            waitbar(count/scanCount)
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
            
            if strcmp(obj.Datasets{1}.Format, 'MS profile')
                % 1. Do if a scan is profile mode
                obj.Datasets{1}.Path2Dat{fln} = fullfile(obj.Path2Fin, rndStr) ;

                structInfo.title        = ['Profile scan #', num2str(count)];
                structInfo.traceType    = obj.Datasets{1}.Format;
                structInfo.XLabel       = obj.Datasets{1}.YLabel;
                structInfo.XUnit        = obj.Datasets{1}.YUnit;
                structInfo.YLabel       = obj.Datasets{1}.ZLabel;
                structInfo.YUnit        = obj.Datasets{1}.ZUnit;
                structInfo.Path2Dat     = obj.Datasets{1}.Path2Dat{fln} ;
                structInfo.Variables    = 0;
                structInfo.Precision    = 'single';
                structInfo.Path2Fin     = obj.Path2Fin;
                structInfo.Log          = 'PRFMSSCN DTS=1';
                
                % sort out the whole mz axes and normalised scan from
                % each scan
                indNotZeros = intValue ~= 0;
                
                if isempty(axeMZ)
                    axeMZ = mzValue;
                    axeMZ(:,2) = intValue;
                    axeMZ(indNotZeros,3) = 1;
                    axeMZ(:,4) = intValue;
                else
                    % Normalize axeMS with new axe and check for
                    % missing values
                    [~,ia,ib] = intersect(round(axeMZ(:,1)*10^rdg)/10^rdg,...
                        round(mzValue*10^rdg)/10^rdg);
                    cst = mean((axeMZ(ia) - mzValue(ib))./axeMZ(ia));
                    structInfo.Variables = cst;
                    axe = mzValue/(1-cst);
                    [tf, loc] = ismember(round(axe*10^rdg)/10^rdg, ...
                        round(axeMZ(:,1)*10^rdg)/10^rdg);
                    indZeros = find(tf == 0);
                    if ~isempty(indZeros)
                        [tfc, locc] = ismember(ceil(axe*10^rdg)/10^rdg, ...
                            ceil(axeMZ(:,1)*10^rdg)/10^rdg);
                        loc(indZeros) = locc(indZeros);
                        tf(indZeros) = tfc(indZeros);
                    end
                    
                    if ~isempty(find(tf))
                        axeMZ(loc(tf), 2) = ...
                            axeMZ(loc(tf), 2) + intValue(loc ~=0);
                        freqAtt = intValue(loc ~=0) > 0;
                        axeMZ(loc(tf), 3) = ...
                            axeMZ(loc(tf), 3) + freqAtt;
                         axeMZ(loc(tf), 4) = max([axeMZ(loc(tf), 4), ...
                             intValue(loc ~=0)], [], 2);
                    end
                    
                    if  ~isempty(find(~tf))
                        MZ2add = axe(~tf);
                        MZ2add(:,2) = intValue(~tf);
                        indNotZeros = MZ2add(:,2) ~= 0;
                        MZ2add(indNotZeros,3) = 1 ;
                        MZ2add(:,4) = intValue(~tf);
                        axeMZ = [axeMZ; MZ2add];
                        axeMZ = sortrows(axeMZ,1);
                    end
                end
                
                % Filter spikes if needed
                if obj.Options.RemSpks
                end
                
                % reduced MS data by keeping only one trailing zero
                provMat = [MS(2:end, 2); 0];
                provMat(:,2) = MS(:, 2);
                provMat(:,3) = [0; MS(1:end-1, 2)];
                MS = MS(sum(provMat, 2) > 0, :);
                
                obj.Datasets{1}.ListOfScans{count} = Trace(structInfo, MS);
                s = dir(structInfo.Path2Dat);
                if s.bytes > obj.Options.MaxFileSize;
                    [~, rndStr] = fileparts(tempname);
                    fln = fln + 1;
                    obj.Datasets{1}.Path2Dat{fln}= fullfile(obj.Path2Fin, rndStr) ;
                    structInfo.Path2Dat     = obj.Datasets{1}.Path2Dat{fln} ;
                end
                    
                
            elseif strcmp(obj.Datasets{1}.Format, 'MS centroid')
                % 2. Do if a scan is centroid mode
                
                structInfo.title        = ['Centroid scan #', num2str(count)];
                structInfo.traceType    = obj.Datasets{1}.Format;
                structInfo.XLabel       = obj.Datasets{1}.YLabel;
                structInfo.XUnit        = obj.Datasets{1}.YUnit;
                structInfo.YLabel       = obj.Datasets{1}.ZLabel;
                structInfo.YUnit        = obj.Datasets{1}.ZUnit;
                structInfo.link2file    = path2dat;
                structInfo.Variables    = 0;
                structInfo.Precision    = 'single';
                structInfo.Path2Finnee  = obj.Path;
                structInfo.Log          = 'CTRMSSCN DTS=1';
                obj.Datasets{1}.ListOfScans{count} = Trace(structInfo, MS);
                axeMZ = [0 0 0 ];
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
if exist('h', 'var'), close(h); end

fclose(fidRead);

% iv. Ending with scan list, creating Finnee and Dataset

structInfo.XLabel     	= obj.Datasets{1}.XLabel;
structInfo.XUnit      	= obj.Datasets{1}.XUnit;
structInfo.Log      	= 'AXETIME DTS=1';
obj.Datasets{1}.TimeAxe	= Axe(structInfo, allProfiles(:,1));

structInfo.XLabel       = obj.Datasets{1}.YLabel;
structInfo.XUnit        = obj.Datasets{1}.YUnit;
structInfo.Log          = 'AXEMZ DTS=1';
obj.Datasets{1}.MzAxe 	= Axe(structInfo, axeMZ(:,1));

structInfo.title    	= 'Base Peak Profiles';
structInfo.traceType 	= 'ion profile';
structInfo.Log       	= 'BPP DTS=1';
structInfo.XLabel       = obj.Datasets{1}.XLabel;
structInfo.XUnit        = obj.Datasets{1}.XUnit;
structInfo.YLabel   	= obj.Datasets{1}.ZLabel;
structInfo.YUnit        = obj.Datasets{1}.ZUnit;
obj.Datasets{1}.BPP  	= Trace(structInfo, [allProfiles(:,1), allProfiles(:,3)]);

structInfo.title        = 'Total Ion profiles';
structInfo.Log       	= 'TIP DTS=1';
obj.Datasets{1}.TIP  	= Trace(structInfo, [allProfiles(:,1), allProfiles(:,2)]);

structInfo.title        = 'Total Ion Scan';
structInfo.Log          = 'TIS DTS=1';
structInfo.traceType 	= 'MS profile';
structInfo.XLabel   	= obj.Datasets{1}.YLabel;
structInfo.XUnit    	= obj.Datasets{1}.YUnit;
structInfo.YLabel     	= obj.Datasets{1}.ZLabel;
structInfo.YUnit        = obj.Datasets{1}.ZUnit;
obj.Datasets{1}.TIS   	= Trace(structInfo, [axeMZ(:,1), axeMZ(:,2)]);

structInfo.title     	= 'Frequency Ion Scan';
structInfo.Log          = 'FIS DTS=1';
structInfo.traceType 	= 'MS profile';
obj.Datasets{1}.FIS 	= Trace(structInfo, [axeMZ(:,1), axeMZ(:,3)]);

structInfo.title     	= 'Base Ion Scan';
structInfo.Log       	= 'BIS DTS=1';
structInfo.traceType   	= 'MS profile';
obj.Datasets{1}.BIS  	= Trace(structInfo, [axeMZ(:,1), axeMZ(:,4)]);
end