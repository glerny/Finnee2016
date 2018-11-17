%% DESCRIPTION
%
%% INPUT PARAMETERS
%
%% Copyright
% BSD 3-Clause License
% Copyright 2016-2017 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function exportAs( obj, dts, varargin)
%% CORE OF THE FUNCTION

% 1- Initialisation and options

% Check the options and create the Finnee object
options = checkVarargin(varargin{:});
filterIndex = 1;

if isempty(options.FolderOut)
    options.FolderOut = uigetdir(obj.Path2Fin);
end

if isempty (options.FileName)
    txtStg              = 'Export Dataset';
    ext                 = '.mzML';
    path                = options.FolderOut;
    name                = obj.FileID;
    propal              = fullfile(path, [name, '_dataset', num2str(dts), ext]);
    
    [options.FileName, options.FolderOut, filterIndex] = uiputfile(ext, txtStg, propal);
    if ~ischar((options.FileName)) && ~ischar(options.FolderOut)
        error('myApp:argChk', 'User cancel file selection');
    end
end

switch filterIndex
    case 1
        ExportAsmzML
end

%% CHECKVARARGIN

    function ExportAsmzML
        fidRead  = fopen(obj.FileIn, 'r');
        FO = fullfile(options.FolderOut, options.FileName);
        fidWrite = fopen(FO, 'w');
        
        for ii = 1:length(obj.MZMLDump)
            line = obj.MZMLDump{ii};
            % line        = fgetl(fidRead);
            fprintf(fidWrite, '%s', line);
            fprintf(fidWrite, '\n');
        end
        line = '<spectrumList count="%s" defaultDataProcessingRef="exportation">';
        fprintf(fidWrite, line, num2str(length(obj.Datasets{dts}.ListOfScans)));
        fprintf(fidWrite, '\n');
        
        for ii = 1:length(obj.Datasets{dts}.ListOfScans)
            MS = obj.Datasets{dts}.ListOfScans{ii}.Data;
            
            if isempty(MS)
                MS = [0 0];
            end
            mzArray = MS(:,1);
            mzArray = typecast(mzArray,'uint8');
            mzArray = zlibencode(mzArray);
            mzArray =  base64encode(mzArray);
            ItArray = MS(:,2);
            ItArray = typecast(ItArray,'uint8');
            ItArray = zlibencode(ItArray);
            ItArray = base64encode(ItArray);
            
            line = '<spectrum index="%s" id="scan=%s" defaultArrayLength="%s">';
            fprintf(fidWrite, line, num2str(ii-1), num2str(ii),  num2str(size(MS,1)));
            fprintf(fidWrite, '\n');
            
            line = '<cvParam cvRef="MS" accession="MS:1000130" name="positive scan"/>'; %TO RECORD
            fprintf(fidWrite, line);
            fprintf(fidWrite, '\n');
            
            line = '<cvParam cvRef="MS" accession="MS:1000579" name="MS1 spectrum"/>'; %TO RECORD
            fprintf(fidWrite, line);
            fprintf(fidWrite, '\n');
            
            switch obj.Datasets{dts}.Format
                case 'profile'
                    line = '<cvParam cvRef="MS" accession="MS:1000128" name="profile spectrum"/>'; %TO RECORD
                case 'centroid'
                    line = '<cvParam cvRef="MS" accession="MS:1000128" name="centroid spectrum"/>'; %TO RECORD
            end
            fprintf(fidWrite, line);
            fprintf(fidWrite, '\n');
            
            line = '<cvParam cvRef="MS" accession="MS:1000285" name="total ion current" value="%s"/>';
            fprintf(fidWrite, line, num2str(sum(MS(:,2))));
            fprintf(fidWrite, '\n');
            
            line = '<cvParam cvRef="MS" accession="MS:1000504" name="base peak m/z" value="%s" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>';
            fprintf(fidWrite, line, num2str(max(MS(:,2))));
            fprintf(fidWrite, '\n');
            
            line = '<cvParam cvRef="MS" accession="MS:1000505" name="base peak intensity" value="0.00000000" unitCvRef="MS" unitAccession="MS:1000131" unitName="number of detector counts"/>';
            fprintf(fidWrite, line);
            fprintf(fidWrite, '\n');
            
            line = '<cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="1"/>';
            fprintf(fidWrite, line);
            fprintf(fidWrite, '\n');
            
            line = '<cvParam cvRef="MS" accession="MS:1000527" name="highest observed m/z" value="%s" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>';
            fprintf(fidWrite, line, num2str(max(MS(:,1))));
            fprintf(fidWrite, '\n');
            
            line = '<cvParam cvRef="MS" accession="MS:1000528" name="lowest observed m/z" value="%s" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>';
            fprintf(fidWrite, line, num2str(min(MS(:,1))));
            fprintf(fidWrite, '\n');
            
            line = '<scanList count="1">';
            fprintf(fidWrite, line);
            fprintf(fidWrite, '\n');
            
            line = '<cvParam cvRef="MS" accession="MS:1000795" name="no combination"/>';
            fprintf(fidWrite, line);
            fprintf(fidWrite, '\n');
            
            line = '<scan instrumentConfigurationRef="instrument">';
            fprintf(fidWrite, line);
            fprintf(fidWrite, '\n');
            
            line = '<cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="%s" unitCvRef="UO" unitAccession="UO:0000031" unitName="minute"/>';
            fprintf(fidWrite, line, num2str(obj.Datasets{dts}.AxisX.Data(ii)));
            fprintf(fidWrite, '\n');
            
            line = '</scan>';
            fprintf(fidWrite, line);
            fprintf(fidWrite, '\n');
            
            line = '</scanList>';
            fprintf(fidWrite, line);
            fprintf(fidWrite, '\n');
            
            line = '<binaryDataArrayList count="2">';
            fprintf(fidWrite, line);
            fprintf(fidWrite, '\n');
            
            line = '<binaryDataArray encodedLength="%s">';
            fprintf(fidWrite, line, num2str(size(mzArray,2)));
            fprintf(fidWrite, '\n');
            
            line = '<cvParam cvRef="MS" accession="MS:1000574" name="zlib compression"/>';
            fprintf(fidWrite, line);
            fprintf(fidWrite, '\n');
            
            line = '<cvParam cvRef="MS" accession="MS:1000514" name="m/z array" value="" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>';
            fprintf(fidWrite, line);
            fprintf(fidWrite, '\n');
            
            line = '<cvParam cvRef="MS" accession="MS:1000523" name="64-bit float"/>';
            fprintf(fidWrite, line);
            fprintf(fidWrite, '\n');
            
            line = '<binary>%s</binary>';
            fprintf(fidWrite, line, mzArray);
            fprintf(fidWrite, '\n');
            
            line = '</binaryDataArray>';
            fprintf(fidWrite, line);
            fprintf(fidWrite, '\n');
            
            line = '<binaryDataArray encodedLength="%s">';
            fprintf(fidWrite, line, num2str(size(ItArray,2)));
            fprintf(fidWrite, '\n');
            
            line = '<cvParam cvRef="MS" accession="MS:1000574" name="zlib compression"/>';
            fprintf(fidWrite, line);
            fprintf(fidWrite, '\n');
            
            line = '<cvParam cvRef="MS" accession="MS:1000515" name="intensity array" value="" unitCvRef="MS" unitAccession="MS:1000131" unitName="number of detector counts"/>';
            fprintf(fidWrite, line);
            fprintf(fidWrite, '\n');
            
            line = '<cvParam cvRef="MS" accession="MS:1000523" name="64-bit float"/>';
            fprintf(fidWrite, line);
            fprintf(fidWrite, '\n');
            
            line = '<binary>%s</binary>';
            fprintf(fidWrite, line, ItArray);
            fprintf(fidWrite, '\n');
            
            line = '</binaryDataArray>';
            fprintf(fidWrite, line);
            fprintf(fidWrite, '\n');
            
            line = '</binaryDataArrayList>';
            fprintf(fidWrite, line);
            fprintf(fidWrite, '\n');
            
            line = '</spectrum>';
            fprintf(fidWrite, line);
            fprintf(fidWrite, '\n');
        end
        fprintf(fidWrite, '</spectrumList>');
        fprintf(fidWrite, '\n');
        fprintf(fidWrite, '</run>');
        fprintf(fidWrite, '\n');
        fprintf(fidWrite, '</mzML>');
        fprintf(fidWrite, '\n');
        
        if obj.MZMLDump{2}(1:5) == '<inde'
            fprintf(fidWrite, '</indexedmzML>');
        end
        
        
        fclose(fidWrite);
        
    end


    function options = checkVarargin(varargin)
        % CHECKVARARGIN is used to check the optional paramters
        % and create the options parameter.
        
        % 1- Defaults optional parameters
        options.FolderOut = '';
        options.FileName = '';
        
        % 2- Decipher varargin and update options when relevamt
        input = @(x) find(strcmpi(varargin,x),1);
        
        tgtIx = input('folder');
        if ~isempty(tgtIx)
            options.FolderOut = varargin{tgtIx +1};
        end
        
        tgtIx = input('name');
        if ~isempty(tgtIx)
            options.FileName = [varargin{tgtIx +1}, '.mzML'];
        end
        
        tgtIx = input('ext');
        if ~isempty(tgtIx)
            ext              = varargin{tgtIx +1};
            options.FileName = [obj.FileID, '_', ext,  '.mzML'];
        end
        
        tgtIx = input('mzLim');
        if ~isempty(tgtIx)
            mzLim        = varargin{tgtIx +1};
            options.YLim = [min(mzLim) max(mzLim)];
        end
        
    end

    function [m, n] = findFirstInCellArray(cellArray, pattern)
        m = [];
        n = [];
        
        for ii = 1 : max(size(cellArray));
            if any(strfind(cellArray{ii}, pattern))
                m = ii;
                n = strfind(cellArray{ii}, pattern);
                break
            end
        end
    end



end

