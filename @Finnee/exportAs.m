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

if isempty(options.FileOut)
    txtStg              = 'Export Dataset';
    ext                 = {'*.mzML'};
    [pathstr,name,extO] = fileparts(obj.FileIn) ;
    propal              = fullfile(pathstr, [name, '_reduced']);
    
    [fileName, pathName, filterIndex] = uiputfile(ext, txtStg, propal);
    if ~ischar(fileName) && ~ischar(pathName)
        error('myApp:argChk', 'User cancel file selection');
    end
    options.FileOut = fullfile(pathName, fileName);
end

switch filterIndex
    case 1
        ExportAsmzML
end

%% CHECKVARARGIN

    function ExportAsmzML
        fidRead  = fopen(obj.FileIn, 'r');
        fidWrite = fopen(options.FileOut, 'w');
        
        line        = fgetl(fidRead);
        boolRun     = false;
        oneSpectrum = {};
        
        while ischar(line)
            if any(strfind(line, '<run '))
                boolRun = true;
                fprintf(fidWrite, line);
                fprintf(fidWrite, '\n');
                headingLine = fgetl(fidRead);
                line = fgetl(fidRead);
                
            elseif any(strfind(line, '</run>'))
                boolRun = false;
            end
            
            if boolRun
                if any(strfind(line, '</spectrum>'))
                    oneSpectrum{end+1} = line;
                    break
                else
                    oneSpectrum{end+1} = line;
                end
            else
                line = strrep(line, '\', '\\');
                line = strrep(line, '%', '%%');
                line = strrep(line, '''', '''''');
                fprintf(fidWrite, line);
                fprintf(fidWrite, '\n');
            end
            line = fgetl(fidRead);
        end
        
        fclose(fidRead);
        printNewScans
        fprintf(fidWrite, '      </spectrumList>');
        fprintf(fidWrite, '\n');
        fprintf(fidWrite, '    </run>');
        fprintf(fidWrite, '\n');
        fprintf(fidWrite, '  </mzML>');
        fprintf(fidWrite, '\n');
        fprintf(fidWrite, '</indexedmzML>');
        fclose(fidWrite);
        
        function printNewScans
            AxisX     = obj.Datasets{dts}.AxisX.Data;
            
            nbrScans = size(AxisX,1);
            
            % 1. Edit heading
            
            IdS = strfind(headingLine, 'count="') + 6;
            IdE =  strfind(headingLine(IdS:end), '"');
            
            NewHeading = [headingLine(1:IdS(1)), num2str(nbrScans), ...
                headingLine(IdS(1)+IdE(2)-1:end)];
            fprintf(fidWrite, NewHeading);
            fprintf(fidWrite, '\n');
            
            for ii = 1:length(AxisX)
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
                    ItArray =  base64encode(ItArray);
                
                
                [m, ~] = findFirstInCellArray(oneSpectrum, '<spectrum');
                currentLine = oneSpectrum{m};
                msp = m;
                IdS = strfind(currentLine, 'index="') + 6;
                IdE =  strfind(currentLine(IdS:end), '"');
                currentLine = [currentLine(1:IdS(1)), num2str(ii-1), ...
                    currentLine(IdS(1)+IdE(2)-1:end)];
                
                IdS = strfind(currentLine, 'cycle=') + 5;
                IdE =  strfind(currentLine(IdS:end), ' ');
                currentLine = [currentLine(1:IdS(1)), num2str(ii), ...
                    currentLine(IdS(1)+IdE(1)-1:end)];
                
                IdS = strfind(currentLine, 'defaultArrayLength="') + 19;
                IdE =  strfind(currentLine(IdS:end), '"');
                currentLine = [currentLine(1:IdS(1)), num2str(size(MS,1)), ...
                    currentLine(IdS(1)+IdE(2)-1:end)];
                oneSpectrum{m} = currentLine;
                
                [m, ~] = findFirstInCellArray(oneSpectrum, '"scan start time"');
                currentLine = oneSpectrum{m};
                IdS = strfind(currentLine, 'value="') + 6;
                IdE =  strfind(currentLine(IdS:end), '"');
                currentLine = [currentLine(1:IdS(1)), num2str(AxisX(ii)), ...
                    currentLine(IdS(1)+IdE(2)-1:end)];
                oneSpectrum{m} = currentLine;
                
                [m, ~] = findFirstInCellArray(oneSpectrum, 'name="m/z array"');
                currentLine = oneSpectrum{m+1};
                IdS = strfind(currentLine, '<binary>') + 7;
                IdE =  strfind(currentLine, '</binary>');
                currentLine = [currentLine(1:IdS(1)), mzArray, ...
                    currentLine(IdE(1):end)];
                oneSpectrum{m+1} = currentLine;
                
                currentLine = oneSpectrum{m-3};
                IdS = strfind(currentLine, 'encodedLength="') + 14;
                IdE =  strfind(currentLine(IdS:end), '"');
                currentLine = [currentLine(1:IdS(1)), num2str(size(mzArray,2)), ...
                    currentLine(IdS(1)+IdE(2)-1:end)];
                oneSpectrum{m-3} = currentLine;
                
                [m, ~] = findFirstInCellArray(oneSpectrum, 'name="intensity array"');
                currentLine = oneSpectrum{m+1};
                IdS = strfind(currentLine, '<binary>') + 7;
                IdE =  strfind(currentLine, '</binary>');
                currentLine = [currentLine(1:IdS(1)), ItArray, ...
                    currentLine(IdE(1):end)];
                oneSpectrum{m+1} = currentLine;
                
                currentLine = oneSpectrum{m-3};
                IdS = strfind(currentLine, 'encodedLength="') + 14;
                IdE =  strfind(currentLine(IdS:end), '"');
                currentLine = [currentLine(1:IdS(1)), num2str(size(ItArray,2)), ...
                    currentLine(IdS(1)+IdE(2)-1:end)];
                oneSpectrum{m-3} = currentLine;
                
                for jj = 1:size(oneSpectrum, 2)
                    fprintf(fidWrite,'%s', oneSpectrum{jj});
                    fprintf(fidWrite, '\n');
                end
                
            end
            
            
        end
    end


    function options = checkVarargin(varargin)
        % CHECKVARARGIN is used to check the optional paramters
        % and create the options parameter.
        
        % 1- Defaults optional parameters
        options.FileOut     = '';
        options.Overwrite   = false;
        options.FileFormat  = 'mzML';
        options.XLim        = [0 inf];
        options.YLim        = [0 inf];
        
        % 2- Decipher varargin and update options when relevamt
        input = @(x) find(strcmpi(varargin,x),1);
        
        tgtIx = input('fileOut');
        if ~isempty(tgtIx);
            options.FileOut = varargin{tgtIx +1};
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

