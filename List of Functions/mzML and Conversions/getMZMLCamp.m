function [LDR, curLine] = getMZMLCamp(curLine, fidRead)
% Get the data_label and associated data-set from mzML line(s)
curLine = strtrim(curLine); %remove leading and trailing white spaces
if ~strcmp(curLine(1), '<')
    error('myApp:argChk', 'Incorrect data')
end
% find all opening bracket
LDR.label = '';
LDR.attributes = {};
LDR.text = {};
LDR.open = 1;
% find all opening bracket
indSt = strfind(curLine, '<');
if length(indSt) == 1
    % No intercalling text
    if strcmp(curLine(end), '>') %classical case (finish)
        if curLine(end-1) == '/'
            LDR.open = 0;
        else
            LDR.open = 1;
        end
        
        str2decode = curLine(2:end-1);
        ind2att = strfind(str2decode, '="');
        if isempty(ind2att)
            LDR.label = str2decode;
        else
            ind2blk = strfind(str2decode, ' ');
            ind2dbc =  strfind(str2decode, '"');
            LDR.label = str2decode(1:ind2blk(1)-1);
            for ii = 1:length(ind2att)
                istt = find(ind2blk < ind2att(ii), 1, 'last');
                iend = find(ind2dbc > ind2att(ii), 1, 'first');
                istt = ind2blk(istt);
                iend = ind2dbc(iend+1);
                cutstr = str2decode(istt+1:iend-1);
                indEqu = strfind(cutstr, '=');
                LDR.attributes{end+1}.label = cutstr(1:indEqu(1)-1);
                LDR.attributes{end}.field{1} = cutstr(indEqu(1)+2:end);
            end
        end
    else        %case where field as entry in it
        str2decode = curLine(2:end);
        ind2att = strfind(str2decode, '="');
        if isempty(ind2att)
            LDR.label = str2decode;
        else
            ind2blk = strfind(str2decode, ' ');
            LDR.label = str2decode(1:ind2blk(1)-1);
            for ii = 1:length(ind2att)-1
                istt = find(ind2blk < ind2att(ii), 1, 'last');
                iend = find(ind2blk > ind2att(ii), 1, 'first');
                istt = ind2blk(istt);
                if isempty(iend)
                    iend = length(cutstr);
                else
                    iend = ind2blk(iend);
                end
                cutstr = str2decode(istt+1:iend-1);
                indEqu = strfind(cutstr, '=');
                LDR.attributes{end+1}.label = cutstr(1:indEqu(1)-1);
                indGlm = strfind(cutstr, '"');
                LDR.attributes{end}.field{1} = ...
                    cutstr(indGlm(1)+1:indGlm(2)-1);
            end
            ii = ii +1;
            istt = find(ind2blk < ind2att(ii), 1, 'last');
            istt = ind2blk(istt);
            iend = length(str2decode);
            cutstr = str2decode(istt+1:iend-1);
            indEqu = strfind(cutstr, '=');
            LDR.attributes{end+1}.label = cutstr(1:indEqu(1)-1);
            indGlm = strfind(cutstr, '"');
            LDR.attributes{end}.field{1} = ...
                cutstr(indGlm(1)+1:end);
            while 1
                curLine = fgetl(fidRead);
                if strcmp(curLine(end-1:end), '/>')
                    LDR.attributes{end}.field{end+1} = curLine(1:end-2);
                    LDR.open = 0;
                    break
                else
                    LDR.attributes{end}.field{end+1} = curLine;
                    LDR.open = 1;
                end
            end
        end
        
    end
elseif length(indSt) == 2 % More than two starting bracket, check for text
    indEn = strfind(curLine, '>');
    if strcmp(curLine(indSt(1)+1:indEn(1)-1), 'binary')
        LDR.label = 'binary';
        LDR.attributes = {};
        LDR.text = curLine(indEn(1)+1:indSt(2)-1);
    else
        SHOULDNOTHAPPENWITHmzMLFILES
    end
    
else
    SHOULDNOTHAPPENWITHmzMLFILES
end

if ~isempty(LDR.text)
    LDR.open = 0;
end
end
