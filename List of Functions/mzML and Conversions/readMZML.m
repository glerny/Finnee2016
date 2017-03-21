%% DESCRIPTION
%
%% Copyright 
% BSD 3-Clause License
% Copyright 2016-2017 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


function mzMLDump = readMZML( fidRead, stopLabel)

frewind(fidRead);
mzMLDump.fLine = fgetl(fidRead);
mzML{1}        = mzMLDump.fLine;
curLine        = fgetl(fidRead);
mzML{2}        = curLine;
[LDR, ~]       = getMZMLCamp(curLine, fidRead);
[outStrc, ~, mzML] = addChildNode(LDR, fidRead, true, stopLabel, mzML);
mzMLDump = mzML;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SUB FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [outStrc, SB, mzML] = addChildNode(LDR, fidRead, SB, stopLabel, mzML)
outStrc = struct;
outStrc.(LDR.label).attributes = LDR.attributes;
outStrc.(LDR.label).text= LDR.text;

if strcmp(stopLabel, LDR.label)
    SB = false;
    return
end

if ~SB
    return
end
        
if LDR.open
    endLabel = ['/', LDR.label];
    parentLabel = LDR.label;
    while SB
        curLine     = fgetl(fidRead);
        mzML{end+1} = curLine;
        [LDR, ~]    = getMZMLCamp(curLine, fidRead);
        
        if strcmp(stopLabel, LDR.label)
            SB = false;
            break
        end
        
        if strcmp(endLabel, LDR.label)
            break
        end
        if ~isfield(outStrc.(parentLabel), LDR.label)
            outStrc.(parentLabel).(LDR.label) = {};
        end
        [s, SB, mzML] = addChildNode(LDR, fidRead, SB, stopLabel, mzML);
        
        outStrc.(parentLabel).(LDR.label){end+1} = s.(LDR.label);
    end
end
end


