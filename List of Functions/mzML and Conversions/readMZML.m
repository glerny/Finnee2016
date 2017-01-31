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
curLine = fgetl(fidRead);
[LDR, ~] = getMZMLCamp(curLine, fidRead);
[outStrc, ~] = addChildNode(LDR, fidRead, true, stopLabel);
mzMLDump.struct = outStrc;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SUB FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [outStrc, SB] = addChildNode(LDR, fidRead, SB, stopLabel)
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
        curLine = fgetl(fidRead);
        [LDR, ~] = getMZMLCamp(curLine, fidRead);
        
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
        [s, SB] = addChildNode(LDR, fidRead, SB, stopLabel);
        
        outStrc.(parentLabel).(LDR.label){end+1} = s.(LDR.label);
    end
end
end


