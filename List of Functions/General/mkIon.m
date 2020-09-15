function [mf, ea] = mkIon(ea, SglAdduct)
mf = '';

if size(ea, 1) == 1
    eaN = table();
    eaN.number = ea.number;
    eaN.percentage = ea.percentage;
    eaN.element{1} = ea.element;
    
    ea = eaN;
    clear eaN;
end

test = SglAdduct.C;
ele  = 'C';
if test ~= 0
    if any(strcmp(ele, ea.element))
        ea.number(strcmp(ele, ea.element)) = ...
            ea.number(strcmp(ele, ea.element)) + test;
    else
        warning off
        ea.element{end+1} = ele;
        ea.number(end)    = test;
        warning on
    end
end

test = SglAdduct.H;
ele  = 'H';
if test ~= 0
    if any(strcmp(ele, ea.element))
        ea.number(strcmp(ele, ea.element)) = ...
            ea.number(strcmp(ele, ea.element)) + test;
    else
        warning off
        ea.element{end+1} = ele;
        ea.number(end)    = test;
        warning on
    end
end

test = SglAdduct.O;
ele  = 'O';
if test ~= 0
    if any(strcmp(ele, ea.element))
        ea.number(strcmp(ele, ea.element)) = ...
            ea.number(strcmp(ele, ea.element)) + test;
    else
        warning off
        ea.element{end+1} = ele;
        ea.number(end)    = test;
        warning on
    end
end

test = SglAdduct.N;
ele  = 'N';
if test ~= 0
    if any(strcmp(ele, ea.element))
        ea.number(strcmp(ele, ea.element)) = ...
            ea.number(strcmp(ele, ea.element)) + test;
    else
        warning off
        ea.element{end+1} = ele;
        ea.number(end)    = test;
        warning on
    end
end

test = SglAdduct.S;
ele  = 'S';
if test ~= 0
    if any(strcmp(ele, ea.element))
        ea.number(strcmp(ele, ea.element)) = ...
            ea.number(strcmp(ele, ea.element)) + test;
    else
        warning off
        ea.element{end+1} = ele;
        ea.number(end)    = test;
        warning on
    end
end

test = SglAdduct.Na;
ele  = 'Na';
if test ~= 0
    if any(strcmp(ele, ea.element))
        ea.number(strcmp(ele, ea.element)) = ...
            ea.number(strcmp(ele, ea.element)) + test;
    else
        warning off
        ea.element{end+1} = ele;
        ea.number(end)    = test;
        warning on
    end
end

test = SglAdduct.K;
ele  = 'K';
if test ~= 0
    if any(strcmp(ele, ea.element))
        ea.number(strcmp(ele, ea.element)) = ...
            ea.number(strcmp(ele, ea.element)) + test;
    else
        warning off
        ea.element{end+1} = ele;
        ea.number(end)    = test;
        warning on
    end
end

test = SglAdduct.Cu;
ele  = 'Cu';
if test ~= 0
    if any(strcmp(ele, ea.element))
        ea.number(strcmp(ele, ea.element)) = ...
            ea.number(strcmp(ele, ea.element)) + test;
    else
        warning off
        ea.element{end+1} = ele;
        ea.number(end)    = test;
        warning on
    end
end

test = SglAdduct.F;
ele  = 'F';
if test ~= 0
    if any(strcmp(ele, ea.element))
        ea.number(strcmp(ele, ea.element)) = ...
            ea.number(strcmp(ele, ea.element)) + test;
    else
        warning off
        ea.element{end+1} = ele;
        ea.number(end)    = test;
        warning on
    end
end

test = SglAdduct.Br;
ele  = 'Br';
if test ~= 0
    if any(strcmp(ele, ea.element))
        ea.number(strcmp(ele, ea.element)) = ...
            ea.number(strcmp(ele, ea.element)) + test;
    else
        warning off
        ea.element{end+1} = ele;
        ea.number(end)    = test;
        warning on
    end
end

test = SglAdduct.Cl;
ele  = 'Cl';
if test ~= 0
    if any(strcmp(ele, ea.element))
        ea.number(strcmp(ele, ea.element)) = ...
            ea.number(strcmp(ele, ea.element)) + test;
    else
        warning off
        ea.element{end+1} = ele;
        ea.number(end)    = test;
        warning on
    end
end

if any(ea.number < 0)
    mf = '';
    ea = table();
else
    ea.percentage = nan(size(ea.percentage));
    for ii = 1:size(ea, 1)
        if ea.number(ii) > 0
            mf = [mf [ea.element{ii} num2str(ea.number(ii))]];
        end
    end
    
    switch SglAdduct.Charge
        case 1
            mf = [mf, '(+1)'];
        case 2
            mf = [mf, '(+2)'];
        case 3
            mf = [mf, '(+3)'];
        case 4
            mf = [mf, '(+4)'];
        case 5
            mf = [mf, '(+5)'];
        case -1
            mf = [mf, '(-1)'];
        case -2
            mf = [mf, '(-2)'];
        case -3
            mf = [mf, '(-3)'];
        case -4
            mf = [mf, '(-4)'];
        case -5
            mf = [mf, '(-5)'];
    end
end

