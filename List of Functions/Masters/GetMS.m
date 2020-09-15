function MS = GetMS(mf, cAdduct, n)
url1     = 'http://www.chemcalc.org/chemcalc/mf';
url2     = 'http://www.chemcalc.org/chemcalc/em';
MS = [];

cmf = [];
for ll = 1:n
    cmf = [cmf, char(mf)];
end


if char(cAdduct.mf) == ''''
else
    cmf = [cmf, char(cAdduct.mf)];
end

switch cAdduct.Charge
    case 1
        cmf = [cmf, '(+1)'];
        
    case 2
        cmf = [cmf, '(+2)'];
        
    case 3
        cmf = [cmf, '(+3)'];
end

try
    data = webread(url1, 'mf', cmf, 'resolution', '0.0000001',...
        'isotopomers', 'jcamp,xy');
    
    MS = str2num( regexprep(data.xy, '[\n\r]+',';'));
    MS(MS(:,2) < 0.01, :) = [];
    
catch ME
    
    disp(ME)
    disp(mf)
end

end