function [MS, flag, cdata] = getIsoEnv_ChemCal(mf, Options)
Thres = 0.001;
Res   = 0.0000001;

url  = 'http://www.chemcalc.org/chemcalc/mf';
flag = 0;
ips = strfind(mf, '+');
ins = strfind(mf, '-');

if ~isempty(ips)
    mf = [mf(1:ips-1), '(',  mf(ips:end), ')'];
end
if ~isempty(ins)
    mf = [mf(1:ins-1), '(',  mf(ins:end), ')'];
end

while 1
    try
        cdata = webread(url, 'mf',  mf, 'threshold', Thres, 'resolution', Res, 'isotopomers', 'xy');
        MS = str2num( regexprep(cdata.xy, '[\n\r]+',';'));
        break
    catch EM
        disp(EM) 
    end
end

