function dtsOut = doCtr(dtsIn, par4ctr, path2Dir, MaxFileSize, m)
% par4ctr           : is a strcuture witmyFinnh all oprions need for centroid
%   .baseline       : if a baseline correction is need (see AnalyzeThis)
%                   e.g. 'None'
%   .peakpicking    : peak picking method
%                    e.g. 'LmMm:1' 
%   .threshold      : min intensity for a peak to be considered
%   .record         : A (area), I (intensity)

dtsOut = dtsIn;
dtsOut.Title = 'Centroid dataset';
dtsOut.Format ='MS centroid';
dtsOut.LAST = Trace;
dtsOut.CreatedFrom = NaN;
dtsOut.Log = '';
if nargin == 3
    MaxFileSize = 1000000000;
    m = NaN;
end

if nargin == 4
    m = NaN;
end

axeX   = dtsIn.TimeAxe.Data;
dataPrf(:,1) = axeX;

% DATA 4 each scan
structInfo.title        = '';
structInfo.traceType    = 'MS centroid';
structInfo.XLabel       = dtsIn.YLabel;
structInfo.XUnit        = dtsIn.YUnit;
structInfo.YLabel       = dtsIn.ZLabel;
structInfo.YUnit        = dtsIn.ZUnit;
structInfo.Variables    = 0;
structInfo.Precision    = 'single';
structInfo.Path2Fin     = path2Dir;
structInfo.Log          = ['CTRMSSCN DTS=', num2str(m),...
    ' BSL=', par4ctr.baseline,...
    ' MTD=', par4ctr.peakpicking,...
    ' THR=', num2str( par4ctr.threshold),...
    ' RCD=', par4ctr.record];

[~, rndStr] = fileparts(tempname);
fln = 1;
dtsOut.Path2Dat{fln} = fullfile(path2Dir, rndStr);
structInfo.Path2Dat  = dtsOut.Path2Dat{fln};

h = waitbar(0,'Calculating centroids');
for ii = 1:length(axeX)
    waitbar(ii/length(axeX))
    
    XMS = dtsIn.xpend(dtsIn.ListOfScans{ii}, false);
    assignin('base', 'XMS', XMS)
    AT = AnalyzeThis(XMS, 'baseline', par4ctr.baseline, ...
        'peakpicking', par4ctr.peakpicking, ...
        'threshold', par4ctr.threshold);
    assignin('base', 'AT', AT)
    
    switch par4ctr.record
        case 'I'
            MS = [ AT.PeakList.data(:,4), AT.PeakList.data(:,7)];
        case'A'
            MS = [ AT.PeakList.data(:,4), AT.PeakList.data(:,3)];
    end
    
    dataPrf(ii,2) = sum(MS(:,2));
    dataPrf(ii,3) = max(MS(:,2));
    structInfo.title        = ['Centroid scan #', num2str(ii)];
    structInfo.Variables    = [];
    dtsOut.ListOfScans{ii}  = Trace(structInfo, MS);
    
    s = dir(structInfo.Path2Dat);
    if s.bytes > MaxFileSize;
        [~, rndStr] = fileparts(tempname);
        fln = fln + 1;
        dtsOut.Path2Dat{fln}    = fullfile(path2Dir, rndStr) ;
        structInfo.Path2Dat     = dtsOut.Path2Dat{fln} ;
    end
end

try
    close(h)
catch
end

structInfo.XLabel 	= dtsOut.XLabel;
structInfo.XUnit  	= dtsOut.XUnit;
structInfo.Log     	= ['AXETIME DTS=', num2str(m)];
dtsOut.TimeAxe      = Axe(structInfo, dataPrf(:,1));

% structInfo.XLabel 	= dtsOut.YLabel;
% structInfo.XUnit 	= dtsOut.YUnit;
% structInfo.Log  	= ['AXEMZ DTS=', num2str(m)];
dtsOut.MzAxe    	= Axe();

structInfo.title   		= 'Base Peak Profiles';
structInfo.traceType    = 'ion profile';
structInfo.Log          = ['BPP DTS=', num2str(m)];
structInfo.XLabel       = dtsOut.XLabel;
structInfo.XUnit      	= dtsOut.XUnit;
structInfo.YLabel     	= dtsOut.ZLabel;
structInfo.YUnit        = dtsOut.ZUnit;
dtsOut.BPP              = Trace(structInfo, [dataPrf(:,1), dataPrf(:,3)]);

structInfo.title    	= 'Total Ion profiles';
structInfo.Log      	= ['TIP DTS=', num2str(m)];
dtsOut.TIP           	= Trace(structInfo, [dataPrf(:,1), dataPrf(:,2)]);

% structInfo.title    	= 'Total Ion Scan';
% structInfo.Log          = ['TIS DTS=', num2str(m)];
% structInfo.traceType 	= 'MS profile';
% structInfo.XLabel    	= dtsOut.YLabel;
% structInfo.XUnit      	= dtsOut.YUnit;
% structInfo.YLabel     	= dtsOut.ZLabel;
% structInfo.YUnit    	= dtsOut.ZUnit;
dtsOut.TIS              = Trace();

% structInfo.title        = 'Frequency Ion Scan';
% structInfo.Log      	= ['FIS DTS=', num2str(m)];
% structInfo.traceType  	= 'MS profile';
dtsOut.FIS              = Trace();

% structInfo.title        = 'Base Ion Scan';
% structInfo.Log          = ['BIS DTS=', num2str(m)];
% structInfo.traceType   	= 'MS profile';
dtsOut.BIS              = Trace();

end

