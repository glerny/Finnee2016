%% Description
% Description and more information @ 
% https://github.com/glerny/Finnee2016/wiki/Baseline-and_noise-correction
%
%% Copyright 
% BSD 3-Clause License
% Copyright 2016-2017 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function dtsOut = doCtr(dtsIn, algo, path2Dir, MaxFileSize, m)

% 1. Initiation

dtsOut             = dtsIn;
dtsOut.Title       = 'Datasetd converted to centroid mode';
dtsOut.Format      = 'MS centroid';
dtsOut.LAST        = Trace;
dtsOut.CreatedFrom = NaN;
dtsOut.Log         = '';
if nargin == 3
    MaxFileSize = 1000000000;
    m           = NaN;
end

if nargin == 4
    m = NaN;
end

axeX        = dtsIn.TimeAxe.Data;
axeY        = dtsIn.MzAxe.Data;
[~, rndStr] = fileparts(tempname);
fln         = 1;
dtsOut.Path2Dat{fln} = fullfile(path2Dir, rndStr);
structInfo.Path2Dat  = dtsOut.Path2Dat{fln};

% DATA 4 each scan
structInfo.title        = '';
structInfo.traceType    = dtsOut.Format;
structInfo.XLabel       = dtsOut.YLabel;
structInfo.XUnit        = dtsOut.YUnit;
structInfo.YLabel       = dtsOut.ZLabel;
structInfo.YUnit        = dtsOut.ZUnit;
structInfo.Variables    = 0;
structInfo.Precision    = 'single';
structInfo.Path2Fin     = path2Dir;
strAlg                  = algo.name;
for ii = 1:length(algo.prm)
    strAlg = [strAlg,':', num2str(algo.prm{ii})];
end

structInfo.Log          = ['CTRMSSCN DTS=', num2str(m),...
    ' ALG=', upper(strAlg)];

dataPrf(:,1) = axeX;
dataPrf(:,3) = 0;

% 2. performing centroid for each scans

h = waitbar(0,'Correcting spectra');
for ii = 1:length(axeX)
    waitbar(ii/length(axeX))
    structInfo.title        = ['Centroid scan #', num2str(ii)];
    structInfo.Variables    = 0;
    
    MS = dtsIn.ListOfScans{ii};
    if isempty(MS.Data)
        dtsOut.ListOfScans{ii}  = Trace(structInfo, []);
        
    else
        switch algo.name
            case 'localmax'
                MSr = accuMassByLocMax(MS.Data, algo.prm{1});
        end
        dataPrf(ii,2)           = sum(MSr(:,2));
        dataPrf(ii,3)           = max(MSr(:,2));
        structInfo.title        = ['Centroid scan #', num2str(ii)];
        structInfo.Variables    = 0;
        dtsOut.ListOfScans{ii}  = Trace(structInfo, MSr);
        
        s = dir(structInfo.Path2Dat);
        if s.bytes > MaxFileSize;
            [~, rndStr] = fileparts(tempname);
            fln = fln + 1;
            dtsOut.Path2Dat{fln}    = fullfile(path2Dir, rndStr) ;
            structInfo.Path2Dat     = dtsOut.Path2Dat{fln} ;
        end
    end
end

try
    close(h)
catch
end

structInfo.XLabel   	= dtsOut.XLabel;
structInfo.XUnit    	= dtsOut.XUnit;
structInfo.Log      	= ['AXETIME DTS=', num2str(m)];
dtsOut.TimeAxe          = Axe(structInfo, dataPrf(:,1));

structInfo.XLabel   	= dtsOut.YLabel;
structInfo.XUnit     	= dtsOut.YUnit;
structInfo.Log  	    = ['AXEMZ DTS=', num2str(m)];
dtsOut.MzAxe    	    = Axe;

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

structInfo.title    	= 'Total Ion Scan';
structInfo.Log          = ['TIS DTS=', num2str(m)];
structInfo.traceType 	= 'MS profile';
structInfo.XLabel    	= dtsOut.YLabel;
structInfo.XUnit      	= dtsOut.YUnit;
structInfo.YLabel     	= dtsOut.ZLabel;
structInfo.YUnit    	= dtsOut.ZUnit;
dtsOut.TIS              = Trace;

structInfo.title        = 'Frequency Ion Scan';
structInfo.Log      	= ['FIS DTS=', num2str(m)];
structInfo.traceType  	= 'MS profile';
dtsOut.FIS              = Trace;

structInfo.title        = 'Base Ion Scan';
structInfo.Log          = ['BIS DTS=', num2str(m)];
structInfo.traceType   	= 'MS profile';
dtsOut.BIS              = Trace;

end

