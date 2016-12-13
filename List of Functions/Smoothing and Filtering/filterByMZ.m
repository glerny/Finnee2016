function dtsOut = filterByMZ(dtsIn, filter, path2Dir, MaxFileSize, m)

dtsOut = dtsIn;
dtsOut.Title = 'Profile dataset';
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

allProfiles = dtsIn.TimeAxe.Data;
allProfiles(:,3) = 0;
axeMZ = dtsIn.MzAxe.Data;
axeMZ(:,4) = 0;
[~, rndStr] = fileparts(tempname);
fln = 1;


h = waitbar(0,'processing scans');
for ii = 1:length(allProfiles(:,1))
    dtsOut.Path2Dat{fln} = fullfile(path2Dir, rndStr) ;
    
    structInfo.title        = ['Profile scan #', num2str(ii)];
    structInfo.traceType    = dtsOut.Format;
    structInfo.XLabel       = dtsOut.YLabel;
    structInfo.XUnit        = dtsOut.YUnit;
    structInfo.YLabel       = dtsOut.ZLabel;
    structInfo.YUnit        = dtsOut.ZUnit;
    structInfo.Path2Dat     = dtsOut.Path2Dat{fln} ;
    structInfo.Variables    = 0;
    structInfo.Precision    = 'single';
    structInfo.Path2Fin     = path2Dir;
    structInfo.Log          = ['PRFMSSCN DTS=', num2str(m),...
        ' FMZ=1'];
    
    waitbar(ii/length(allProfiles(:,1)));
    XMS = xpend(dtsIn, dtsIn.ListOfScans{ii});
    
    % FilterOut
    XMS(filter, 2) = 0;
    
    
    allProfiles(ii, 2) = sum(XMS(:,2));
    allProfiles(ii, 3) = max(XMS(:,2));
    axeMZ(:,2) = axeMZ(:,2) + XMS(:,2);
    axeMZ(:,4) = max([axeMZ(:,4), XMS(:,2)], [], 2);
    indNotZero = XMS(:,2) > 0;
    axeMZ(:,3) = axeMZ(:,3) + indNotZero;
                
    % remove trailing zeros
     provMat = [XMS(2:end, 2); 0];
     provMat(:,2) = XMS(:, 2);
     provMat(:,3) = [0; XMS(1:end-1, 2)];
     MS = XMS(sum(provMat, 2) > 0, :);
    
    dtsOut.ListOfScans{ii} = Trace(structInfo, MS);
    s = dir(structInfo.Path2Dat);
    if s.bytes > MaxFileSize;
        [~, rndStr] = fileparts(tempname);
        fln = fln + 1;
        dtsOut.Path2Dat{fln}    = fullfile(path2Dir, rndStr) ;
        structInfo.Path2Dat     = dtsOut.Path2Dat{fln} ;
    end
    
    
end

if exist('h','var'), close(h); end

%Remove trailing zeros in axeMZ
provMat = [axeMZ(2:end, 2); 0];
provMat(:,2) = axeMZ(:, 2);
provMat(:,3) = [0; axeMZ(1:end-1, 2)];
axeMZ = axeMZ(sum(provMat, 2) > 0, :);

structInfo.XLabel     	= dtsOut.XLabel;
structInfo.XUnit      	= dtsOut.XUnit;
structInfo.Log      	= ['AXETIME DTS=', num2str(m)];
dtsOut.TimeAxe          = Axe(structInfo, allProfiles(:,1));

structInfo.XLabel       = dtsOut.YLabel;
structInfo.XUnit        = dtsOut.YUnit;
structInfo.Log          = ['AXEMZ DTS=', num2str(m)];
dtsOut.MzAxe            = Axe(structInfo, axeMZ(:,1));

structInfo.title    	= 'Base Peak Profiles';
structInfo.traceType 	= 'ion profile';
structInfo.Log       	= ['BPP DTS=', num2str(m)];
structInfo.XLabel       = dtsOut.XLabel;
structInfo.XUnit        = dtsOut.XUnit;
structInfo.YLabel   	= dtsOut.ZLabel;
structInfo.YUnit        = dtsOut.ZUnit;
dtsOut.BPP              = Trace(structInfo, [allProfiles(:,1), allProfiles(:,3)]);

structInfo.title        = 'Total Ion profiles';
structInfo.Log       	= ['TIP DTS=', num2str(m)];
dtsOut.TIP              = Trace(structInfo, [allProfiles(:,1), allProfiles(:,2)]);

structInfo.title        = 'Total Ion Scan';
structInfo.Log          = ['TIS DTS=', num2str(m)];
structInfo.traceType 	= 'MS profile';
structInfo.XLabel   	= dtsOut.YLabel;
structInfo.XUnit    	= dtsOut.YUnit;
structInfo.YLabel     	= dtsOut.ZLabel;
structInfo.YUnit        = dtsOut.ZUnit;
dtsOut.TIS              = Trace(structInfo, [axeMZ(:,1), axeMZ(:,2)]);

structInfo.title     	= 'Frequency Ion Scan';
structInfo.Log          = ['FIS DTS=', num2str(m)];
structInfo.traceType 	= 'MS profile';
dtsOut.FIS              = Trace(structInfo, [axeMZ(:,1), axeMZ(:,3)]);

structInfo.title     	= 'Base Ion Scan';
structInfo.Log       	= ['BIS DTS=', num2str(m)];
structInfo.traceType   	= 'MS profile';
dtsOut.BIS              = Trace(structInfo, [axeMZ(:,1), axeMZ(:,4)]);
end

