function dtsOut = remSpikes(dtsIn, SpkSz, path2Dat, m)

dtsOut = dtsIn;
dtsOut.Title = 'Profile dataset';
dtsOut.LAST = Trace;
dtsOut.CreatedFrom = NaN;
dtsOut.Log = '';
if nargin == 3
    m = NaN;
end

allProfiles = dtsIn.TimeAxe.Data;
allProfiles(:,3) = 0;
axeMZ = dtsIn.MzAxe.Data;
axeMZ(:,4) = 0;

h = waitbar(0,'processing scans');
for ii = 1:length(allProfiles(:,1))
    
    structInfo.title        = ['Profile scan #', num2str(ii)];
    structInfo.traceType    = dtsOut.Format;
    structInfo.XLabel       = dtsOut.YLabel;
    structInfo.XUnit        = dtsOut.YUnit;
    structInfo.YLabel       = dtsOut.ZLabel;
    structInfo.YUnit        = dtsOut.ZUnit;
    structInfo.link2file    = path2Dat ;
    structInfo.Variables    = 0;
    structInfo.Precision    = 'single';
    structInfo.Path2Finnee  = '';
    structInfo.Log          = ['PRFMSSCN DTS=', num2str(m),...
        ' SPK=', num2str(SpkSz)];
    
    waitbar(ii/length(allProfiles(:,1)));
    XMS = xpend(dtsIn, dtsIn.ListOfScans{ii});
    
    % find spikes
    if SpkSz == 2
        findZeros = find(XMS(:,2) == 0);
        ind2null = findZeros(diff(findZeros) > 2 ...
            & diff(findZeros) <=3);
        XMS(ind2null+1, 2) = 0;
    end
    if SpkSz == 1 || SpkSz == 2  
        findZeros = find(XMS(:,2) == 0);
        ind2null = findZeros(diff(findZeros) > 1 ...
            & diff(findZeros) <=2);
        XMS(ind2null+1, 2) = 0;
    end
    
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
end

if exist('h','var'), close(h); end
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

