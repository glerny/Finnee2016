function dtsOut = doBslCorPrf(dtsIn, par4bas, path2Dir, MaxFileSize, m)
tic
dtsOut = dtsIn;
dtsOut.Title = 'Baseline corrected profile dataset';
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

Link2CorPrf = tempname;
nbrProf = sum(par4bas.obj.IndMax);
lengthTm = length(par4bas.obj.TimeAxe.Data);
fidRead = fopen(par4bas.obj.Link2SelPrf, 'rb');
PrfMat = fread(fidRead, [nbrProf, lengthTm], par4bas.obj.Precision);

h = waitbar(0,'Correcting profiles');
for ii = 1:nbrProf
    waitbar(ii/nbrProf)
    Prf = (PrfMat(ii, :))';
    
    switch par4bas.type
        case 'None'
            parameter = par4bas.parameter;
            bsl = zeros(length(yy(:,2)), 1);
            yy = Prf - bsl;
            indNeg = yy <= 0;
            yy(indNeg) = 0;
            
        case 'doPF'
            parameter = par4bas.parameter;
            [baseline, ~ ] = doPF(yy, parameter(1));
            bsl = round(baseline);
            yy = Prf - bsl;
            indNeg = yy <= 0;
            yy(indNeg) = 0;
            
        case 'arPLS'
            parameter = par4bas.parameter;
            notZeros = find(Prf ~= 0);
            [baseline, ~] =  doArPLS(Prf(notZeros),...
                parameter(1), parameter(2));
            bsl = zeros(length(Prf), 1);
            bsl(notZeros) = round(baseline);
            yy = Prf - bsl;
            indNeg = yy <= 0;
            yy(indNeg) = 0;
            
        case 'arPLS2'
            parameter = par4bas.parameter;
            notZeros = find(Prf ~= 0);
            [baseline, ~] =  doArPLS2(Prf(notZeros),...
                parameter(1));
            bsl = zeros(length(Prf), 1);
            bsl(notZeros) = round(baseline);
            yy = Prf - bsl;
            indNeg = yy <= 0;
            yy(indNeg) = 0;
    end
    PrfMat(ii, :) = yy';
end
fclose(fidRead);
try
    close(h)
catch
end

axeX   = dtsIn.TimeAxe.Data;
axeY   = dtsIn.MzAxe.Data;
M4C    = zeros(length(axeY), 2*par4bas.wdz +1);
wdz    = par4bas.wdz;
noise  = par4bas.noise;
[~, rndStr] = fileparts(tempname);
fln = 1;

% DATA 4 each scan
structInfo.title        = '';
structInfo.traceType    = dtsIn.Format;
structInfo.XLabel       = dtsIn.YLabel;
structInfo.XUnit        = dtsIn.YUnit;
structInfo.YLabel       = dtsIn.ZLabel;
structInfo.YUnit        = dtsIn.ZUnit;
structInfo.Variables    = 0;
structInfo.Precision    = 'single';
structInfo.Path2Fin     = path2Dir;
structInfo.Log          = ['PRFMSSCN DTS=', num2str(m),...
    ' CNS=', num2str(par4bas.noise)];


for ii = 1:2*wdz+1
    XMS = dtsIn.xpend(dtsIn.ListOfScans{ii});
    Corr = PrfMat(:,ii);
    XMS(par4bas.obj.IndMax, 2) = Corr;
    M4C(:, ii) = XMS(:,2);
end

dataMZ(:,1)  = axeY;
dataMZ(:,4)  = 0;
dataPrf(:,1) = axeX;
dataPrf(:,3) = 0;

for ii = [1: wdz, (length(axeX) - wdz):length(axeX)]
    dtsOut.Path2Dat{fln} = fullfile(path2Dir, rndStr);
    structInfo.Path2Dat  = dtsOut.Path2Dat{fln};
    
    XMS = dtsIn.xpend(dtsIn.ListOfScans{ii});
    Corr = PrfMat(:,ii);
    XMS(par4bas.obj.IndMax, 2) = Corr;
    dataMZ(:,2)    = dataMZ(:,2) + XMS(:,2);
    iNZ = XMS(:,2) > 0;
    dataMZ(iNZ, 3) = dataMZ(iNZ, 3) + 1;
    dataMZ(:,4)    = max([dataMZ(:,4), XMS(:,2)], [], 2);
    dataPrf(ii,2)  = sum(XMS(:,2));
    dataPrf(ii,3)  = max(XMS(:,2));
    
    provMat                 = [XMS(2:end, 2); 0];
    provMat(:,2)            = XMS(:, 2);
    provMat(:,3)            = [0; XMS(1:end-1, 2)];
    XMS                     = XMS(sum(provMat, 2) > 0, :);
    structInfo.title        = ['Profile scan #', num2str(ii)];
    structInfo.Variables    = dtsIn.ListOfScans{ii}.Variables;
    dtsOut.ListOfScans{ii}  = Trace(structInfo, XMS);
    
    s = dir(structInfo.Path2Dat);
    if s.bytes > MaxFileSize;
        [~, rndStr] = fileparts(tempname);
        fln = fln + 1;
        dtsOut.Path2Dat{fln}    = fullfile(path2Dir, rndStr) ;
        structInfo.Path2Dat     = dtsOut.Path2Dat{fln} ;
    end
    
end


h = waitbar(0,'Correcting spectra');
for ii = 2*wdz+2:length(axeX)
    waitbar(ii/length(axeX))
    
    dtsOut.Path2Dat{fln} = fullfile(path2Dir, rndStr);
    structInfo.Path2Dat  = dtsOut.Path2Dat{fln};
    
    MSCur = M4C(:, wdz+1);
    
    D = zeros(length(MSCur), (2*wdz+1)^2);
    for jj =1: 2*wdz+1
        id1 = max(wdz+2-jj, 1);
        id2 = max(jj-wdz, 1);
        D(id1:end-id2+1,(2*wdz+1)*(jj-1)+1:(2*wdz+1)*jj) = M4C(id2:end-id1+1, :);
    end
    ind2cut = max(D, [], 2) < 3*noise;
    
    MSCur(ind2cut) = 0;
    XMS = axeY;
    XMS(:,2) = MSCur;
    dataMZ(:,2)    = dataMZ(:,2) + XMS(:,2);
    iNZ = XMS(:,2) > 0;
    dataMZ(iNZ, 3) = dataMZ(iNZ, 3) + 1;
    dataMZ(:,4)    = max([dataMZ(:,4), XMS(:,2)], [], 2);
    dataPrf(ii-wdz-1,2)     = sum(XMS(:,2));
    dataPrf(ii-wdz-1,3)     =  max(XMS(:,2));
    
    provMat                 = [XMS(2:end, 2); 0];
    provMat(:,2)            = XMS(:, 2);
    provMat(:,3)            = [0; XMS(1:end-1, 2)];
    XMS                     = XMS(sum(provMat, 2) > 0, :);
    structInfo.title        = ['Profile scan #', num2str(ii-wdz-1)];
    structInfo.Variables    = dtsIn.ListOfScans{ii-wdz-1}.Variables;
    dtsOut.ListOfScans{ii-wdz-1}  = Trace(structInfo, XMS);
    
    s = dir(structInfo.Path2Dat);
    if s.bytes > MaxFileSize;
        [~, rndStr] = fileparts(tempname);
        fln = fln + 1;
        dtsOut.Path2Dat{fln}    = fullfile(path2Dir, rndStr) ;
        structInfo.Path2Dat     = dtsOut.Path2Dat{fln} ;
    end
    
    XMS = dtsIn.xpend(dtsIn.ListOfScans{ii});
    Corr = PrfMat(:,ii);
    XMS(par4bas.obj.IndMax, 2) = Corr;
    M4C = [M4C(:, 2:end), XMS(:,2)];
end
try
    close(h)
catch
end

%Remove trailing zeros in dataMZ
provMat = [dataMZ(2:end, 2); 0];
provMat(:,2) = dataMZ(:, 2);
provMat(:,3) = [0; dataMZ(1:end-1, 2)];
dataMZ = dataMZ(sum(provMat, 2) > 0, :);


structInfo.XLabel 	= dtsOut.XLabel;
structInfo.XUnit  	= dtsOut.XUnit;
structInfo.Log     	= ['AXETIME DTS=', num2str(m)];
dtsOut.TimeAxe      = Axe(structInfo, dataPrf(:,1));

structInfo.XLabel 	= dtsOut.YLabel;
structInfo.XUnit 	= dtsOut.YUnit;
structInfo.Log  	= ['AXEMZ DTS=', num2str(m)];
dtsOut.MzAxe    	= Axe(structInfo, dataMZ(:,1));

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
dtsOut.TIS              = Trace(structInfo, [dataMZ(:,1), dataMZ(:,2)]);

structInfo.title        = 'Frequency Ion Scan';
structInfo.Log      	= ['FIS DTS=', num2str(m)];
structInfo.traceType  	= 'MS profile';
dtsOut.FIS              = Trace(structInfo, [dataMZ(:,1), dataMZ(:,3)]);

structInfo.title        = 'Base Ion Scan';
structInfo.Log          = ['BIS DTS=', num2str(m)];
structInfo.traceType   	= 'MS profile';
dtsOut.BIS              = Trace(structInfo, [dataMZ(:,1), dataMZ(:,4)]);

end

