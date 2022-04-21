function BestIOns = myAnnotation(MSSpectra, DatabaseIons)

Rp = 60000;% Estimated Resolving power of the MS instrument
Thrs_mz = 1;% ppm
MzInt = 3.3;
Thres_coeff = 0.002;
MSSpectra = sortrows(MSSpectra, 1);
BestIOns = table;

% 1.Find all isotopics envelopes.
MSSpectra(2:end, 3) = diff(MSSpectra(:, 1)) > MzInt;
MSSpectra(1, 3) = 1;

iS = 1;
while 1
    iE = find(MSSpectra(iS+1:end, 3) == 1, 1, 'first');
    if iE > 1
        mySpectra = MSSpectra(iS:iS+iE-1, 1:2);
        
        tgtMz = mySpectra(mySpectra(:,2)==max(mySpectra(:,2)), 1);
        tgtMz(2) = tgtMz(1)-tgtMz(1)*Thrs_mz/1000000;
        tgtMz(3) = tgtMz(1)+tgtMz(1)*Thrs_mz/1000000;
        % Possibility in DatabaseIons
        ix = find(DatabaseIons.baseIonMass >= tgtMz(2) ...
            & DatabaseIons.baseIonMass <= tgtMz(3));
        
        for ii = 1:length(ix)
            testSpectra = DatabaseIons.IsotopicPattern{ix(ii)};
            
            mzLim = [min([mySpectra(:,1); testSpectra(:,1)])...
                max([mySpectra(:,1); testSpectra(:,1)])];
            step  = (mzLim(1)/(Rp*2.354*5));
            
            % Creation of the Master mz axis
            aMZ   = mzLim(1)-4*mzLim(1)/(Rp*2.354);
            while 1
                aMZ(end+1) = aMZ(end)+step;
                
                if aMZ(end) > mzLim(2)+4*mzLim(2)/(Rp*2.354)
                    break
                end
                step  = (aMZ(end)/(Rp*2.354*5));
            end
            
            XY4Ref  =  centr2prof(mySpectra, Rp, aMZ);
            XY2algn = centr2prof(testSpectra, Rp, aMZ);
            [Shift, Coeff] = doAlignMSProfile(XY4Ref, XY2algn, Thrs_mz);
            if Coeff <= Thres_coeff
                TargetSpectra{1} = mySpectra;
                BestIOns = [BestIOns; [DatabaseIons(ix(ii), :), table(Shift, Coeff, TargetSpectra)]];
            end
        end
        
        iS = iS + iE; 
    else
       iS = iS + 1; 
    end
    if iS >= size(MSSpectra, 1)
        break
    end
end

    function XYout = centr2prof(XYin, Rp, axis)
        XYout = zeros(size(axis));
        for lp1 = 1: size(XYin, 1)
            XYout = XYout + XYin(lp1, 2)*...
                pdf('normal', axis, XYin(lp1, 1), (XYin(lp1, 1)/Rp)/2.354);
        end
        XYout = [axis; XYout];
    end
end