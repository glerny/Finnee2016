function Match = MatchSpectra(MS, RefSpectra, Motif, ppm)
Rp = 60000;

Match.Ix        = [];
Match.Adducts   = {};
Match.ReptMotif = [];
Match.Charge    = [];
Match.mf        = {};
Match.Ions      = {};
Match.mz        = [];
Match.MSTheo    = {};
Match.MSComp    = {};
Match.Shift     = [];
Match.Coeff     = [];
Match.Ratio     = [];

Spectra2Test    = find(strcmp(Motif, RefSpectra.Motif));

for ii = 1:size(Spectra2Test, 1)
    SPE = [];
    mz  = RefSpectra.BPmz(Spectra2Test(ii));
    
    if any(abs(MS(:,1) - mz)/mz*1000000 <= ppm)
        %
        IP = find(abs(MS(:,1) - mz)/mz*1000000 <= ppm);
        if length(IP) >= 2
            [~, iIP] = min(abs(MS(IP,1) - mz));
            IP = IP(iIP);
        end
        SPE = [SPE; MS(IP, :)];
        
        
        MSth =  RefSpectra.MS{Spectra2Test(ii)};
        
        %mk axis
        sigAve = (sum(MSth(:,1).*MSth(:,2))/sum(MSth(:,2))/Rp)/2.354;
        Axis   = (min(MSth(:,1))-4*sigAve:sigAve/20:...
            max(MSth(:,1))+4*sigAve)';
        
        XYtheo = centr2prof(MSth, Rp, Axis);
        XYexp  = centr2prof(MS, Rp, Axis);
        [Shift, Coef] = doAlignMSProfile(XYtheo, XYexp);
        XYexp(:, 3) = interp1(XYtheo(:,1) - Shift, XYtheo(:,2), XYexp(:,1));
        XYexp(isnan(XYexp)) = 0;
        Ratio  = XYexp(:,3)\XYexp(:,2);
        
        Match.Ix(end+1, 1)        = IP;
        Match.Adducts{end+1, 1}   = RefSpectra.Adducts{Spectra2Test(ii)};
        Match.ReptMotif(end+1, 1) = RefSpectra.Repetition(Spectra2Test(ii));
        Match.Charge(end+1, 1)    = RefSpectra.Charge(Spectra2Test(ii));
        Match.mf{end+1, 1}        = RefSpectra.mf{Spectra2Test(ii)};
        Match.Ions{end+1, 1}      = MS(IP, :);
        Match.mz(end+1, 1)        = mz;
        Match.MSTheo{end+1, 1}    = MSth;
        Match.MSComp{end+1, 1}    = XYexp;
        Match.Shift(end+1, 1)     = Shift;
        Match.Coeff(end+1, 1)     = Coef;
        Match.Ratio(end+1, 1)     = Ratio;
    end
    
end
Match = struct2table(Match);

    function XYout = centr2prof(XYin, Rp, axis)
        XYout = zeros(size(axis));
        for lp1 = 1: size(XYin, 1)
            XYout = XYout + XYin(lp1, 2)*...
                pdf('normal', axis, XYin(lp1, 1), (XYin(lp1, 1)/Rp)/2.354);
        end
        XYout = [axis, XYout];
    end
end

