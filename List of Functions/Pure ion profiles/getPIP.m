function LstPIP = getPIP(LoPts, ThrMZ, ThrIt, minPts, axe, infoDts)
LoPts  = sortrows(LoPts, 1);
FstCut = [1; find(diff(LoPts(:,1)) > ThrMZ)+1; size(LoPts,1)+1];
LstPIP    = {};
for ii = 1:length(FstCut)-1
    LFC    = LoPts(FstCut(ii):FstCut(ii+1)-1, :);
    LFC    = sortrows(LFC, 3);
    SndCut = [1; find(diff(LFC(:,3)) > 1)+1; size(LFC,1)+1];
    for jj = 1:length(SndCut) -1
        cPIP    = LFC(SndCut(jj):SndCut(jj+1)-1, :);
        if size(unique(cPIP(:,3)), 1) >= minPts && max(cPIP(:,2)) >= ThrIt
            cPIP       = sortrows(cPIP, 1);
            ttd        = find(diff(cPIP(:,1)) > ThrMZ)+1;
            if ~isempty(ttd)
                LstPIP = [LstPIP, getPIP(cPIP, ThrMZ, ThrIt, minPts, axe, infoDts)];
            else
                cPIP          = sortrows(cPIP, 3);
                LstPIP{end+1} = PIP(cPIP, axe, infoDts);
            end
        end
    end
end
end

