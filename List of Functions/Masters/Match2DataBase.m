function BestMatch = Match2DataBase(DataBase, FOM, shift, error, threshold1, threshold2)

BestMatch.Metabolite   = {};
BestMatch.ListAdducts  = {};
BestMatch.FirstAdducts = {};
BestMatch.ReptMotif    = [];
BestMatch.Charge       = [];
BestMatch.ListIons     = {};
BestMatch.mf1          = {};
BestMatch.me           = {};
BestMatch.mfall        = {};
BestMatch.Coeff1       = [];
BestMatch.Coeff2       = [];
BestMatch.shift        = {};
BestMatch.Spectra      = {};
BestMatch.allMatch     = {};
BestMatch.Ratio        = [];
BestMatch.nbrAdducts   = [];

Rp      = 60000;
mzMin = min(FOM.mean_AccMass);
mzMaz = max(FOM.mean_AccMass);
step  = (mzMin/(Rp*2.354*5));
aMZ   = mzMin-4*mzMin/(Rp*2.354);

while 1
    aMZ(end+1) = aMZ(end)+step;
    
    if aMZ(end) > mzMaz+4*mzMaz/(Rp*2.354)
        break
    end
    step  = (aMZ(end)/(Rp*2.354*5));
end

FOM = sortrows(FOM,'mean_Area','descend');
MS = [FOM.mean_AccMass, FOM.mean_Area];
MS(:,1) = MS(:, 1) + shift*MS(:, 1)/1000000;

IX = find(abs(DataBase.BPmz - MS(1, 1))/MS(1, 1)*1000000 < error);
M = unique(DataBase.Motif(IX));

if isempty(M)
    return
end

for ii = 1: size(M, 1)
    Match = MatchSpectra(MS,  DataBase, M{ii}, error);
    Match(Match.Coeff > threshold2, :) = [];
    if ~any(Match.Coeff <= threshold1)
        continue
    end
    
    BestMatch.Metabolite{end+1, 1} = M{ii};
    
    Lst = unique(Match.Ix);
    Id4Best = [];
    for jj = 1:length(Lst)
        LstJj =  find(Match.Ix == Lst(jj));
        [~, Imin] = min(Match.Coeff(LstJj));
        Id4Best = [Id4Best; LstJj(Imin)];
    end
    
    %% Find finalcomp
    XYexp  = centr2prof(MS, Rp, aMZ');
    shift = sum(Match.Shift(Id4Best).*MS(Match.Ix(Id4Best), 2))/sum(MS(Match.Ix(Id4Best), 2));
    
    
    for jj = 1:length(Id4Best)
        MS2Ali = Match.MSTheo{Id4Best(jj)};
        MS2Ali(:,1) = MS2Ali(:,1) - shift;
        XY2Ali = centr2prof(MS2Ali, Rp, aMZ');
        Ratio =  XY2Ali((XY2Ali(:,2) ~= 0), 2)\XYexp((XY2Ali(:,2) ~= 0), 2);
        XYexp(:, jj+2) = Ratio*XY2Ali(:,2);
    end
    
    xy = (1:size(XYexp, 1))';
    xy(:,2) = sum(XYexp(:, 2:end), 2);
    xy = trailRem(xy, 2);
    XYexp = XYexp(xy(:,1), :);
    
    BestMatch.ListAdducts{end+1, 1} = Match.Adducts(Id4Best);
    [~, IP] = max(MS(Match.Ix(Id4Best), 2));
    BestMatch.FirstAdducts{end+1, 1} = char(Match.Adducts{Id4Best(IP)});
    BestMatch.ReptMotif(end+1, 1) = Match.ReptMotif(Id4Best(IP));
    BestMatch.Charge(end+1, 1) = Match.Charge(Id4Best(IP));
    BestMatch.ListIons{end+1, 1} = [cell2mat(Match.Ions(Id4Best)), ...
        Match.ReptMotif(Id4Best), Match.Coeff(Id4Best), ...
        Match.mz(Id4Best), Match.Ratio(Id4Best)];
    BestMatch.mf1{end+1, 1}      = Match.mf{Id4Best(IP)};
    BestMatch.me{end+1, 1}       = Match.mz(Id4Best(IP));
    BestMatch.mfall{end+1, 1}    = Match.mf(Id4Best);
    BestMatch.Coeff1(end+1, 1)   = Match.Coeff(Id4Best(IP));
    BestMatch.Coeff2(end+1, 1)   = pdist2( XYexp(:,2)', sum(XYexp(:,3:end), 2)' , 'correlation');
    BestMatch.shift{end+1, 1}    = shift;
    BestMatch.Spectra{end+1, 1}  = XYexp;
    BestMatch.allMatch{end+1, 1} = Match;
    BestMatch.Ratio(end+1, 1)    = Match.Ratio(Id4Best(IP));
    BestMatch.nbrAdducts(end+1, 1) =   length(Id4Best);
    
end
BestMatch = struct2table(BestMatch);
