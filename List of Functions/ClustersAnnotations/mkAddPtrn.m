function AdductsPattern = mkAddPtrn(motif, Data, Adducts, NeutralLoss)
nmax = 3;

AdductsPattern.n = [];
AdductsPattern.Adducts = {};
AdductsPattern.NeutralLoss = {};
AdductsPattern.Z = [];
AdductsPattern.mf = {};
AdductsPattern.mi = [];
NeutralLoss(NeutralLoss.inESI ~= 1, :) = [];

for ii = 1:size(Adducts)
    
    for jj = 1:nmax
        
        mf = '';
        for kk = 1:jj
            mf = [mf motif];
        end
        mf = [mf Adducts.mf{ii}];
        
        switch Adducts.Charge(ii)
            case 1
                mf = [mf '(+1)'];
                
            case 2
                mf = [mf '(+2)'];
                
            case 3
                mf = [mf '(+3)'];
                
            case 4
                mf = [mf '(+4)'];
                
            case 5
                mf = [mf '(+5)'];
        end
        
        
        [mf2, mi, a] = DoCheckMf(mf, Data);
        if any(a.nbElements < 0)
            continue
        else
            
            AdductsPattern.n(end+1, 1) = kk;
            AdductsPattern.Adducts{end+1, 1} = Adducts.IonName{ii};
            AdductsPattern.NeutralLoss{end+1, 1} = '';
            AdductsPattern.Z(end+1, 1) = Adducts.Charge(ii);
            AdductsPattern.mf{end+1, 1} = mf2;
            AdductsPattern.mi(end+1, 1) = mi;
            
            for ll = 1:size(NeutralLoss, 1)
                [mf3, mi, a] = DoCheckMf([mf2 NeutralLoss.mf{ll}], Data);
                if any(a.nbElements < 0)
                    continue
                else
                    
                    AdductsPattern.n(end+1, 1) = kk;
                    AdductsPattern.Adducts{end+1, 1} = Adducts.IonName{ii};
                    AdductsPattern.NeutralLoss{end+1, 1} = NeutralLoss.IonName{ll};
                    AdductsPattern.Z(end+1, 1) = Adducts.Charge(ii);
                    AdductsPattern.mf{end+1, 1} = mf3;
                    AdductsPattern.mi(end+1, 1) = mi;
                end
                
            end
            
        end
        
        
    end
    
end
AdductsPattern = struct2table(AdductsPattern);

