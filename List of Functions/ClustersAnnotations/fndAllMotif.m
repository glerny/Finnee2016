function TblMotif = fndAllMotif(mf, Data, AdductsRules)
nmax = 3;
[~, ~, struct, Z] = DoCheckMf(mf, Data);
TblMotif.motif   = {};
TblMotif.mi   = [];
TblMotif.a    = {};
AdductsRules(AdductsRules.BaseAdduct ~= 1, :) = [];
I2Chck = find(AdductsRules.Charge == Z);
for ii = 1:length(I2Chck)
    a = struct;
    cAdducts = AdductsRules(I2Chck(ii), :);
    VN = cAdducts.Properties.VariableNames;
    
    for jj = 1:size(a, 1)
        I2Ele = find(strcmp(VN, a.Symbol{jj}));
        if length(I2Ele) == 1
            a.nbElements(jj) = a.nbElements(jj) - cAdducts.(VN{I2Ele});
        end
    end
    
    if any(a.nbElements < 0)
        continue
    else
        Motif = [];
        for jj = 1:size(a, 1)
            if a.nbElements(jj) > 0
                Motif = [Motif a.Symbol{jj} num2str(a.nbElements(jj))];
            end
        end
        
        [TblMotif.motif{end+1, 1}, TblMotif.mi(end+1, 1), TblMotif.a{end+1, 1}, Z]...
            = DoCheckMf(Motif, Data);
        
        for kk = 2:nmax
            
            if all(mod(a.nbElements, kk) == 0)
                Motif = [];
                
                for jj = 1:size(a, 1)
                    if a.nbElements(jj) > 0
                        Motif = [Motif a.Symbol{jj} num2str(a.nbElements(jj)/kk)];
                    end
                end
                
                [TblMotif.motif{end+1, 1}, TblMotif.mi(end+1, 1), TblMotif.a{end+1, 1}, Z]...
                    = DoCheckMf(Motif, Data);
            end
            
        end
    end
    
end
TblMotif = struct2table(TblMotif);
