function [Results, AdductsList] = getFormulae_2(Isotopomer)

%% 1- Intro
polarity = '+';
url1     = 'http://www.chemcalc.org/chemcalc/mf';
url2     = 'http://www.chemcalc.org/chemcalc/em';
mfRange1 = 'C0-100 H0-100 N0-10 O0-10';
mfRange2 = 'C0-100 H0-100 N0-10 O0-10 S0-5 Cl0-5 F0-5 Br0-5';
ste      = 5; %ppm
nN       = 1:3;

opts = spreadsheetImportOptions("NumVariables", 14);

% Specify sheet and range
opts.Sheet = "Adducts";
% opts.DataRange = "A2:N36";

% Specify column names and types
opts.VariableNames = ["Ionname", "Charge", "C", "H", "O", "N", "S", "Na", "K", "Cu", "F", "Cl", "mass", "mf"];
opts.VariableTypes = ["string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "string"];
opts = setvaropts(opts, [1, 14], "WhitespaceRule", "preserve");
opts = setvaropts(opts, [1, 14], "EmptyFieldRule", "auto");

% Import the data
AdductsList = readtable("D:\OneDrive - Universidade do Porto\MSressources\POSADDUCTS.xlsx", opts, "UseExcel", false);
AdductsList(1, :) = [];

clear opts
Isotopomer  = sortrows(Isotopomer, 'Area', 'descend');

%% 2- Get possible mf for each possible M

formulae = string([]);
Results.Formulae = {};

for ii = 1:size(AdductsList, 1)
    cAdducts = AdductsList(ii, :);
    z = cAdducts.Charge;
    mAdducts = cAdducts.mass;
    
    for n = nN
        
        M  = (z*Isotopomer.AccurateMass(1) - mAdducts)/nN(n);
        st = M*ste/1000000;
        
        if M > 1
            cData = webread(url2, 'monoisotopicMass', M,'mfRange', ...
                mfRange1, 'massRange', st);
            if cData.numberResults > 0
                C = struct2table(cData.results);
                % C is the table of possible formulae for M calucate from
                % mz
                
                for jj = 1:size(C, 1)
                    
                    %% 3- Based on this M defined with sepecific adduct and repetition of motif find and match all possible adducts
                    if size(C, 1) == 1
                        Adducts = GetAdducts(M, C.mf, AdductsList);
                        
                    else
                        Adducts = GetAdducts(M, C.mf{jj}, AdductsList);
                        
                    end
                    
                    % Check for the adducts existence,
                    ix = Adducts.AdductMass >= min(Isotopomer.AccurateMass - ste*Isotopomer.AccurateMass/1000000) & ...
                        Adducts.AdductMass <= max(Isotopomer.AccurateMass + ste*Isotopomer.AccurateMass/1000000);
                    Adducts = Adducts(ix, :);
                    Link2Adducts = cell(size(Isotopomer, 1), 1);
                    
                    for kk = 1:size(Adducts, 1)
                        [STE, iSTE] = min((abs(Isotopomer.AccurateMass - Adducts.AdductMass(kk)))/Adducts.AdductMass(kk));
                        if STE*1000000 < ste
                            Link2Adducts{iSTE} = [Link2Adducts{iSTE}, kk];
                        end
                    end
                    Link2Adducts{1} = find(Adducts.AdductID == ii & Adducts.AdductnM == n);
                    ix = cell2mat(cellfun(@size,Link2Adducts,'uni',false));
                    
                    ix = find(ix(:,1) >= 1);
                    
                    MStheo = [];
                    for kk = 1:length(ix)
                        if  size(Link2Adducts{ix(kk)}, 2) > 1
                            testAdducts = Adducts(Link2Adducts{ix(kk)}, :);%!!!'DOTHETEST'
                        end
                        
                        adc = Link2Adducts{ix(kk)};
                        adc = adc(1); %!!!'DOTHETEST'
                        
                        MS = GetMS(Adducts.Adductmf{adc});
                        
                        if isempty(MS)
                            continue;
                        end
                        
                        MS(:,3)     = adc;
                        MS(:, 4:6) = nan;
                        MStheo = [MStheo; MS];
                    end
                    
                    MS = MStheo;
                    
                    if isempty(MS)
                        continue;
                    end
                    for kk = 1:size(MS, 1)
                        [STE, iSTE] = min((abs(Isotopomer.AccurateMass -  MS(kk, 1))/MS(kk, 1)));
                        if STE*1000000 < ste
                            MS(kk,4:5) = [Isotopomer.AccurateMass(iSTE), Isotopomer.Area(iSTE)];
                            MS(kk,6)   = iSTE;
                        end
                        
                    end
                    
                    if isnan(MS(1, 4))
                        disp('')
                    end
                    
                    % Check for double input
                    List = unique(MS(~isnan(MS(:,6)),6));
                    for kk = 1:length(List)
                        if sum(MS(:,6) == List(kk)) > 1
                            IxX = find(MS(:,6) == List(kk));
                            
                            %check if single or multiple Adducts
                            if length(unique(MS(IxX, 3))) == 1
                                [~, ML] = max(MS(IxX, 2));
                                IxX(ML) = [];
                                MS(IxX, 4:6) = NaN;
                            else
                                listAdducts = unique(MS(IxX, 3));
                                
                                % Find most aboundant adducts
                                I2Max = 0;
                                V2Max = 0;
                                
                                for ll = 1:length(listAdducts)
                                    if MS(MS(:,3) == listAdducts(ll) & MS(:,2) == 100, 5) >= V2Max
                                        I2Max = listAdducts(ll);
                                        
                                    end
                                end
                                IxX(MS(IxX, 3) == I2Max) = [];
                                MS(IxX, 4:6) = NaN;
                                
                                IxX = find(MS(:,6) == List(kk));
                                if length(IxX) > 1
                                    [~, ML] = max(MS(IxX, 2));
                                    IxX(ML) = [];
                                    MS(IxX, 4:6) = NaN;
                                    
                                end
                                
                            end
                        end
                    end
                    
                    
                    fAdd = MS(1, 3);
                    if ~isnan(MS(1, 3))
                        
                        MSfAdd = MS(MS(:,3) == fAdd, :);
                        [~, ord2] = sort(MSfAdd(:,2), 'descend','MissingPlacement','last');
                        [~, ord4] = sort(MSfAdd(:,5), 'descend','MissingPlacement','last');
                        
                        % Test.Formulae
                        Test.mzBP_theo      = MSfAdd(ord2(1), 1);
                        Test.mzBP_exp       = MSfAdd(ord4(1), 4);
                        Test.error1         = mean(MSfAdd(:,1) - MSfAdd(:,4), 'omitnan');
                        Test.error2         = mean((MSfAdd(:,1) - MSfAdd(:,4)).*(MSfAdd(:,2).*MSfAdd(:,5)), 'omitnan')/...
                            mean(MSfAdd(:,2).*MSfAdd(:,5), 'omitnan');
                        Test.DeltaMZ21_theo = MSfAdd(ord2(1), 1) -  MSfAdd(ord2(2), 1);
                        Test.DeltaMZ21_exp  = MSfAdd(ord4(1), 4) -  MSfAdd(ord4(2), 4);
                        Test.RatioMZ21_theo = MSfAdd(ord2(2), 2) / MSfAdd(ord2(1), 2);
                        Test.RatioMZ21_exp  = MSfAdd(ord4(2), 5) / MSfAdd(ord4(1), 5);
                        Test.Mw_theo        = sum(MSfAdd(:,1).*MSfAdd(:,2), 'omitnan')/...
                            sum(MSfAdd(:,2), 'omitnan');
                        Test.Mw_exp = sum(MSfAdd(:,4).*MSfAdd(:,5), 'omitnan')/...
                            sum(MSfAdd(:,5), 'omitnan');
                        
                        
                        if (Test.mzBP_theo - Test.mzBP_exp)/Test.mzBP_exp*1000000 < ste && ...
                                (Test.DeltaMZ21_theo - Test.DeltaMZ21_exp) < 0.01 && ...
                                2*abs(Test.RatioMZ21_theo - Test.RatioMZ21_exp)/(Test.RatioMZ21_theo + Test.RatioMZ21_exp) < 0.5
                            % Keep it!
                            
                            if size(C, 1) == 1
                                Results.Formulae{end+1}.mf  = C.mf;
                            else
                                Results.Formulae{end+1}.mf  = C.mf{jj};
                            end
                            Results.Formulae{end}.nM    = nN(n);
                            Results.Formulae{end}.nAdd  = length(unique(MS(:,3)));
                            Results.Formulae{end}.Match = sum(~isnan(MS(:,4)));
                            Results.Formulae{end}.PercMatch = sum(~isnan(MS(:,4)))/length(Isotopomer.AccurateMass)*100;
                            Results.Formulae{end}.Test = Test;
                            Results.Formulae{end}.MS = MS;
                            Results.Formulae{end}.Adducts = Adducts(unique(MS(:, 3)), :);
                            
                            % Based on MS check for other adducts...
                            
                        end
                    else
                        disp('checkme')
                    end
                end
                
            end
        end
        
    end
end
%% - NESTED FUNCTIONS

    function MS = GetAdducts(M, mf, Adducts)
        MS = {};
        
        for GAii = 1:size(Adducts, 1)
            for GAn = nN
                
                cmf = [];
                for GMjj = 1:GAn
                    cmf = [cmf, char(mf)];
                end
                
                cmf = [cmf, char(Adducts.mf{GAii})];
                cmf = VerifyFormulae(cmf);
                
                if isempty(cmf), continue; end
                
                
                switch Adducts.Charge(GAii)
                    
                    case 1
                        cmf = [cmf, '(+1)'];
                        
                    case 2
                        cmf = [cmf, '(+2)'];
                        
                    case 3
                        cmf = [cmf, '(+3)'];
                end
                
                mz = (GAn*M + Adducts.mass(GAii))/Adducts.Charge(GAii);
                
                if isempty(MS)
                    MS.ID(1)           = 1;
                    MS.AdductID(1)     = GAii;
                    MS.AdductnM(1)     = GAn;
                    MS.AdductName{1}   = Adducts.Ionname{GAii};
                    MS.AdductCharge(1) = Adducts.Charge(GAii);
                    MS.AdductMass(1)   = mz;
                    MS.Adductmf{1}     = cmf;
                else
                    MS.ID(end+1, 1)           = MS.ID(end) +1;
                    MS.AdductID(end+1, 1)     = GAii;
                    MS.AdductnM(end+1, 1)     = GAn;
                    MS.AdductName{end+1, 1}   = Adducts.Ionname{GAii};
                    MS.AdductCharge(end+1, 1) = Adducts.Charge(GAii);
                    MS.AdductMass(end+1, 1)   = mz;
                    MS.Adductmf{end+1, 1}     = cmf;
                end
            end
            
        end
        
        MS = struct2table(MS);
        
    end

    function cmp = VerifyFormulae(cmp)
    end

    function MS = GetMS(mf)
        thresMS   = 0.1;
        url     = 'http://www.chemcalc.org/chemcalc/mf';
        
        try
            data = webread(url, 'mf', mf, 'resolution', '0.0000001',...
                'isotopomers', 'jcamp,xy');
            
            MS = str2num( regexprep(data.xy, '[\n\r]+',';'));
            MS(MS(:,2) < thresMS, :) = [];
            
        catch ME
            
            if strcmp(ME.identifier, 'MATLAB:nonExistentField')
                MS = [];
            else
                
                disp(ME)
                disp(mf)
            end
        end
        
    end


end

