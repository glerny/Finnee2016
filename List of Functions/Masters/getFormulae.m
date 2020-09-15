function [formulae, allData, AdductsList] = getFormulae(Isotopomer)

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
opts.DataRange = "A2:N36";

% Specify column names and types
opts.VariableNames = ["Ionname", "Charge", "C", "H", "O", "N", "S", "Na", "K", "Cu", "F", "Cl", "mass", "mf"];
opts.VariableTypes = ["string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "string"];
opts = setvaropts(opts, [1, 14], "WhitespaceRule", "preserve");
opts = setvaropts(opts, [1, 14], "EmptyFieldRule", "auto");

% Import the data
AdductsList = readtable("D:\OneDrive - Universidade do Porto\MSressources\POSADDUCTS.xlsx", opts, "UseExcel", false);

clear opts
Isotopomer  = sortrows(Isotopomer, 'Area', 'descend');

formulae = string([]);
for ii = 1:size(AdductsList, 1)
    cAdducts = AdductsList(ii, :);
    z = cAdducts.Charge;
    mAdducts = cAdducts.mass;
    
    for n = nN
        
        M  = (z*Isotopomer.AccurateMass(1) - mAdducts)/nN(n);
        st = Isotopomer.AccurateMass(1)*ste/1000000;
        
        if M > 1
            cData = webread(url2, 'monoisotopicMass', M,'mfRange', ...
                mfRange1, 'massRange', st);
            if cData.numberResults > 0
                C = struct2table(cData.results);
                formulae = vertcat(formulae, string(C.mf));
                
            end
        end
        
    end
end

allData = {};
formulae = unique(formulae);
length(formulae)
for ii = 1:length(formulae)
    ii
    
    data = webread(url1, 'mf', formulae{ii});
    M = data.em;
    El = struct2table(data.parts.ea);
    
    MStheo = [];
    for jj = 1:size(AdductsList, 1)
        for kk = nN
            
            mf = [];
            if (nN(kk)*M + AdductsList.mass(jj))/AdductsList.Charge(jj) >= ...
                    min( (Isotopomer.AccurateMass - st*Isotopomer.AccurateMass/1000000)) && ...
                    (nN(kk)*M + AdductsList.mass(jj))/AdductsList.Charge(jj) <= ...
                    max( (Isotopomer.AccurateMass + st*Isotopomer.AccurateMass/1000000))
                
                for ll = 1:nN(kk)
                    mf = [mf, char(formulae{ii})];
                end
                
                if any(table2array(AdductsList(jj, 3:12)) <0)
                    data = webread(url1, 'mf', mf);
                    mf = [mf, char(AdductsList.mf(jj))];
                    
                    % NEED TO DO THE TEST TO AVOID THE TRY CATCH
                elseif char(AdductsList.mf(jj)) == ''''
                else
                    mf = [mf, char(AdductsList.mf(jj))];
                end
                
                switch AdductsList.Charge(jj)
                    case 1
                        mf = [mf, '(+1)'];
                        
                    case 2
                        mf = [mf, '(+2)'];
                        
                    case 3
                        mf = [mf, '(+3)'];
                        
                end
                
                try
                    data = webread(url1, 'mf', mf, 'resolution', '0.0000001',...
                        'isotopomers', 'jcamp,xy');
                    
                    MS = str2num( regexprep(data.xy, '[\n\r]+',';'));
                    MS(MS(:,2) < 0.5, :) = [];
                    MS(:,3) = ii;
                    MS(:,4) = jj;
                    MS(:,5) = kk;
                    
                    MStheo = [MStheo; MS];
                catch ME
                    disp(ME)
                    disp(mf)
                end
            end
            
        end
    end
    
    disp('testhere')
    % Test the formulae for best fit
    MSexp = [Isotopomer.AccurateMass, Isotopomer.Area];
    MSexp(:,5) = 0;
    MSmer = [MStheo; MSexp];
    MSmer = sortrows(MSmer, 1);
    allData{ii} = MSmer;
    assignin('base', 'allData', allData)
    
end




disp('testhere')








end

