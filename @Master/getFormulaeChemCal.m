function obj = getFormulaeChemCal(obj, dmzst)

url1 = 'http://www.chemcalc.org/chemcalc/mf';
url2 = 'http://www.chemcalc.org/chemcalc/em';
mfRange = 'C0-100H0-100N0-10O0-10F0-3Cl0-3Br0-1S0-3';

obj.HACA.MonoisotopicMasses.Predicted = {};
for ii = 1:size(obj.HACA.MonoisotopicMasses.Summary, 1)
    TestForm = [];
    mass = obj.HACA.MonoisotopicMasses.Summary(ii, 3);
    IEnv = obj.HACA.MonoisotopicMasses.ListIons{ii};
    TestForm(1, :) = [mass, obj.HACA.MonoisotopicMasses.Summary(ii, 4)];
    Data4Mw = IEnv.IsotopicEnvelope(:,6:10);
    
    M0 = sum(Data4Mw(:,4));
    sM0 = sqrt(sum(Data4Mw(:,5).^2));
    Data4Mw(:,6) = Data4Mw(:,2).*Data4Mw(:,4);
    Data4Mw(:,7) = Data4Mw(:,6).*sqrt((Data4Mw(:,3)./Data4Mw(:,2)).^2 + ...
        (Data4Mw(:,5)./Data4Mw(:,4)).^2); 
    sM1 = sqrt(sum(Data4Mw(:,7).^2));
    M1 = sum(Data4Mw(:,2).*Data4Mw(:,4));
    TestForm(2, 1) = M1/M0;
    TestForm(2, 2) = M1/M0*sqrt((sM0/M0)^2 + (sM1/M1)^2);
    formulae = {};
    
    massRange = mass*dmzst/1000000;
    data = webread(url2, 'monoisotopicMass', mass,'mfRange', mfRange, ...
        'massRange', massRange);
    
    if data.numberResults > 0
        for jj = 1:data.numberResults
            data2 = webread(url1, 'mf', data.results(jj).mf, 'isotopomers', 'jcamp,xy');
            TestForm(1, jj+2) = data2.em;
            TestForm(2, jj+2) = data2.mw;
            formulae{jj} = data2.mf;
        end
    end
    obj.HACA.MonoisotopicMasses.Predicted{ii}.data = TestForm;
    obj.HACA.MonoisotopicMasses.Predicted{ii}.formulae = formulae;
end


myMaster = obj;
try
    save(fullfile(obj.Path, obj.Name), 'myMaster')
catch
    uisave('myMaster')
end

end
