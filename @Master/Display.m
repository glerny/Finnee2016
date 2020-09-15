function obj = Display(obj, tgtQC, mode, Filter, TgtFolder)


if nargin == 1
    thtQC = 3;
    mode = 'QC';
    Filter = 'RSD15;freq75;';
    TgtFolder =  uigetdir();
    
elseif nargin == 2
    mode = 'QC';
    Filter = 'RSD15;freq75;';
    TgtFolder =  uigetdir();
    
elseif nargin == 3
    Filter = 'RSD15;freq75;';
    TgtFolder =  uigetdir();
    
elseif nargin == 4
    TgtFolder =  uigetdir();
    
elseif nargin > 5 || nargin == 0 
    error('')
end

%% 1- Parameters and calculating the ROIs(to be changed and add to varadin)
switch tgtQC
    case 1
        ProfQC = obj.QC.Method1.Profiles;
        FOMQC  = obj.QC.Method1.FOM;
        
    case 2
        ProfQC = obj.QC.Method2.Profiles;
        FOMQC  = obj.QC.Method2.FOM;
        
    case 3
        ProfQC = obj.QC.Method3.Profiles;
        FOMQC  = obj.QC.Method3.FOM;
        
    otherwise
        error
end

fileName = fullfile(TgtFolder, 'ReadMe.txt');
fileID   = fopen(fileName,  'wt');

fprintf(fileID, '%s\n', datetime);
fprintf(fileID, '\nINFORMATION MASTER OBJECT\n');
fprintf(fileID, '\tPath: %s\n', obj.Path);
fprintf(fileID, '\tName: %s\n', obj.Name);
fprintf(fileID, '\tTarget Analysis: %i\n', tgtQC);
fprintf(fileID, '\tFilter: %s\n', Filter);
repl = max(size(obj.QC.Files));
fprintf(fileID, '\tNumber of files: %i\n', repl);
for ii = 1:repl
    fprintf(fileID, '\t\t%i. %s\n', ii, obj.QC.Files{ii});
end


Filter   = dcdFilter(Filter);
IF = FOMQC.RSD_Area <= Filter.RSD & FOMQC.nbrDetec >= Filter.freq*max(FOMQC.nbrDetec)/100;

fprintf(fileID, '\nDETECTED FEATURES \n');
fprintf(fileID, '\tTotal number of features: %i\n', size(IF, 1));
fprintf(fileID, '\t\t %.2f < Time(min) < %.2f \n', min(FOMQC.mean_M1), max(FOMQC.mean_M1));
fprintf(fileID, '\t\t %.4f < mz < %.4f \n', min(FOMQC.mean_AccMass), max(FOMQC.mean_AccMass));
fprintf(fileID, '\t\t %.0f < Area (a.u.) < %.0f \n', min(FOMQC.mean_Area), max(FOMQC.mean_Area));

fprintf(fileID, '\n\tFeatures after filtering: %i\n', sum(IF));
fprintf(fileID, '\t\t %.2f < Time(min) < %.2f \n', min(FOMQC.mean_M1(IF)), max(FOMQC.mean_M1(IF)));
fprintf(fileID, '\t\t %.4f < mz < %.4f \n', min(FOMQC.mean_AccMass(IF)), max(FOMQC.mean_AccMass(IF)));
fprintf(fileID, '\t\t %.0f < Area (a.u.) < %.0f \n', min(FOMQC.mean_Area(IF)), max(FOMQC.mean_Area(IF)));
fclose(fileID)

disp('doit')
IF = find(IF);
for ii = 1:length(IF)
    cProf = ProfQC{IF(ii)};
    cFOM = FOMQC(IF(ii), :);
    figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off')
    plot(cProf(:,1), cProf(:,2), 'k')
    title(sprintf('Feature #%i: Target Time:%.2f min| Target Mass:%.4f m/z',...
        cFOM.IDFeature, cFOM.mean_M1, cFOM.mean_AccMass))
    xlabel('time /min')
    ylabel('Intensity /a.u.')
    is = find(cProf(:,1) >= cFOM.Lm1, 1, 'first');
    text(cProf(is,1),cProf(is,2),'\downarrow t_i','FontSize', ...
        14,'FontWeight','bold', 'Color', 'red', 'VerticalAlignment', 'bottom')
    
    ie = find(cProf(:,1) <= cFOM.Lm2, 1, 'last');
    text(cProf(ie,1),cProf(ie,2),'\downarrow t_f','FontSize', ...
        14,'FontWeight','bold', 'Color', 'red', 'VerticalAlignment', 'bottom')
    
    saveas(gcf,  fullfile(TgtFolder, ['figure', num2str(ii) '.jpg']))
    close(gcf)
    
end

    function FltOut =  dcdFilter(FltIn)
        Stoppers = strfind(FltIn, ';');
        if length(Stoppers) ~= 2
            error('Filter incorrect')
        end
        
        IxF = strfind(FltIn, 'RSD');
        FltOut.RSD  = str2double(FltIn(IxF + 3:Stoppers(1)-1));
        IxF = strfind(FltIn, 'freq');
        FltOut.freq = str2double(FltIn(IxF + 4:Stoppers(2)-1));
    end

end
