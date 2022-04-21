function obj = doHierarchicalClustering(obj, QCAnalysis, Filter, varargin)

%% 1- INITIALISATION AND OPTIONS;
ThrDifTime = 2;
MaxProf    = 2000;


switch QCAnalysis
    case 1
        FOM      = obj.QC.Method1.FOM;
        Profiles = obj.QC.Method1.Profiles;
        Cnew = {'IDFeature', 'nbrDetec', 'mean_M1', 'std_M1', 'mean_AccMass', ...
            'str_AccMass', 'mean_Area', 'std_Area', 'RSD_Area', 'Mean_PearsonCoef', 'Min_PearsonCoef'};
        FOM.Properties.VariableNames = Cnew;
        
    case 2
        FOM      = obj.QC.Method2.FOM;
        Profiles = obj.QC.Method2.Profiles;
        
    case 3
        FOM      = obj.QC.Method3.FOM;
        Profiles = obj.QC.Method3.Profiles;
        
    otherwise
        error
end
%% 2- Filter & Partition
AxisTime = obj.QC.Axis.AxisX.Data;
AxisMz   = obj.QC.Axis.AxisY;
FOM = sortrows(FOM,'IDFeature');
Filter   = dcdFilter(Filter);
Id2Rem   = FOM.nbrDetec < Filter.freq*max(FOM.nbrDetec)/100 | FOM.RSD_Area > Filter.RSD;
FOM(Id2Rem, :) = [];
ThrDifTime = ThrDifTime*mean(FOM.std_M1, 'omitnan');
FOM = sortrows(FOM, 'mean_M1');

Partitions = ones(size(FOM.IDFeature));

if size(Partitions, 1) <= MaxProf
    Cont = false;
    LstPart = 1;
else
    Cont = true;
end
i20 = [];
test = [0; diff(FOM.mean_M1)];
while Cont
    
    [Cut, IdCut] = max(test); 
    i20 = [i20, IdCut];
    
    if Cut < ThrDifTime, break, end
    Partitions(IdCut:end) = Partitions(IdCut:end)+1;
    
    LstPart = unique(Partitions);
    Cont = false;
    
    for ii = 1:length(LstPart)
        if sum(Partitions == LstPart(ii)) > MaxProf
            Cont = true;
            break
        end
    end
    test(i20) = 0;
end

if Cont
    Biggest = 0;
    LstPart = unique(Partitions);
    for ii = 1:length(LstPart)
        Biggest = max(Biggest, sum(Partitions == LstPart(ii)));
    end
    warning('Partitionning stop, biggest partition is %i',  Biggest)
else
    fprintf('\n\t Partitionning Done\n')
end

%% 3 - Hierarchical clustering
PartHACA = {};
for ii = 1:length(LstPart)
    id2H = Partitions == LstPart(ii);
    cFOM =  FOM(id2H, :);
    
    cProf  = zeros(length(AxisTime), length(cFOM.IDFeature));
    for jj = 1:length(cFOM.IDFeature)
        cXY = Profiles{cFOM.IDFeature(jj)};
        [~, ~, ib] = intersect(cXY(:,1), AxisTime, 'rows', 'legacy');
        cProf(ib, jj) = cXY(:,2);
    end
    
    cFOM = addvars(cFOM, (1:length(cFOM.IDFeature))', 'NewVariableNames','IDHA');
    Ix   = cFOM.IDHA;
    Z    = [];
    
    % DO the HACA
    text = ['creating cluster tree for partition ', num2str(ii), ' out of ', num2str(length(LstPart))];
    h = waitbar(0,'1','Name', text);
    
    while 1
        R    = corrcoef(cProf);
        R(isnan(R)) = 0;
        max_CC = max(max(triu(R,1)));
        
        if max_CC == 0
            break
        end
        
        if isvalid(h)
            waitbar(ii/length(R(:,1)), h, sprintf('%.3f', max_CC))
        else
            h = waitbar(0,'1','Name',text);
        end

        [linMin, colMin] = find(triu(R,1) == max_CC);
        
        intLin = Ix(linMin(1));
        intCol = Ix(colMin(1));
       
        Z(end+1, 1)     = min(intLin, intCol);
        Z(end,   2)     = max(intLin, intCol);
        Z(end,   3)     = max(Ix)+1;
        Z(end,   4)     = max_CC;
        cProf(:, end+1) = sum(cProf(:, [linMin(1), colMin(1)]), 2);
        cProf(:,  [linMin(1), colMin(1)]) = [];
        Ix(end+1) = max(Ix) + 1;
        Ix([linMin(1), colMin(1)]) = [];
        
        if size(cProf, 2) == 1
            break
        end
        
    end
    
    PartHACA{ii}.FOM = cFOM;
    PartHACA{ii}.Z   = Z;
    
    if isvalid(h)
        close(h)
    end
end

switch QCAnalysis
    case 1
        obj.QC.Method1.HACA = PartHACA;
        
    case 2
        obj.QC.Method2.HACA = PartHACA;
        
    case 3
        obj.QC.Method3.HACA = PartHACA;
        
    otherwise
        error
end

myMaster = obj;
save(fullfile(obj.Path, obj.Name), 'myMaster')


%% SUB FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CHECKVARARGIN
    function options = checkVarargin(CluLvl, varargin)
        % CHECKVARARGIN is used to check the input paramters
        % and create the options parameter.
        
        % 1- Verify validity of CluLvl
        if CluLvl > 1 || CluLvl < 0.8
            error('CluLvl should be between 0.8 and 1')
        end
        
        options.IntThreshold = 0;
        options.minPIPS      = 2;
        
        % 2- Decipher varargin
        input = @(x) find(strcmpi(varargin,x),1);
        tgtIx = input('IntThreshold');
        if ~isempty(tgtIx)
            options.IntThreshold = varargin{tgtIx +1};
        end
        
        tgtIx = input('minPIPS');
        if ~isempty(tgtIx)
            options.minPIPS = varargin{tgtIx +1};
        end
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

