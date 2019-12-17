function obj = ClusterPlot(obj, CluLvl, varargin)

% 1- Initialisation and options
narginchk(2, inf)

options = checkVarargin(CluLvl, varargin);


X       = obj.ClusteredQC.TimeAxis;
% Build the matrix of profiles and initialize the HACA;
MatOfPrf = zeros(length(X), length(obj.QCFiles));

for ii = 1:length(obj.ClusteredQC.Profiles)
    MatOfPrf(ismember(X, obj.ClusteredQC.Profiles{ii}(:,1)), ii) = ...
        obj.ClusteredQC.Profiles{ii}(:,2);
    Id2PIP(ii,1) = ii;
    Id2PIP(ii,2) = max(obj.ClusteredQC.Profiles{ii}(:,2));
    HACA{ii}.LstPIP    = ii;
    HACA{ii}.coco      = 1;
end
R = corrcoef(MatOfPrf);

% DO the HACA
h = waitbar(0,'1','Name','creating cluster tree...');
for ii = 2:length(R(:,1))
    max_CC = max(max(triu(R,1)));
    if isvalid(h)
        waitbar(ii/length(R(:,1)), h, sprintf('%.3f => %.3f', max_CC, CluLvl))
    else
        h = waitbar(0,'1','Name','creating cluster tree...');
    end
    if max_CC < CluLvl, break; end
    [linMin, colMin] = find(triu(R,1) == max_CC);
    intLin = Id2PIP(linMin(1), 2);
    intCol = Id2PIP(colMin(1), 2);
    if intLin > intCol  % take only one
        indMax = linMin(1); indMin = colMin(1);
    else
        indMax = colMin(1); indMin = linMin(1);
    end
    Id2PIP1 = Id2PIP(indMax, 1);
    Id2PIP2 = Id2PIP(indMin, 1);
    
    HACA{Id2PIP1}.coco   = max_CC;
    HACA{Id2PIP2}.coco   = 0;
    HACA{Id2PIP1}.LstPIP = [HACA{Id2PIP1}.LstPIP, HACA{Id2PIP2}.LstPIP];
    HACA{Id2PIP2}.LstPIP = [];
    R(:, indMin)         = [];
    R(indMin, :)         = [];
    Id2PIP(indMin, :)    = [];
end
nHACA = {};

% Calculate Cluster's figure of merits
FOMPIP   = [];
count    = 0;
Summary  = [];
Profiles = {};
for ii = 1:length(HACA)
    if length(HACA{ii}.LstPIP) >= 1
        cFOM  = [];
        count = count + 1;
        
        FOMPIP(end+1, 1) = count;
        FOMPIP(end, 2)   = length(HACA{ii}.LstPIP);
        FOMPIP(end, 3)   = HACA{ii}.coco;
        nHACA{end+1}     = HACA{ii};
        
        cP = X;
        for jj = 1:length(HACA{ii}.LstPIP)
            cP(:, jj+1) = MatOfPrf(:, HACA{ii}.LstPIP(jj));
        end
        Profiles{count} = cP;
        XX = corrcoef([mean(cP(:, 2:end), 2), cP(:, 2:end)]);
        XX = XX(2:end, 1);
        
        for jj = 1:length(HACA{ii}.LstPIP) 
            cC = obj.ClusteredQC.Cluster{HACA{ii}.LstPIP(jj)};
            Summary(end+1, 1) = HACA{ii}.LstPIP(jj);
            Summary(end,   2) = count;
            Summary(end,   3) = size(cC, 1);
            Summary(end,   4) = mean(cC(:,3));
            Summary(end,   5) = std(cC(:,3));
            Summary(end,   6) = mean(cC(:,6));
            Summary(end,   7) = std(cC(:,6));
            Summary(end,   8) = mean(cC(:,8));
            Summary(end,   9) = std(cC(:,8));
            Summary(end,  10) = XX(jj);
        end
    end
end
obj.HACA.FOM      = FOMPIP;
obj.HACA.Summary  = Summary;
obj.HACA.Profiles = Profiles;
uisave('obj', 'myMaster.mat')

%

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



end

