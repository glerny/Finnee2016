function [BPP, TIP] = ClusterPlot(obj, CluLvl, varargin)

% 1- Initialisation and options
narginchk(2, inf)

options = checkVarargin(CluLvl, varargin);
X       = obj.AxisX.Data;
infoX   = obj.AxisX.InfoAxis;
infoZ   = obj.AxisZ.InfoAxis;

% Build the matrix of profiles and initialize the HACA;
MatOfPrf = zeros(length(X), length(obj.LstPIP));

for ii = 1:length(obj.LstPIP)
    cPIP   = obj.LstPIP{ii};
    IdSrt  = find(X == cPIP.x(1));
    IdEnd  = find(X == cPIP.x(end));
    MatOfPrf(IdSrt:IdEnd, ii) = cPIP.y;
    Id2PIP(ii,1) = ii;
    Id2PIP(ii,2) = max(cPIP.y);
    HACA{ii}.LstPIP    = ii;
    HACA{ii}.coco      = 1;
end
R = corrcoef(MatOfPrf);

% DO the HACA
for ii = 2:length(R(:,1))
    max_CC = max(max(triu(R,1)));
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

% Calculate Cluster's figure of merits
FOMPIP      = [];
forPrf(:,1) = X;
forPrf(:,3) = 0;
for ii = 1:length(HACA)
    if length(HACA{ii}.LstPIP) > 1
        FOMPIP(end+1, 1) = ii;
        FOMPIP(end, 2)   = length(HACA{ii}.LstPIP);
        FOMPIP(end, 3)   = HACA{ii}.coco;
        provFOM = [];
        maxI    = 0;
        for jj = 1: length(HACA{ii}.LstPIP)
            cPIP     = obj.LstPIP{HACA{ii}.LstPIP(jj)};
            provFOM  = [provFOM; cPIP.FOM];
            IdSrt    = find(X == cPIP.x(1));
            IdEnd    = find(X == cPIP.x(end));
            prf      = zeros(length(X), 1);
            
            prf(IdSrt:IdEnd) = cPIP.y;
            forPrf(:,2)      = forPrf(:,2) + prf;
            forPrf(:,3)      = max(forPrf(:,3), prf);
            maxI             = max(maxI, max(cPIP.y));
        end
        FOMPIP(end, 4) = mean( provFOM(:,4 ));
        [~, Id]        = max(provFOM(:,3 ));
        FOMPIP(end, 5) = provFOM(Id,9 );
        FOMPIP(end, 6) = sum(provFOM(:,3));
        FOMPIP(end, 7) = maxI;
    end
end
InfoTrc.Title  = 'Total Ion Profile';
InfoTrc.FT     = sprintf('From HACA @ %s', obj.Path2PkL);
InfoTrc.TT     ='SEP';
InfoTrc.AxisX  = Axis(infoX);
InfoTrc.AxisY  = Axis(infoZ); 
InfoTrc.Loc    = 'inTrace';
InfoTrc.AdiPrm = {};

TIP = Trace(InfoTrc, forPrf(:, [1 2]));
InfoTrc.Title = 'Base Peak profile';
BPP = Trace(InfoTrc, forPrf(:, [1 3]));

% Cluster plot
gui4ClustPlot(obj, CluLvl, HACA, FOMPIP, options, BPP, obj.LstPIP)







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

