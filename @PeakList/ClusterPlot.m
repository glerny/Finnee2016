function [BPP, TIP, nHACA] = ClusterPlot(obj, CluLvl, varargin)

% 1- Initialisation and options
narginchk(2, inf)

options = checkVarargin(CluLvl, varargin);
rpts    = length(obj.Path2Fin);
R       = [];
for jj = 1:rpts;
    X       = obj.AxisX{jj}.Data;
    
    % Build the matrix of profiles and initialize the HACA;
    MatOfPrf = zeros(length(X), length(obj.LstPIP{jj}));
    
    for ii = 1:length(obj.LstPIP{jj})
        cPIP   = obj.LstPIP{jj}{ii};
        IdSrt  = find(X == cPIP.x(1));
        IdEnd  = find(X == cPIP.x(end));
        MatOfPrf(IdSrt:IdEnd, ii) = cPIP.y;
        Id2PIP(ii,1) = ii;
        Id2PIP(ii,2) = max(cPIP.y);
        HACA{ii}.LstPIP    = ii;
        HACA{ii}.coco      = 1;
    end
    if isempty(R)
        R = corrcoef(MatOfPrf);
    else
        R = R + corrcoef(MatOfPrf);
    end
    BPP{jj}(:,1) = obj.AxisX{jj}.Data;
    BPP{jj}(:,4) = 0;
    TIP{jj}(:,1) = obj.AxisX{jj}.Data;
    TIP{jj}(:,4) = 0;
end

R = R/rpts;

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
delete(h)
assignin('base', 'HACA', HACA)
nHACA = {};

% Calculate Cluster's figure of merits
FOMPIP = [];
count  = 0;
for ii = 1:length(HACA)
    if length(HACA{ii}.LstPIP) == 1
        for jj = 1:rpts
            cPIP = obj.LstPIP{jj}{HACA{ii}.LstPIP};
            IdS  = cPIP.IdS;
            y    = cPIP.y;
            IdE  = IdS + length(y) - 1;
            BPP{jj}(IdS:IdE,2) = max(BPP{jj}(IdS:IdE,2), y);
            BPP{jj}(IdS:IdE,4) = max(BPP{jj}(IdS:IdE,4), y);
            TIP{jj}(IdS:IdE,2) = TIP{jj}(IdS:IdE,2) + y;
            TIP{jj}(IdS:IdE,4) = TIP{jj}(IdS:IdE,4) + y;
        end
    elseif length(HACA{ii}.LstPIP) > 1
        cFOM  = [];
        count = count + 1;
        
        FOMPIP(end+1, 1) = count;
        FOMPIP(end, 2)   = length(HACA{ii}.LstPIP);
        FOMPIP(end, 3)   = HACA{ii}.coco;
        nHACA{end+1}     = HACA{ii};
        
        provFOM = [];
        maxI    = 0;
        cFOM    = [];
        for jj = 1: length(HACA{ii}.LstPIP)
            for kk = 1:rpts
                cPIP = obj.LstPIP{kk}{HACA{ii}.LstPIP(jj)};
                IdS  = cPIP.IdS;
                y    = cPIP.y;
                IdE  = IdS + length(y) - 1;
                BPP{kk}(IdS:IdE,2) = max(BPP{kk}(IdS:IdE,2), y);
                BPP{kk}(IdS:IdE,3) = max(BPP{kk}(IdS:IdE,3), y);
                TIP{kk}(IdS:IdE,2) = TIP{kk}(IdS:IdE,2) + y;
                TIP{kk}(IdS:IdE,3) = TIP{kk}(IdS:IdE,3) + y;
                cFOM = [cFOM; cPIP.FOM];
            end
            
            if sum(cFOM(:,1)) > maxI
                maxI = sum(cFOM(:,1));
                FOMPIP(end, 4) = mean(cFOM(:, 4));
                FOMPIP(end, 5) = mean(cFOM(:, 9));
                FOMPIP(end, 6) = mean(cFOM(:, 3));
                FOMPIP(end, 7) = mean(cFOM(:, 1));
            end
        end
    end
end

InfoTrc = obj.BPP{1}.InfoTrc;

for ii = 1:rpts
    InfoTrc.Title  = 'Total Ion Profile';
    nTIP{ii} = Trace(InfoTrc, TIP{ii}(:, [1 3]));
    InfoTrc.Title = 'Base Peak profile';
    nBPP{ii} = Trace(InfoTrc, BPP{ii}(:, [1 3]));
    figure, subplot(2,1,1)
    plot(BPP{ii}(:,1), BPP{ii}(:,2), 'k')
    hold on
    plot(BPP{ii}(:,1), BPP{ii}(:,3), 'r')
    hold off
    subplot(2,1,2)
    plot(BPP{ii}(:,1), BPP{ii}(:,4), 'r')
end

% Cluster plot
gui4ClustPlot(obj, CluLvl, nHACA, FOMPIP, options, nBPP)
% 
% 
% 
% 
% 
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

