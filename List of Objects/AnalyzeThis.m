%% DESCRIPTION
%
%% Copyright 
% BSD 3-Clause License
% Copyright 2016-2017 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


classdef AnalyzeThis
    
    properties
        %% Options for fileInfo
        DataIn              = []
        % First column: axe, second column dependent values
        
        BaselineMethod  	= 'ArPLS:10E6:0.001'
        % Possible methods:
        %   'None'      	No additional parameter
        % 	'PF:p1'         Polynomial fitting, p1 is the degree of the
        %                   polynomial (p1=0: constant; p1=1: linear...) 
        %   'ArPLS:p1:p2'   Baseline correction using asymmetrically 
        %                   reweighted penalized least squares smoothing, 
        %                   (doi: 10.1039/c4an01061b), 
        %                   p1: lambda(10E2-10E9); p2: alpha
        
        SmoothingMethod     = 'None'
        % Possible methods:
        %   'None'          No additional parameter
        
        PeakPickingMethod 	= 'LmMm:2:10'
        % Possible methods:
        %   'LmMm:p1:p2'    Peak detected as a local Maximum between 2
        %                   local minima, p1 is the number of neighbourgh 
        %                   used to detect local maximum or minimum (1:10),
        %                   p2 is an intensity threshold below witch peaks
        %                   will not be calculated
        
        DeconvolutionMethod = struct('name', 'None')
        % Possible methods:
        %   'None'          No additional parameter
        
        
    end
    
    properties (Dependent)
        Baseline
        PeakList
        Deconvolution
    end
    
    methods
        function obj = AnalyzeThis(dataIn, varargin)
             if nargin  <1
                error('error')
            end
            obj.DataIn = dataIn;
            
            if nargin  < 2
                return
            end
            
            input = @(x) find(strcmpi(varargin,x),1);
            
            tgtIx = input('baseline');
            if ~isempty(tgtIx)
                obj.BaselineMethod = varargin{tgtIx +1};
            end
            
            tgtIx = input('peakpicking');
            if ~isempty(tgtIx)
                obj.PeakPickingMethod = varargin{tgtIx +1};
            end
        end
        
        function baseline = get.Baseline(obj)
            basDef =  regexp(obj.BaselineMethod, ':', 'split');
            
            switch basDef{1}
                case 'None'
                    [m, ~] = size(obj.DataIn);
                    XY = obj.DataIn;
                    baseline.bckgPts = false(m,1);
                    baseline.noise = 4*std(XY(:,2));
                    baseline.vals = zeros(m, 1);
                    
                case 'PF'
                    XY = obj.DataIn;
                    [z, bslPts] = doPF(XY, str2double(basDef{2}));
                    baseline.bckgPts = bslPts;
                    baseline.noise = 4*std(XY(bslPts,2) - z(bslPts));
                    baseline.vals = z;
                    
                case 'ArPLS'
                    XY = obj.DataIn;
                    lambda = str2double(basDef{2});
                    ratio = str2double(basDef{3});
                    [m, ~] = size(obj.DataIn);
                    
                    iNZ = XY(:,2) ~= 0;
                    [z, bslPts] = doArPLS(XY(iNZ, 2), lambda, ratio);
                   baseline.bckgPts = false(m,1);
                    baseline.bckgPts(iNZ) = bslPts;
                    baseline.vals = zeros(m, 1);
                    baseline.vals(iNZ) = z;
                    baseline.noise = 4*std(nonzeros(XY(baseline.bckgPts,2) ...
                        - baseline.vals(baseline.bckgPts)));
                    
                case 'ArPLS2'
                    XY = obj.DataIn;
                    lambda = str2double(basDef{2});
                    [m, ~] = size(obj.DataIn);
                    
                    iNZ = XY(:,2) ~= 0;
                    [z, bslPts] = doArPLS2(XY(iNZ, 2), lambda);
                    baseline.bckgPts = false(m,1);
                    baseline.bckgPts(iNZ) = bslPts;
                    baseline.vals = zeros(m, 1);
                    baseline.vals(iNZ) = z;
                    baseline.noise = 4*std(nonzeros(XY(baseline.bckgPts,2) ...
                        - baseline.vals(baseline.bckgPts)));

            end
        end
        
        function peakList = get.PeakList(obj)
            basDef =  regexp(obj.PeakPickingMethod, ':', 'split');
            
            switch basDef{1}
                case 'LmMm'
                    
                    if str2double(basDef{3}) == Inf
                        thrs = 3*obj.Baseline.noise;
                    else
                        thrs = str2double(basDef{3});
                    end
             
                    yOri = obj.DataIn(:,2) - obj.Baseline.vals;
                    X = spread(yOri, str2double(basDef{2}));
                    
                    locMax = find(X(:,((1+end)/2)) == max(X, [], 2) & ...
                        max(X, [], 2) >= thrs);
                    locMin =  find(X(:,((1+end)/2)) == min(X, [], 2));
                    IdPair = zeros(length(locMax), 2);
                    for jj = 1:length(locMax)
                        IX = find(locMin < locMax(jj), 1, 'last');
                        if isempty(IX), IX = 1; end
                        IdPair(jj, 1) = locMin(IX);
                        
                        IX = find(locMin > locMax(jj), 1, 'first');
                        if isempty(IX), IX = length(locMin); end
                        IdPair(jj, 2) = locMin(IX);
                        
                    end
                    
                    [~, m, ~] = unique(IdPair(:,1));
                    IdPair = IdPair(m, :);
                    I2R = IdPair(:,2) - IdPair(:,1) < 2;
                    IdPair(I2R, :) = [];
                    
                    FOM = zeros(length(IdPair(:,2)), 8);
                    
                    peakList.headings = {'1:PeakStart', '2:PeakEnd', ...
                        '3:M0', '4:M1', '5:M2', '6:M3', ...
                        '7:Imax', '8:@Imax'};
                    for jj = 1:length(IdPair(:,2))
                        xc = obj.DataIn(IdPair(jj,1): IdPair(jj,2), 1);
                        yc = yOri(IdPair(jj,1): IdPair(jj,2));
                        FOM(jj,1) = xc(1);
                        FOM(jj,2) = xc(end);
                        FOM(jj,3) = trapz(xc, yc);
                        FOM(jj,4) = trapz(xc, xc.*yc)/FOM(jj,3);
                        FOM(jj,5) = trapz(xc, (xc-FOM(jj,4)).^2.*yc)/FOM(jj,3);
                        FOM(jj,6) = trapz(xc, (xc-FOM(jj,4)).^3.*yc)/FOM(jj,3);
                        [FOM(jj,7), indAtMax] = max(yc);
                        FOM(jj,8) = xc(indAtMax);
                    end
                    
            end
            peakList.data = FOM;
        end
        
        function plotAnalysis(obj)
            
                figure('Name', 'Result AnalyzeThis');
                hold on
                lim2Plot(:,1) = ...
                    [obj.PeakList.data(:,1); obj.PeakList.data(:,2)];
                for jj = 1:length(lim2Plot(:,1))
                    lim2Plot(jj,2) = ...
                        obj.DataIn(obj.DataIn(:,1) == lim2Plot(jj,1),2);
                    lim2Plot(jj,3) = ...
                       obj.Baseline.vals(obj.DataIn(:,1) == lim2Plot(jj,1));
                end
                
                plot(obj.DataIn(:,1), obj.DataIn(:,2));
                plot(obj.DataIn(:,1), obj.Baseline.vals(:), 'r')
                stem(lim2Plot(:,1), lim2Plot(:,2), 'k', 'Marker', 'none')
                stem(lim2Plot(:,1), lim2Plot(:,3), 'w', 'Marker', 'none')
                
                fprintf('peakNbr \t%s \t%s \t%s \t%s \t%s \t%s \t%s \t%s', obj.PeakList.headings{:});
                fprintf('\n');
                for jj = 1:length(obj.PeakList.data(:,1))
                    fprintf('\n');
                    fprintf('#:%i \t%.4g \t%.4g \t%.4g \t%.4g \t%.4g \t%.4g \t%.4g \t%.4g', [jj, obj.PeakList.data(jj,:)]);
                    
                end
                fprintf('\n');
                
                
                hold off

        end
        
    end
    
    
    
    
    
end

