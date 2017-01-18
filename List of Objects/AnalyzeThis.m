classdef AnalyzeThis
    
    properties
        %% Options for fileInfo
        DataIn              = []
        % First colum: axe; following traces (can be more than 1 but should
        % be related
        
        BaselineMethod  	= 'ArPLS:10E6:0.001'
        % Possible methods:
        %   'None'      	No additional parameter
        % 	'Polyn:p1'      Polynomial fitting, p1 is the degree of the
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
        %   'LmMm:p1'       Peak detected as a local Maxime between 2
        %                   local minimaux, p1 is the number of neighbourgh 
        %                   used to detect local maximum or minimum (1:10),
        %                   p2 is a intensity threshold below witch peaks
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
            obj.DataIn = dataIn;
            
            if nargin  < 2
                return
            end
            
            input = @(x) find(strcmpi(varargin,x),1);
            
            tgtIx = input('baseline');
            if ~isempty(tgtIx)
                str2dec = varargin{tgtIx +1};
                C =  regexp(str2dec, ':', 'split');
                switch lower(C{1})
                    case 'none'
                        obj.BaselineMethod = struct('name', 'None');
                        
                    case 'constant'
                        obj.BaselineMethod = struct('name', 'Constant');
                        
                    case 'linear'
                        obj.BaselineMethod = struct('name', 'Linear');
                        
                    case 'arpls'
                        if length(C) < 2
                            C(2) = '10E6';
                        end
                        
                        obj.BaselineMethod = struct(...
                            'name', 'ArPLS',...
                            'param1', str2double(C(2)) );
                end
            end
            
            tgtIx = input('peakpicking');
            if ~isempty(tgtIx)
                str2dec = varargin{tgtIx +1};
                C =  regexp(str2dec, ':', 'split');
                switch lower(C{1})
                    case 'lmmm'
                        if length(C) < 2
                            C(2) = '2';
                        end
                        
                        obj.PeakPickingMethod = struct(...
                            'name', 'LmMm',...
                            'param1', str2double(C(2)) );
                end
            end
            
            tgtIx = input('threshold');
            if ~isempty(tgtIx)
                obj.Threshold = varargin{tgtIx +1};
            end
        end
        
        function baseline = get.Baseline(obj)
            
            baseline.definition = obj.BaselineMethod;
            baseline.noise(1) = NaN;
            
            switch obj.BaselineMethod.name
                case 'ArPLS'
                    [~, n] = size(obj.DataIn);
                    results(:,1) = obj.DataIn(:,1);
                    lambda = obj.BaselineMethod.param1;
                    
                    for ii = 2:n
                        yori = obj.DataIn(:,ii);
                        indNotZeros = yori ~= 0;
                        [z, bslPts, ~] = doArPLS(yori(indNotZeros), lambda);
                        results(indNotZeros, ii) = z;
                        baseline.noise(ii) = 4*std(z(bslPts));
                    end
                    
                case 'None'
                    [~, n] = size(obj.DataIn);
                    results = zeros(size(obj.DataIn));
                    results(:,1) = obj.DataIn(:,1);
                    
                    for ii = 2:n
                        yori = obj.DataIn(:,ii);
                        baseline.noise(ii) = 4*std(yori);
                    end
            end
            baseline.values = results;
        end
        
        function peakList = get.PeakList(obj)
            [~, n] = size(obj.DataIn);
            peakList(n - 1) = struct();
            for ii = 2:n;
                if isempty(obj.Threshold)
                    thrs = 3* obj.Baseline.noise(ii);
                else
                    thrs = obj.Threshold;
                end
                
                yOri = obj.DataIn(:,ii) - obj.Baseline.values(:,ii);
                
                switch obj.PeakPickingMethod.name
                    case 'LmMm'
                        X = spread(yOri, obj.PeakPickingMethod.param1);
                        
                        locMax = find(X(:,((1+end)/2)) == max(X, [], 2) & max(X, [], 2) >= thrs);
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
            end
            peakList.data = FOM;      
        end
        
        function plotAnalysis(obj)
            [~, n] = size(obj.DataIn);
            for ii = 2:n
                figure('Name', 'Result AnalyzeThis');
                hold on
                lim2Plot(:,1) = ...
                    [obj.PeakList.data(:,1); obj.PeakList.data(:,2)];
                for jj = 1:length(lim2Plot(:,1))
                    lim2Plot(jj,2) = ...
                        obj.DataIn(obj.DataIn(:,1) == lim2Plot(jj,1),2);
                    lim2Plot(jj,3) = ...
                       obj.Baseline.values(obj.DataIn(:,1) == lim2Plot(jj,1),2);
                end
                
                plot(obj.DataIn(:,1), obj.DataIn(:,2));
                plot(obj.Baseline.values(:,1), obj.Baseline.values(:,2), 'r')
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
    
    
    
    
    
end

