%% DESCRIPTION
%
%% Copyright
% BSD 3-Clause License Copyright 2016-2017 G. Erny (guillaume@fe.up,pt),
% FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef AnalyzeThis
    
    properties
        TraceIn       % A trace object that contain all the information
        % about the profile to analyxe
        
        SmoothMethod % A sting with the method name followed by the
        % necessary parameters e.g. 'methodName:para1:para2:...:paran'
        % Possible methods:
        %   'None'          No additional parameter
        
        BaseMethod    % A sting with the method name followed by the
        % necessary parameters e.g. 'methodName:pm1:pm2:...:pmn' Possible
        % methods:
        %   'None'
        %       No additional parameter
        % 	'PolyFit:pm1'
        %       Polynomial fitting, pm1 is the degree of the polynomial
        %       (pm1=0: constant; pm1=1: linear...)
        %   'ArPLS:pm1:pm2'
        %       Baseline correction using asymmetrically reweighted
        %       penalized least squares smoothing, (doi:
        %       10.1039/c4an01061b), pm1: lambda(10E2-10E9); pm2: alpha
        %   'ArPLS2:pm1'
        %       Modified ArPLS method, pm1: lambda (10E2-10E9
        
        PeakPicking   % A sting with the method name followed by the
        % necessary parameters e.g. 'methodName:pm1:pm2:...:pamn
        % Possible methods:
        %   'localMax:pm1:pm2'
        %       Detect local Maxima, p1 is the number of neighbourgh used
        %       to detect the local maximum (1:10), p2 is an intensity
        %       threshold below witch peaks will not be processed
        
        Func4Deconv % A sting with the method name followed by the
        % necessary parameters e.g. 'methodName:para1:para2:...:paran'
        % Possible methods:
        %   'None'          No additional parameter
        
        Options             % Structure with all the options
        
    end
    
    properties (Dependent)
        XY
        SmoothData
        Baseline
        PeakList
        Fitting
    end
    
    methods
        function obj = AnalyzeThis(dataIn, varargin)
            narginchk(1, inf)
            
            if ~isa(dataIn, 'Trace')
                error('the input parameter for AnalyzeThis should be a non-empty Trace');
            end
            if isempty(dataIn.Data)
                error('the input parameter for AnalyzeThis should be a non-empty Trace');
            end
            
            % Set default parameters
            obj.Func4Deconv      = 'None';
            obj.Options.NegPts   = true;
            obj.Options.rounding = false;
            obj.Options.XLim     = [0 inf];
            InfoTrc              = dataIn.InfoTrc;
            InfoTrc.Loc          = 'inTrace';
            obj.TraceIn          = Trace(InfoTrc, dataIn.Data);
            obj.BaseMethod       = 'ArPLS2:10E6';
            obj.SmoothMethod     = 'None';
            obj.PeakPicking      = ['localMax:2:', ...
                num2str(3*obj.Baseline.noise, '%.0f')];
            
            %%% Decipher varargin
            input = @(x) find(strcmpi(varargin,x),1);
            
            tgtIx = input('SmoothMethod');
            if ~isempty(tgtIx)
                obj.SmoothMethod = varargin{tgtIx +1};
            end
            
            tgtIx = input('BaseMethod');
            if ~isempty(tgtIx)
                obj.BaseMethod   = varargin{tgtIx +1};
            end
            
            tgtIx = input('PeakPicking');
            if ~isempty(tgtIx)
                obj.PeakPicking  = varargin{tgtIx +1};
            end
            
            tgtIx = input('Func4Deconv');
            if ~isempty(tgtIx)
                obj.Func4Deconv  = varargin{tgtIx +1};
            end
            
            
        end
        
        function xy = get.XY(obj)
            xy      = obj.SmoothData;
            IdS     = find(xy(:,1) <=  obj.Options.XLim(1), 1, 'last');
            if isempty(IdS), IdS = 1; end
            IdE     = find(xy(:,1) >=  obj.Options.XLim(2), 1, 'first');
            if isempty(IdE), IdE = size(xy,1); end
            xy      = xy(IdS:IdE, :);
            xy(:,2) = xy(:,2) - obj.Baseline.vals;
            
            if obj.Options.NegPts == false
                xy(xy(:,2) < 0, 2) = 0;
            end
            
            if obj.Options.rounding
                xy(:,2) = round(xy(:,2));
            end
            
            
            
        end
        
        function obj = set.SmoothMethod(obj, value)
            basDef =  strsplit(value, ':');
            
            switch lower(basDef{1})
                case 'none'
                    obj.SmoothMethod = 'None';
                    
                otherwise
                    error('%s is not a recognised smoothing method', basDef{1})
            end
        end
        
        function smoothData = get.SmoothData(obj)
            basDef =  strsplit(obj.SmoothMethod, ':');
            
            switch basDef{1}
                case 'None'
                    XY         = obj.TraceIn.Data;
                    smoothData = XY;
            end
        end
        
        function obj = set.BaseMethod(obj, value)
            basDef =  strsplit(value, ':');
            
            switch lower(basDef{1})
                case 'none'
                    obj.BaseMethod = 'None';
                    
                case 'polyfit'
                    if length(basDef) == 1;
                        basDef{2} = '2';
                    end
                    if int8(str2double(basDef{2})) < 0
                        warning('The parameter for the polynomial fitting should be positive or null, the value has been set to 0');
                        basDef{2} = '0';
                    elseif int8(str2double(basDef{2})) > 5
                        warning('The parameter for the polynomial fitting should be lower than 5, the value has been set to 5');
                        basDef{2} = '5';
                    else
                        basDef{2} = num2str(int8(str2double(basDef{2})));
                    end
                    obj.BaseMethod = ['PolyFit:', basDef{2}];
                    
                case 'arpls'
                    if length(basDef) == 1;
                        basDef{2} = '10E6';
                    end
                    if length(basDef) == 2;
                        basDef{3} = '0.1';
                    end
                    
                    if double(str2double(basDef{2})) < 10E2
                        warning('Lambda should be greater than 10E2, the value has been set to 10E2');
                        basDef{2} = '10E2';
                    elseif double(str2double(basDef{2})) > 10E9
                        warning('Lambda should be be lower than 10E9, the value has been set to 10E9');
                        basDef{2} = '10E9';
                    else
                        basDef{2} = num2str(double(str2double(basDef{2})));
                    end
                    
                    if double(str2double(basDef{3})) < 0
                        warning('Alpha should be positive, the value has been set to 0');
                        basDef{3} = '0';
                    elseif double(str2double(basDef{3})) > 1
                        warning('Alpha should be lower than 1, the value has been set to 1');
                        basDef{2} = '1';
                    else
                        basDef{3} = num2str(double(str2double(basDef{3})));
                    end
                    
                    obj.BaseMethod = ['ArPLS:', basDef{2}, ':', basDef{3}];
                    
                case 'arpls2'
                    if length(basDef) == 1;
                        basDef{2} = '10E6';
                    end
                    
                    if double(str2double(basDef{2})) < 10E2
                        warning('Lambda should be greater than 10E2, the value has been set to 10E2');
                        basDef{2} = '10E2';
                    elseif double(str2double(basDef{2})) > 10E9
                        warning('Lambda should be be lower than 10E9, the value has been set to 10E9');
                        basDef{2} = '10E9';
                    else
                        basDef{2} = num2str(double(str2double(basDef{2})));
                    end
                    
                    obj.BaseMethod = ['ArPLS2:', basDef{2}];
                    
                otherwise
                    error('%s is not a recognised baseline correction method', basDef{1})
            end
        end
        
        function baseline = get.Baseline(obj)
            basDef =  strsplit(obj.BaseMethod, ':');
            XY = obj.SmoothData;
            IdS     = find(XY(:,1) <=  obj.Options.XLim(1), 1, 'last');
            if isempty(IdS), IdS = 1; end
            IdE     = find(XY(:,1) >=  obj.Options.XLim(2), 1, 'first');
            if isempty(IdE), IdE = size(XY,1); end
            XY      = XY(IdS:IdE, :);
            
            switch basDef{1}
                case 'None'
                    baseline.bckgPts = false(XY,1);
                    baseline.noise = 4*std(XY(:,2));
                    baseline.vals = zeros(m, 1);
                    
                case 'PolyFit'
                    [z, bslPts] = doPF(XY, str2double(basDef{2}));
                    baseline.bckgPts = bslPts;
                    baseline.noise = 4*std(XY(bslPts,2) - z(bslPts));
                    baseline.vals = z;
                    
                case 'ArPLS'
                    lambda = str2double(basDef{2});
                    ratio = str2double(basDef{3});
                    [m, ~] = size(XY);
                    iNZ = XY(:,2) ~= 0;
                    [z, bslPts] = doArPLS(XY(iNZ, 2), lambda, ratio);
                    baseline.bckgPts = false(m,1);
                    baseline.bckgPts(iNZ) = bslPts;
                    baseline.vals = zeros(m, 1);
                    baseline.vals(iNZ) = z;
                    baseline.noise = 4*std(nonzeros(XY(baseline.bckgPts,2) ...
                        - baseline.vals(baseline.bckgPts)));
                    
                case 'ArPLS2'
                    lambda = str2double(basDef{2});
                    [m, ~] = size(XY);
                    iNZ = XY(:,2) ~= 0;
                    [z, bslPts] = doArPLS2(XY(iNZ, 2), lambda);
                    baseline.bckgPts = false(m,1);
                    baseline.bckgPts(iNZ) = bslPts;
                    baseline.vals = zeros(m, 1);
                    baseline.vals(iNZ) = z;
                    baseline.noise = 4*std(nonzeros(XY(baseline.bckgPts,2) ...
                        - baseline.vals(baseline.bckgPts)));
            end
            
%             %             if obj.Options.plot == true
%                             s             = subplot(2,1, 1);
%                             s.Parent.Name = ['Baseline Method: ', obj.BaseMethod];
%                             plot(XY(:,1), XY(:,2), 'k');
%                             title(['Original profile: Noise estimated = ', ...
%                                 num2str(baseline.noise, '%.0f'), ' ', obj.TraceIn.AxeY.Unit]);
%                             xlabel([obj.TraceIn.AxeX.Label, ' / ', obj.TraceIn.AxeX.Unit]);
%                             ylabel([obj.TraceIn.AxeY.Label, ' / ', obj.TraceIn.AxeY.Unit]);
%                             hold on
%                             plot(XY(baseline.bckgPts,1), ...
%                                 XY(baseline.bckgPts,2), 'or');
%                             plot(XY(:,1), baseline.vals, 'r')
%                             hold off
%             
%                             subplot(2,1, 2);
%                             y_ =  round(XY(:,2) -  baseline.vals);
%                             y_ (y_ < 0 ) = 0;
%                             plot(XY(:,1), y_, 'k');
%                             title('Corrected profile');
%                             xlabel([obj.TraceIn.AxeX.Label, ' / ', obj.TraceIn.AxeX.Unit]);
%                             ylabel([obj.TraceIn.AxeY.Label, ' / ', obj.TraceIn.AxeY.Unit]);
%             %
%             %             end
        end
        
        function obj = set.PeakPicking(obj, value)
            basDef =  strsplit(value, ':');
            
            switch lower(basDef{1})
                case 'none'
                    obj.PeakPicking = 'None';
                    
                case 'localmax'
                    if length(basDef) == 1;
                        basDef{2} = '2';
                    end
                    
                    if int8(str2double(basDef{2})) < 1
                        warning('The first paramter for LocalMax should be greater than 1');
                        basDef{2} = '1';
                    elseif int8(str2double(basDef{2})) > 10
                        warning('The first paramter for LocalMax should be lowe than 10');
                        basDef{2} = '10';
                    else
                        basDef{2} = num2str(int8(str2double(basDef{2})));
                    end
                    
                    if double(str2double(basDef{3})) < 0
                        warning('The intensity threshold should be positive, it has has been set to 0');
                        basDef{3} = '0';
                    else
                        basDef{3} = num2str(double(str2double(basDef{3})));
                    end
                    
                    obj.PeakPicking = ['LocalMax:', basDef{2}, ':', basDef{3}];
                    
                otherwise
                    error('%s is not a recognised peak picking method', basDef{1})
            end
        end
        
        function peakList = get.PeakList(obj)
            basDef =  regexp(obj.PeakPicking, ':', 'split');
            XY     = obj.XY;
            
            switch basDef{1}
                case 'LocalMax'
                    % Get values
                    XYpl    = LocalMaxima(XY, str2double(basDef{2}), ...
                        str2double(basDef{3}));
                    
                    peakList.headings = {'X@Apex', 'Int@Apex', ...
                        'peakWidth@05'};
                    peakList.H4CTR    = {'AM1', 'IT1', 'OT1'};
                    peakList.Data     = XYpl;
                    
%                     %if obj.Options.plot == true
%                         figure('Name', ...
%                             ['Peak picking Method: ', obj.PeakPicking]);
%                         plot(XY(:,1), XY(:,2), 'k');
%                         xlabel([obj.TraceIn.AxeX.Label, ...
%                             ' / ', obj.TraceIn.AxeX.Unit]);
%                         ylabel([obj.TraceIn.AxeY.Label, ...
%                             ' / ', obj.TraceIn.AxeY.Unit]);
%                         hold on
%                         stem(XYpl(:,1), XYpl(:,2), 'r');
%                         hold off
%                    % end
            end
        end
        
        function obj = set.Func4Deconv(obj, value)
            basDef =  strsplit(value, ':');
            
            switch lower(basDef{1})
                case 'none'
                    obj.Func4Deconv = 'None';
                    
                case 'gauss'
                    obj.Func4Deconv = 'Gauss';
                    
                case 'emg'
                    obj.Func4Deconv = 'EMG';
                    
                case 'pmg2'
                    obj.Func4Deconv = 'PMG2';
                    
                otherwise
                    error('%s is not a recognised peak function',...
                        basDef{1});
            end
        end
        
        function fitting = get.Fitting(obj)
            basDef =  strsplit(obj.Func4Deconv, ':');
            if strcmp(basDef{1}, 'None')
                fitting = {};
                return
            end
            
            prmIn = obj.PeakList.Data;
            if isempty(prmIn)
                fitting = {};
                return
            end
            
            basDef     = regexp(obj.Func4Deconv, ':', 'split');
            XY         = obj.XY;
            prmIn(:,3) = (prmIn(:,3)/2.355);
            
            [fitting.FitPrm, fitting.Int, fitting.model, fitting.SSE] = ...
                doPeakMinimization(XY(:,1), XY(:,2), basDef{1}, ...
                prmIn(:, [1,3]));
        end
    end
end

