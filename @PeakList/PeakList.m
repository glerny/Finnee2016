%% DESCRIPTION
% PeakList is the class that is contain all PIP extracted from a single
% centroid mode dataset
%
%% LIST OF THE CLASS'S PROPERTIES
%
%% LIST OF THE CLASS'S METHODS
%
%% Copyright
% BSD 3-Clause License
% Copyright 2016-2017 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef PeakList
    
    properties
        BPP         % Base peak profile calculated using all the PIPs
        TIP         % Total ion profile calculated using all the PIPs
        LstPIP      % List of all the PIP
        FOM         % summary of all FOMs
        % #PIP|max I|time @ max I|M0|M1|M2|M3|mean(m/z)|std(m/z)|AccMass
        
        sFOM        = nan(1, 11);
        nFOM        = ones(1, 11);
        Type        % 'singleton'  dtsIn is a 'centroid' dataset
        % 'replicates' dtsIn is a array of peaklist
        Path2Fin
        Log2crea
        p4norm
        options     % Options
        AxisX
        AxisY
        AxisZ
        Path2Pkl    % where to save
        Noise
        Ntm         = 0;
        Nmz         = 0;      
    end
    
    methods
        function obj = PeakList(dtsIn, prm1, prm2, prm3, prm4)
            % creator method will group as PIP any series of points that
            % does not differ in their m/z by more than prm2. PIP will be
            % recorded only if it contain at least prm3 whith at least
            % one points of intensity higher the prm1.
            if nargin == 0
                obj.BPP      = {};
                obj.TIP      = {};
                obj.LstPIP   = {};
                obj.FOM      = {};
                obj.Type     = 'blank';
                obj.Path2Fin = {};
                obj.Log2crea = {};
                obj.options  = {};
                obj.AxisX    = {};
                obj.AxisY    = {};
                obj.AxisZ    = {};
                obj.Path2Pkl = '';
                obj.p4norm   = {};
                obj.Noise    = {};
                
            elseif isa(dtsIn, 'Dataset')
                
                if ~strcmp(dtsIn.Format, 'centroid')
                    error('dtsIn should be centroid mode')
                end
                if nargin < 5
                    prm4 = [0 inf];
                end
                obj.Type            = 'singleton';
                obj.options.InfoDts = dtsIn.InfoDts;
                obj.options.ThIt    = prm1;
                obj.options.ThMZ    = prm2;
                obj.options.ThBk    = 1;
                obj.options.minPts  = prm3;
                obj.options.XLim    = prm4;
                obj.Path2Fin{1}     = dtsIn.Path2Fin;
                obj.Log2crea{1}     = dtsIn.Log;
                obj.p4norm          = {};
                obj.Noise           = dtsIn.Noise;
                
                InfoAxis     = dtsIn.AxisX.InfoAxis;
                InfoAxis.Loc = 'inAxis';
                obj.AxisX{1} = Axis(InfoAxis, dtsIn.AxisX.Data);
                obj.AxisY{1} = Axis(dtsIn.AxisY.InfoAxis);
                obj.AxisZ{1} = Axis(dtsIn.AxisZ.InfoAxis);
                [obj.Path2Pkl, ~, ~] = fileparts(dtsIn.InfoDts.Path2Dat{1});
                
                LoPts = [];
                X      = obj.AxisX{1}.Data;
                IdS    = find(X >= prm4(1), 1, 'first');
                IdE    = find(X <= prm4(2), 1, 'last');
                
                for ii = IdS:IdE
                    MS      = dtsIn.ListOfScans{ii}.Data;
                    MS(:,3) = ii;
                    LoPts   = [LoPts ; MS];
                end
                LoPts(LoPts(:,1) == 0, :) = [];
                obj.LstPIP{1}    = getPIP(LoPts, prm2, prm1, prm3, X, obj.options.InfoDts);
                obj.FOM{1}.Headings = {'Id', 'IntMax', 'Tm@IM', 'M0', 'M1',...
                    'M2', 'M3', 'mean(m/z)', 'std(m/z)', 'Acc. Mass', 'DetInt', 'S/N'};
                
                BPP      = obj.AxisX{1}.Data;
                BPP(:,2) = 0;
                TIP      = obj.AxisX{1}.Data;
                TIP(:,2) = 0;
                
                for ii = 1:length(obj.LstPIP{1} )
                    obj.FOM{1}.Data(ii, 1:11) = [ii, obj.LstPIP{1}{ii}.FOM];
                    IdN = findCloser(obj.LstPIP{1}{ii}.FOM(7),obj.Noise.Data(:,1));
                    obj.FOM{1}.Data(ii, 12) = obj.FOM{1}.Data(ii, 11)...
                        /max(obj.Noise.Data(IdN, 2), 15);
                    IdS = obj.LstPIP{1}{ii}.IdS;
                    y   = obj.LstPIP{1}{ii}.y;
                    IdE = IdS + size(y, 1) - 1;
                    
                     % Calculating BPP & TIP
                     TIP(IdS:IdE, 2) = TIP(IdS:IdE, 2) + y;
                     BPP(IdS:IdE, 2) = max(BPP(IdS:IdE, 2), y);  
                end
               
                InfoTrc      = dtsIn.BPP.InfoTrc;
                InfoTrc.Loc  = 'inTrace';
                obj.BPP{1}   = Trace(InfoTrc, BPP);
                InfoTrc      = dtsIn.TIP.InfoTrc;
                InfoTrc.Loc  = 'inTrace';
                obj.TIP{1}   = Trace(InfoTrc, TIP);
                
                myPeakList = obj; %#ok<*NASGU>
                save(fullfile(obj.Path2Pkl, 'myPeakList.mat'), 'myPeakList')
                
            else
                error('dtsIn not recognised')
            end
            
        end
        
        function save(obj)
            %% DESCRIPTION
            myPeakList = obj; %#ok<*NASGU>
            save(fullfile(obj.Path2Pkl, 'myPeakList.mat'), 'myPeakList')
        end
        
    end
end

