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
    end
    
    properties (Hidden)
        options     % Options
        AxisX
        AxisY
        AxisZ
        Path2PkL    % where to save
    end
    
    methods
        function obj = PeakList(dtsIn, ThIt, ThMZ, minPts, XLim) 
            % creator method will group as PIP any series of points that
            % does not differ in their m/z by more than ThMZ. PIP will be
            % recorded only if it contain at least minPts whith at least
            % one points of intensity higher the ThIt.
            if ~strcmp(dtsIn.Format, 'centroid')
                error('dtsIn should be centroid mode')
            end
            if nargin < 5
                XLim = [0 inf];
            end
            
            obj.options.InfoDts = dtsIn.InfoDts;
            obj.options.ThIt    = ThIt;
            obj.options.ThMZ    = ThMZ;
            obj.options.ThBk    = 1;
            obj.options.minPts  = minPts;
            obj.options.XLim    = XLim;
            
            InfoAxis     = dtsIn.AxisX.InfoAxis;
            InfoAxis.Loc = 'inAxis';
            obj.AxisX    = Axis(InfoAxis, dtsIn.AxisX.Data);
            obj.AxisY    = Axis(dtsIn.AxisY.InfoAxis);
            obj.AxisZ    = Axis(dtsIn.AxisZ.InfoAxis);
            InfoTrc      = dtsIn.BPP.InfoTrc;
            InfoTrc.Loc  = 'inTrace';
            obj.BPP      = Trace(InfoTrc, dtsIn.BPP.Data);
            InfoTrc      = dtsIn.TIP.InfoTrc;
            InfoTrc.Loc  = 'inTrace';
            obj.TIP      = Trace(InfoTrc, dtsIn.TIP.Data);
            [obj.Path2PkL, ~, ~] = fileparts(dtsIn.InfoDts.Path2Dat{1});

            LoPts = [];
            X      = obj.AxisX.Data;
            IdS    = find(X >= XLim(1), 1, 'first');
            IdE    = find(X <= XLim(2), 1, 'last');
            X      = X(IdS:IdE);
            
            for ii = IdS:IdE
                MS      = dtsIn.ListOfScans{ii}.Data;
                MS(:,3) = ii;
                LoPts   = [LoPts ; MS];
            end
            LoPts(LoPts(:,1) == 0, :) = [];
            obj.LstPIP       = getPIP(LoPts, ThMZ, ThIt, minPts, X, obj.options.InfoDts);
            obj.FOM.Headings = {'Id', 'IntMax', 'Tm@IM', 'M0', 'M1',...
                'M2', 'M3', 'mean(m/z)', 'std(m/z)', 'Acc. Mass'};
            for ii = 1:length(obj.LstPIP)
                obj.FOM.Data(ii, :) = [ii, obj.LstPIP{ii}.FOM];
            end
            myPeakList = obj; %#ok<*NASGU>
            save(fullfile(obj.Path2PkL, 'myPeakList.mat'), 'myPeakList')
            
        end
        
         function save(obj)
            %% DESCRIPTION
            myPeakList = obj; %#ok<*NASGU>
            save(fullfile(obj.Path2PkL, 'myPeakList.mat'), 'myPeakList')
 
        end
        
       
    end
end

