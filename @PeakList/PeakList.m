classdef PeakList
    
    properties
        options
        AxisX
        AxisY
        AxisZ
        BPP
        TIP
        LstPIP
        FOM
        Path2PkL
    end
    
    methods
        function obj = PeakList(dtsIn, ThIt, ThMZ, ThBk, minPts, XLim)
            if ~strcmp(dtsIn.Format, 'centroid')
                error('dtsIn should be centroid mode')
            end
            obj.options.InfoDts = dtsIn.InfoDts;
            obj.options.ThIt    = ThIt;
            obj.options.ThMZ    = ThMZ;
            obj.options.ThBk    = ThBk;
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
            obj.LstPIP       = getPIP(LoPts, ThMZ, ThBk, ThIt, minPts, X, obj.options.InfoDts);
            obj.FOM.Headings = {'Id', 'IntMax', 'Tm@IM', 'M0', 'M1',...
                'M2', 'M3', 'mean(m/z)', 'std(m/z)'};
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

