classdef PIP
    %PIP Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Data
        x
        IdS
        AxisX
        AxisY
        AxisZ
    end
    
    properties (Hidden)
        deconvolved = false
    end
    
    properties (Dependent)
        y
        FOM
    end
    
    methods
        function obj  = PIP(DataIn, axe, infoDts)
            obj.Data  = DataIn;
            obj.x     = axe(min(DataIn(:,3)):max(DataIn(:,3)));
            obj.IdS   = min(DataIn(:,3));
            obj.AxisX = infoDts.AxisX;
            obj.AxisY = infoDts.AxisY;
            obj.AxisZ = infoDts.AxisZ;
        end
        
        function y   = get.y(obj)
            y = zeros(size(obj.x,1), 1);
            C = unique(obj.Data(:,3));
            
            if length(C) ~= size(obj.Data, 1)
                for jj = 1:length(C)
                    y(C(jj) - min(C) + 1) = ...
                        sum(obj.Data(obj.Data(:,3) == C(jj),2));
                end
            else
                y(obj.Data(:,3) - min(obj.Data(:,3))+1) = obj.Data(:,2);
            end
        end
        
        function FOM = get.FOM(obj)
            [FOM(1), Id2x] = max(obj.y);
            FOM(2) = obj.x(Id2x);
            FOM(3) = trapz(obj.x, obj.y);
            FOM(4) = trapz(obj.x, obj.x.*obj.y)/FOM(3);
            FOM(5) = trapz(obj.x, (obj.x-FOM(4)).^2.*obj.y)/FOM(3);
            FOM(6) = trapz(obj.x, (obj.x-FOM(4)).^3.*obj.y)/FOM(3);
            FOM(7) = mean(obj.Data(:,1));
            FOM(8) = std(obj.Data(:,1));
            FOM(9) = sum(obj.Data(:,1).*obj.Data(:,2))/sum(obj.Data(:,2));
        end
        
        
        function plot(obj, stringTitle)
            if nargin == 1, stringTitle = ''; end
            InputFig = figure(                          ...
                'Visible'          , 'on'             , ...
                'Name'             , stringTitle      , ...
                'Toolbar'          , 'none'           , ...
                'MenuBar'          , 'none'           , ...
                'Units'            , 'normalized'     , ...
                'WindowStyle'      , 'normal'         , ...
                'Resize'           , 'on');
            
            AxisH1 = axes(InputFig ,...
                'Units'            , 'normalized'     , ...
                'OuterPosition'    , [0 0 2/3 1]);
            
            hold on
            yyaxis left
            hP = plot(AxisH1, obj.x, obj.y);
            yyaxis right
            sP = scatter(AxisH1, ...
                obj.x(obj.Data(:,3)- min(obj.Data(:,3)) + 1), obj.Data(:,1));
            
            title('Intensities and accurate masses')
            xlabel([obj.AxisX.Label, ' / ', obj.AxisX.Unit]);
            AxisH1.YAxis(1).Label.String = ...
                [obj.AxisZ.Label, ' / ', obj.AxisZ.Unit];
            AxisH1.YAxis(2).Label.String = ...
                [obj.AxisY.Label, ' / ', obj.AxisY.Unit];

            hold off
            
            String4Edit{1} = '  FIGURES OF MERITS';
            String4Edit{2} = '_____________________';
            String4Edit{3} = '';
            String4Edit{4} = sprintf(['Max Intensity : '...
                , obj.AxisZ.fo, ' %s'], obj.FOM(1),  obj.AxisZ.Unit);
            String4Edit{5} = sprintf(['Time @ IMax  : '...
                , obj.AxisX.fo, ' %s'], obj.FOM(2),  obj.AxisX.Unit);
            String4Edit{6} = sprintf(['Peak Area      : '...
                , '%.3e', ' %s'], obj.FOM(3),  'a.u.');
            String4Edit{7} = sprintf(['Peak center   : '...
                , obj.AxisX.fo, ' %s'], obj.FOM(4),  obj.AxisX.Unit);
            String4Edit{8} = sprintf(['Peak variance: '...
                , '%.3e', ' %s^2'], obj.FOM(5), obj.AxisX.Unit);
            String4Edit{10} = sprintf(['m/z = ', obj.AxisY.fo, ...
                ' +/- '  obj.AxisY.fo], obj.FOM(7), obj.FOM(8));
            String4Edit{11} = sprintf(['Accurate mass = ', obj.AxisY.fo],...
                obj.FOM(9));
            
            
            EditH1 = uicontrol(InputFig ,...
                'Style'              , 'Edit'      , ...
                'Visible'            , 'on'        , ...
                'Units'              , 'normalized', ...
                'Max'                , 2           , ...
                'String'             , String4Edit ,...
                'HorizontalAlignment', 'left', ...
                'Position'           , [2/3 0 1/3 1]);
            
            PushH1 = uicontrol(InputFig ,...
                'Style'   , 'pushbutton', ...
                'Visible' , 'on'        , ...
                'Units'   , 'normalized', ...
                'Max'     , 2           , ...
                'String'  , 'Deconvolve',...
                'Tag'     , 'PushH1    ',...
                'Callback', @pushMe     , ...
                'Position', [2/3 0 1/3 1]);
            wPB1            = PushH1.Extent(3);
            hPB1            = PushH1.Extent(4);
            PushH1.Position = [2/3+0.1*wPB1 0+0.1*hPB1  wPB1  hPB1];
            
            function pushMe(source, ~)
            end
            
        end
    end
    
end

