classdef Trace 
    
    properties (Hidden = true)
        TraceType   = 'continuous'
        XLabel      = 'Time'
        XUnit       = ''
        YLabel      = 'Intensity'
        YUnit       = ''
        Path2Dat    = ''
        Index       = [0 0 0]
        Precision   = 'single'
        Variables   = []
        Path2Fin    = ''
        Log         = ''
        PrecInX     = '%.2f'
        PrecInY     = '%.4f'
    end
    
    properties
        Title     = ''
    end
    
    properties (Dependent)
        Data
    end
    
    methods
        
        function obj = Trace(structInfo, data2write)
            if nargin == 0
            else
                obj.Title      	= structInfo.title;
                obj.TraceType 	= structInfo.traceType;
                obj.XLabel    	= structInfo.XLabel;
                obj.XUnit       = structInfo.XUnit;
                obj.YLabel     	= structInfo.YLabel;
                obj.YUnit       = structInfo.YUnit;
                obj.Path2Dat  	= structInfo.Path2Dat;
                obj.Variables  	= structInfo.Variables;
                obj.Precision  	= structInfo.Precision;
                obj.Path2Fin    = structInfo.Path2Fin;
                obj.Log         = structInfo.Log;
                
                % DATA ARE WRITTEN IN SINGLE PRECISION FOR SPEED OPTIMISATION
                fidWriteDat  = fopen(obj.Path2Dat, 'ab');
                if~isempty(data2write)
                    obj.Index(1) =  ftell(fidWriteDat);
                    fwrite(fidWriteDat, data2write, obj.Precision);
                    obj.Index(2) =  ftell(fidWriteDat);
                    obj.Index(3) = length(data2write(1,:));
                else
                    obj.Index(1) =  ftell(fidWriteDat);
                    obj.Index(2) = ftell(fidWriteDat);
                    obj.Index(3) = 0;
                end
                fclose(fidWriteDat);
            end
        end
        
        function data = get.Data(obj)
            fidReadDat = fopen(obj.Path2Dat, 'rb');
            fseek(fidReadDat,  obj.Index(1), 'bof');
            switch obj.Precision
                case 'single'
                    data = fread(fidReadDat, [( obj.Index(2)- obj.Index(1))/...
                        (4* obj.Index(3)) obj.Index(3)], 'single');
                case 'double'
                    data = fread(fidReadDat, [( obj.Index(2)- obj.Index(1))/...
                        (8* obj.Index(3))  obj.Index(3)], 'double');
            end
            fclose(fidReadDat);
        end
        
        function plot(obj)
            fig = figure('Name', obj.Log);
            dcm_obj = datacursormode(fig);
            set(dcm_obj,'UpdateFcn',@myupdatefcn);
            c = uicontextmenu;
            m1 = uimenu(c, 'Label','ExportObj','Callback',@doExport);
            fig.UIContextMenu = c;
            
            switch obj.TraceType
                case 'MS profile'
                    plot(obj.Data(:,1), obj.Data(:,2));
                    title(obj.Title);
                    xlabel([obj.XLabel, ' / ',obj.XUnit]);
                    ylabel([obj.YLabel, ' / ',  obj.YUnit]);
                    
                case 'MS centroid'
                    stem(obj.Data(:,1), obj.Data(:,2), 'Marker', 'none');
                    title(obj.Title);
                    xlabel([obj.XLabel, ' / ',obj.XUnit]);
                    ylabel([obj.YLabel, ' / ',  obj.YUnit]);
                    
                case 'ion profile'
                    plot(obj.Data(:,1), obj.Data(:,2));
                    title(obj.Title);
                    xlabel([obj.XLabel, ' / ',obj.XUnit]);
                    ylabel([obj.YLabel, ' / ',  obj.YUnit]);
                    
            end
            
            function doExport(~, ~)
                
                assignin('base', 'currentTrace', obj)
            end
            
            function txt = myupdatefcn(~,event_obj)
                pos = get(event_obj,'Position');
                
                switch obj.TraceType
                    case {'MS profile', 'MS centroid'}
                        txt = {['m/z: ',num2str(pos(1), obj.PrecInY)],...
                            ['Intensity: ',num2str(pos(2))]};
                        
                    case 'ion profile'
                        txt = {['Time: ',num2str(pos(1), obj.PrecInX)],...
                            ['Intensity: ',num2str(pos(2))]};
                end
            end
        end
    end
    
end

