classdef Axe 
    
    properties 
        XLabel   	= 'Time'
        XUnit       = ''
        Path2Dat    = ''
        Index       = [0 0]
        Precision   = 'single'
        Path2Fin    = ''
        Log         = ''
    end
    
    properties (Dependent)
        Data
    end
    
    methods
        function obj = Axe(structInfo, data2write)
            if nargin == 0
            else
                obj.XLabel      = structInfo.XLabel;
                obj.XUnit       = structInfo.XUnit;
                obj.Path2Dat    = structInfo.Path2Dat;
                obj.Path2Fin    = structInfo.Path2Fin;
                obj.Log         = structInfo.Log;
                
                % DATA ARE WRITTEN IN SINGLE PRECISION FOR SPEED OPTIMISATION
                fidWriteDat     = fopen(obj.Path2Dat, 'ab');
                obj.Index(1)    =  ftell(fidWriteDat);
                fwrite(fidWriteDat, data2write, obj.Precision);
                obj.Index(2)    =  ftell(fidWriteDat);
                fclose(fidWriteDat);
            end
        end
        
        function data = get.Data(obj)
            fidReadDat  = fopen(obj.Path2Dat, 'rb');
            fseek(fidReadDat,  obj.Index(1), 'bof');
            data        = fread(fidReadDat, ...
                [( obj.Index(2)- obj.Index(1))/4 1], 'single');
            fclose(fidReadDat);
        end
    end
    
end

