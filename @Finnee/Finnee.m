classdef Finnee
    
    properties (Hidden = true)
        MZMLDump
        Options  
        Path2Fin
    end
    
    properties
        FileID
        Datasets
    end
    
    methods
        function obj = Finnee(varargin)
            
            % Check the options and create the Finnee object
            obj.Options = Options4Finnee(varargin);
            
            switch obj.Options.FileFormat
                case 'mzML'
                    obj = domzML2Finnee(obj);
            end
            
            myFinnee = obj;
            save(fullfile(obj.Path2Fin, 'myFinnee.mat'), 'myFinnee')
            
        end
        
        function queryMZMLDump(obj)
        end
        
        function doSave(obj)
            myFinnee = obj;
            save(fullfile(obj.Path2Fin, 'myFinnee.mat'), 'myFinnee')
        end
        
        function toMZML(obj)
        end
        
        function obj = doBaseline(obj)
        end
        
        function obj = doCentroid(obj)
            m = length(obj.Datasets) + 1;
            newDts = doCtr(obj.Datasets{dts}, par4ctr, ...
                obj.Path2Fin, obj.Options.MaxFileSize, m);
            newDts.CreatedFrom = obj.Datasets{dts}.Log;
            obj.Datasets{end+1} = newDts;
            
            myFinnee = obj;
            save(fullfile(obj.Path2Fin, 'myFinnee.mat'), 'myFinnee')
        end
        
        function obj = doFilter(obj)
        end
        
    end
    
end

