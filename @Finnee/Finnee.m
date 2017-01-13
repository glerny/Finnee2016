classdef Finnee

    properties
        MZMLDump
        Path2Fin  
        Options
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
        
        function queryMZML(obj, what)
        end
        
        function save(obj)
            myFinnee = obj;
            save(fullfile(obj.Path2Fin, 'myFinnee.mat'), 'myFinnee')
        end
        
        function toMZML(obj, dts)
        end
        
        function obj = doBaseline(obj, dts, par4bas)
            m = length(obj.Datasets) + 1;
            newDts = doBslCor(obj.Datasets{dts}, par4bas, ...
                obj.Path2Fin, obj.Options.MaxFileSize, m);
            newDts.CreatedFrom = obj.Datasets{dts}.Log;
            newDts.Log = ['MSPROFDTS DTS=', num2str(m), ...
                ' BLC=1 ', ...
                par4bas.type, '=', num2str(par4bas.type),...
                ' NOI=', num2str(par4bas.noise)];
            
            obj.Datasets{end+1} = newDts;
            
            myFinnee = obj;
            save(fullfile(obj.Path2Fin, 'myFinnee.mat'), 'myFinnee')
        end
        
        function obj = doCentroid(obj, dts, par4ctr)
            m = length(obj.Datasets) + 1;
            newDts = doCtr(obj.Datasets{dts}, par4ctr, ...
                obj.Path2Fin, obj.Options.MaxFileSize, m);
            newDts.CreatedFrom = obj.Datasets{dts}.Log;
            obj.Datasets{end+1} = newDts;
            
            myFinnee = obj;
            save(fullfile(obj.Path2Fin, 'myFinnee.mat'), 'myFinnee')
        end
        
        function obj = remNoise(obj, dts, varargin)
            % the implemented filtes are:
            % - 'spikes' - followed by 1 or 2 (length of spikes as pts)
            % - 'filterMS'- followed by a boolean array of the same length
            % as the MZ axes with true for mz calues to remove and false
            % for the one to keep (i.e. filter =
            % myFinnee.Datasets{2}.BIS.Data(:,2) <= 1000)
                
            switch lower(varargin{1})
                case 'spikes'
                    if ~strcmp(obj.Datasets{dts}.Format, 'MS profile')
                        error('Error.\nDataset %i shuld be MS profile as format.\n not %s',...
                            dts, obj.Datasets{dts}.Format)
                    end
                    m = length(obj.Datasets) + 1;
                    newDts = remSpikes(obj.Datasets{dts}, varargin{2}, obj.Path2Dat, m);
                    newDts.CreatedFrom = obj.Datasets{dts}.Log;
                    newDts.Log = ['MSPROFDTS DTS=', num2str(m), ...
                        ' SPK=', num2str( varargin{2})];
                    
                case 'filtermz'
                    if ~strcmp(obj.Datasets{dts}.Format, 'MS profile')
                        error('Error.\nDataset %i shuld be MS profile as format.\n not %s',...
                            dts, obj.Datasets{dts}.Format)
                    end
                    m = length(obj.Datasets) + 1;
                    newDts = filterByMZ(obj.Datasets{dts}, varargin{2}, ...
                        obj.Path2Fin, obj.Options.MaxFileSize, m);
                    newDts.CreatedFrom = obj.Datasets{dts}.Log;
                    newDts.Log = ['MSPROFDTS DTS=', num2str(m), ...
                        ' FMZ=1'];
                    
                otherwise
                    error('Error.\nThe method %s does not exist in this context',...
                        varargin{1})
            end
            
            obj.Datasets{end+1} = newDts;
                        
            myFinnee = obj;
            save(fullfile(obj.Path2Fin, 'myFinnee.mat'), 'myFinnee')
        end
        
    end
    
end