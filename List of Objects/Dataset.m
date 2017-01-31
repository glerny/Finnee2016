%% DESCRIPTION
%
%% Copyright 
% BSD 3-Clause License
% Copyright 2016-2017 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


classdef Dataset
    
    properties
        Title           = ''
        Format          = 'MS profile'
        CreatedFrom     = ''
        XLabel          = 'Time'
        XUnit           = ''
        PrecInX         = '%.2f'
        TimeAxe
        YLabel          = 'Mass'
        YUnit           = ''
        PrecInY         = '%.4f'
        MzAxe
        ZLabel          = 'Intensity'
        ZUnit           = ''
        RndFactor       = 3
        Path2Dat        = {}
        Path2Fin        = ''
        Log             = ''
        BPP
        TIP
        TIS
        FIS
        BIS
        LAST
        ListOfScans
    end
    
    methods
        function obj = Dataset(structInfo)
            if nargin > 1
                obj.Title      	= structInfo.title;
                obj.CreatedFrom = structInfo.CreatedFrom;
                obj.XLabel    	= structInfo.XLabel;
                obj.XUnit    	= structInfo.XUnit;
                obj.YLabel     	= structInfo.YLabel;
                obj.YUnit       = structInfo.YUnit;
                obj.ZLabel     	= structInfo.ZLabel;
                obj.ZUnit       = structInfo.ZUnit;
                obj.Path2Fin    = structInfo.Path2Fin;
                obj.Log         = structInfo.Log;
                
            end
        end
        
        function s = getProfile(obj, mzInt)
            if length(mzInt) == 1
                [mzStart, mzEnd] = deal(mzInt(1));
            else
                mzStart = min(mzInt);
                mzEnd = max(mzInt);
            end
            axeX = obj.MzAxe.Data;
            indMZStt = findCloser(mzStart, axeX);
            indMZEnd = findCloser(mzEnd, axeX);
            
            switch obj.Format
                case 'MS profile'
                    prof(:,1) = obj.TimeAxe.Data;
                    prof(:,2) = 0;
                    
                    h = waitbar(0,'Calculating profile, please wait');
                    for ii = 1:length(prof(:,1))
                        waitbar(ii/length(prof(:,1)))
                        XMS = xpend(obj, obj.ListOfScans{ii});
                        prof(ii,2) = sum(XMS(indMZStt:indMZEnd, 2));
                    end
                    log = 'EIP ';
                    try
                        close(h)
                    catch
                    end
                    
                    structInfo.title = ['Extracted ion Profiles ', ...
                        num2str(axeX(indMZStt), obj.PrecInY),...
                        ' to ', num2str(axeX(indMZEnd), obj.PrecInY), ...
                        ' ', obj.YUnit];
                    
                case 'MS centroid'
                    prof(:,1) = obj.TimeAxe.Data;
                    prof(:,2) = 0;
                    
                    h = waitbar(0,'Calculating profile, please wait');
                    for ii = 1:length(prof(:,1))
                        waitbar(ii/length(prof(:,1)))
                        
                        MS = obj.ListOfScans{ii}.Data;
                        ind2keep = MS(:,1) >= mzStart & MS(:,1) <= mzEnd;
                        
                        if ~isempty(find(ind2keep))
                            prof(ii,2) = sum(MS(ind2keep, 2));
                        else
                            prof(ii,2) = 0;
                        end
                    end
                    log = 'EIP ';
                    close(h);
                    
                    structInfo.title = ['Extracted ion Profiles ', ...
                        num2str(mzStart, obj.PrecInY),...
                        ' to ', num2str(mzEnd, obj.PrecInY), ...
                        ' ', obj.YUnit];
                    
            end
            
            structInfo.traceType    = 'ion profile';
            structInfo.XLabel       = obj.XLabel;
            structInfo.XUnit        = obj.XUnit;
            structInfo.YLabel       = obj.ZLabel;
            structInfo.YUnit        = obj.ZUnit;
            structInfo.Path2Dat     = tempname;
            structInfo.Variables    = 0;
            structInfo.Precision    = 'single';
            structInfo.Path2Fin     = obj.Path2Fin;
            if isempty(obj.Log)
                structInfo.Log = '';
            else
                strLog = decodeLog(obj.Log);
                structInfo.Log = [log, 'DTS=', strLog.DTS, ' MZSTART=', ...
                    num2str(mzStart), ' MZEND=', ...
                    num2str(mzEnd)];
            end
            
            s = Trace(structInfo, prof);
        end
        
        function s = getSpectra(obj, timeInt)
            if length(timeInt) == 1
                [timeStart, timeEnd] = deal(timeInt(1));
            else
                timeStart = min(timeInt);
                timeEnd = max(timeInt);
            end
            axeX       = obj.TimeAxe.Data;
            indTimeStt = findCloser(timeStart, axeX);
            indTimeEnd = findCloser(timeEnd, axeX);
            
            switch obj.Format
                case 'MS profile'
                    data(:,1) = obj.MzAxe.Data;
                    data(:,2) = 0;
                    
                    for ii = indTimeStt:indTimeEnd
                        XMS = xpend(obj, obj.ListOfScans{ii});
                        data(:,2) = data(:,2) + XMS(:,2);
                    end
                    log = 'PRFMS ';
                    
                case 'MS centroid'
                    data = [];
                    
                    for ii = indTimeStt:indTimeEnd
                        MS = obj.ListOfScans{ii}.Data;
                        data = [data; MS];
                        data = sortrows(data, 1);
                    end
                    
                    log = 'CTRMS ';
            end
            
            structInfo.title        = ['MS scan from ', num2str(axeX(indTimeStt), obj.PrecInX),...
                ' to ', num2str(axeX(indTimeEnd), obj.PrecInX), ' ', obj.XUnit];
            switch obj.Format
                case 'MS profile'
                    structInfo.traceType = 'MS profile';
                case 'MS centroid'
                    structInfo.traceType = 'MS centroid';
                otherwise
                    structInfo.traceType = 'MS centroid';
            end
            structInfo.XLabel       = obj.YLabel;
            structInfo.XUnit        = obj.YUnit;
            structInfo.YLabel       = obj.ZLabel;
            structInfo.YUnit        = obj.ZUnit;
            structInfo.Path2Dat     = tempname;
            structInfo.Variables    = 0;
            structInfo.Precision    = 'single';
            structInfo.Path2Fin     = obj.Path2Fin;
            if isempty(obj.Log)
                structInfo.Log = '';
            else
                strLog = decodeLog(obj.Log);
                
                structInfo.Log = [log, 'DTS=', strLog.DTS, ' #START=', ...
                    num2str(indTimeStt), ' #END=', ...
                    num2str(indTimeEnd)];
            end
            s = Trace(structInfo, data);
        end
        
        function XMS = xpend(obj, MS, AxeOr)
            if nargin == 2
                AxeOr = true;
            end
            
            XMS(:,1) = obj.MzAxe.Data;
            if isempty(MS.Data)
                XMS(:,2) = 0;
            else
                cst = MS.Variables;
                rdg = obj.RndFactor;
                axe = MS.Data(:,1)/(1-cst);
                [tf, loc] = ismember(round(axe*10^rdg)/10^rdg, ...
                    round(XMS(:,1)*10^rdg)/10^rdg);
                
                indZeros = find(tf == 0);
                if ~isempty(indZeros)
                    [~, locc] = ismember(ceil(axe*10^rdg)/10^rdg, ...
                        ceil(XMS(:,1)*10^rdg)/10^rdg);
                    loc(indZeros) = locc(indZeros);
                end
                
                indNotNull = loc ~= 0;
                XMS(loc(indNotNull), 2) = MS.Data(indNotNull,2);
            end
            
            if ~AxeOr
                XMS(:,1) = XMS(:,1)*(1-cst);
            end
            
        end
        
        function plot(obj, what)
            figure
            dcm_obj = datacursormode(gcf);
            set(dcm_obj,'UpdateFcn',@myupdatefcn);
            c = uicontextmenu;
            
            switch obj.(what).TraceType
                case 'MS profile'
                    plotObj = plot(obj.(what).Data(:,1), obj.(what).Data(:,2));
                    set(gcf, 'Name', obj.(what).Log);
                    title(obj.(what).Title);
                    xlabel([obj.(what).XLabel, ' / ',obj.(what).XUnit]);
                    ylabel([obj.(what).YLabel, ' / ',  obj.(what).YUnit]);
                    
                    m1 = uimenu(c,'Label','Get Profile At','Callback',@getExtract);
                    m2 = uimenu(c,'Label','Get Profile Between','Callback',@getExtract);
                    
                case 'MS centroid'
                    plotObj = stem(obj.(what).Data(:,1), obj.(what).Data(:,2), ...
                        'Marker', 'none');
                    set(gcf, 'Name', obj.(what).Log);
                    title(obj.(what).Title);
                    xlabel([obj.(what).XLabel, ' / ',obj.(what).XUnit]);
                    ylabel([obj.(what).YLabel, ' / ',  obj.(what).YUnit]);
                    
                    m1 = uimenu(c,'Label','Get Profile At','Callback',@getExtract);
                    m2 = uimenu(c,'Label','Get Profile Between','Callback',@getExtract);
                    
                case 'ion profile'
                    plotObj = plot(obj.(what).Data(:,1), obj.(what).Data(:,2));
                    set(gcf, 'Name', obj.(what).Log);
                    title(obj.(what).Title);
                    xlabel([obj.(what).XLabel, ' / ',obj.(what).XUnit]);
                    ylabel([obj.(what).YLabel, ' / ',  obj.(what).YUnit]);
                    
                    m1 = uimenu(c,'Label','Get Mass Spectrum At','Callback',@getExtract);
                    m2 = uimenu(c,'Label','Get Mass Spectrum Between','Callback',@getExtract);
            end
            m3 =  uimenu(c, 'Label','ExportObj','Callback',@getExtract);
            
            plotObj.UIContextMenu = c;
            
            function getExtract(source, ~)
                switch source.Label
                    case 'Get Mass Spectrum At'
                        [timeInt,~] = ginput(1);
                        obj.LAST = getSpectra(obj, timeInt);
                        plot(obj, 'LAST')
                        
                    case 'Get Mass Spectrum Between'
                        [timeInt,~] = ginput(2);
                        obj.LAST = getSpectra(obj, timeInt);
                        plot(obj, 'LAST')
                        
                    case 'Get Profile At'
                        [mzInt,~] = ginput(1);
                        obj.LAST = getProfile(obj, mzInt);
                        plot(obj, 'LAST')
                        
                    case 'Get Profile Between'
                        [mzInt,~] = ginput(2);
                        obj.LAST = getProfile(obj, mzInt);
                        plot(obj, 'LAST')
                        
                    case 'ExportObj'
                        assignin('base', 'currentObj', obj.LAST)
                        
                end
            end
            
            function txt = myupdatefcn(~,event_obj)
                pos = get(event_obj,'Position');
                
                switch obj.(what).TraceType
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

