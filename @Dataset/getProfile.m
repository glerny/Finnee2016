% BSD 3-Clause License
% Copyright 2016-2017 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function s = getProfile(obj, mzInt, varargin)
options  = checkVarargin( varargin{:});
if length(mzInt) == 1
    [mzStart, mzEnd] = deal(mzInt(1));
else
    mzStart = min(mzInt);
    mzEnd   = max(mzInt);
end
AxisX     = obj.AxisX.Data;
AxisY     = obj.AxisY.Data;

IdStt     = find(AxisX < options.XLim(1), 1, 'last');
if isempty(IdStt), IdStt = 1; end

IdEnd     = find(AxisX > options.XLim(2), 1, 'first');
if isempty(IdEnd), IdEnd = length(AxisX); end

prof(:,1) = obj.AxisX.Data(IdStt:IdEnd);
prof(:,2) = 0;

switch obj.Format
    case 'profile'
        if ~isempty(AxisY)
            indMZStt = findCloser(mzStart, AxisY);
            indMZEnd = findCloser(mzEnd, AxisY);
            
            h = waitbar(0,'Calculating profile, please wait');
            for ii = IdStt:IdEnd
                waitbar(ii/length(prof(:,1)))
                XMS = xpend(obj, obj.ListOfScans{ii});
                
                if isnan(options.setSz)
                    % Modification 21/07/2017
                    % calculate the extracted profile based on the set
                    % limits
                    
                    prof(ii - IdStt +1,2) = sum(XMS.Data(indMZStt:indMZEnd, 2));
                else
                    % Modification 21/07/2017
                    % calculate the extracted profile using the max
                    % intensity between the set values and a set number of
                    % intervals 
                    
                    [ItMax, IdMax] = max(XMS.Data(indMZStt:indMZEnd, 2));
                        Sz    = options.setSz;
                    if ItMax ~=0
                        IdMax = IdMax(1) + indMZStt - 1;
                        prof(ii - IdStt +1,2) = sum(XMS.Data(IdMax - Sz:IdMax + Sz, 2));
                    end
                    
                end
            end
            
            if isnan(options.setSz)
                InfoTrc.Title = ['Extracted ion Profiles from ', ...
                    num2str(AxisY(indMZStt), obj.AxisY.fo),...
                    ' to ', num2str(AxisY(indMZEnd), obj.AxisY.fo), ...
                    ' ', obj.AxisY.Unit];
            else
                InfoTrc.Title = ['Extracted ion Profiles: max between ', ...
                    num2str(mzStart, obj.AxisY.fo)                     ,...
                    ' and ', num2str(mzEnd, obj.AxisY.fo)              , ...
                    ' ', obj.AxisY.Unit                                ,...
                    ' +/- ', num2str(Sz), ' interval(s)'];
            end
        else
            h = waitbar(0,'Calculating profile, please wait');
            for ii = IdStt:IdEnd
                waitbar(ii/length(prof(:,1)))
                MS = obj.ListOfScans{ii};
                indMZStt = findCloser(mzStart, MS.Data(:,1));
                indMZEnd = findCloser(mzEnd,   MS.Data(:,1));
                if isnan(options.setSz)
                    % Modification 21/07/2017
                    % calculate the extracted profile based on the set
                    % limits
                    
                    prof(ii - IdStt +1,2) = sum(MS.Data(indMZStt:indMZEnd, 2));
                else
                    % Modification 21/07/2017
                    % calculate the extracted profile using the max
                    % intensity between the set values and a set number of
                    % intervals 
                    
                    [ItMax, IdMax] = max(MS.Data(indMZStt:indMZEnd, 2));
                        Sz    = options.setSz;
                    if ItMax ~=0
                        IdMax = IdMax(1) + indMZStt - 1;
                        Id1 = max(1, IdMax - Sz);
                        Id2 = min(size(MS.Data, 1), IdMax + Sz);
                        prof(ii - IdStt +1,2) = sum(MS.Data(Id1:Id2, 2));
                    end
                    
                end
            end
            
            if isnan(options.setSz)
                InfoTrc.Title = ['Extracted ion Profiles from ', ...
                    num2str(mzStart, obj.AxisY.fo)             ,...
                    ' to ', num2str(mzEnd, obj.AxisY.fo)       , ...
                    ' ', obj.AxisY.Unit];
            else
                InfoTrc.Title = ['Extracted ion Profiles: max between ', ...
                    num2str(mzStart, obj.AxisY.fo)                     ,...
                    ' and ', num2str(mzEnd, obj.AxisY.fo)              , ...
                    ' ', obj.AxisY.Unit                                ,...
                    ' +/- ', num2str(Sz), ' interval(s)'];
            end
        end
        
    case 'centroid'
        h = waitbar(0,'Calculating profile, please wait');
        for ii = IdStt:IdEnd
            waitbar(ii/length(prof(:,1)))
            MS = obj.ListOfScans{ii}.Data;
            InfoTrc.Title = ['Extracted ion Profiles from ', ...
                num2str(mzStart, obj.AxisY.fo),...
                ' to ', num2str(mzEnd, obj.AxisY.fo),...
                ' ', obj.AxisY.Unit];
            if isempty(MS)
                prof(ii - IdStt +1,2) = 0;
                continue
            end
            
            ind2keep = MS(:,1) >= mzStart & MS(:,1) <= mzEnd;
            if any(ind2keep)
                prof(ii - IdStt +1,2) = sum(MS(ind2keep, 2));
            else
                prof(ii - IdStt +1,2) = 0;
            end
        end
end

try close(h); catch, end %#ok<CTCH>

InfoTrc.TT     = 'SEP';
[~, partial]   = decipherLog(obj.Log);
InfoTrc.FT     = partial{1};
InfoTrc.AxisX  = Axis(obj.AxisX.InfoAxis);
InfoTrc.AxisY  = Axis(obj.AxisZ.InfoAxis);
InfoTrc.Loc    = 'inTrace';
InfoTrc.AdiPrm = {};
InfoTrc.P2Fin  = obj.InfoDts.P2F;

s = Trace(InfoTrc, prof);

%% CHECKVARARGIN
    function options = checkVarargin(varargin)
        % CHECKVARARGIN is used to check the input paramters
        % and create the options parameter.
        
        options.XLim  = [0 inf];
        options.setSz = NaN;
        
        % Decipher varargin
        input = @(x) find(strcmpi(varargin,x),1);
        
        tgtIx = input('XLim');
        if ~isempty(tgtIx)
            XmM = varargin{tgtIx +1};
            options.XLim(1) = min(XmM);
            options.XLim(2) = max(XmM);
        end
        
        tgtIx = input('SetSize');
        if ~isempty(tgtIx)
            options.setSz = varargin{tgtIx +1};
        end
    end
end
