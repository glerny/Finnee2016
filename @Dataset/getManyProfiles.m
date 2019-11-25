% BSD 3-Clause License
% Copyright 2016-2017 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function s = getManyProfiles(obj, mzList, Dx, timeList, interval, varargin)
options  = checkVarargin( varargin{:});
AxisX     = obj.AxisX.Data;
AxisY     = obj.AxisY.Data;

IdStt     = find(AxisX < options.XLim(1), 1, 'last');
if isempty(IdStt), IdStt = 1; end

IdEnd     = find(AxisX > options.XLim(2), 1, 'first');
if isempty(IdEnd), IdEnd = length(AxisX); end

prof(:,1) = obj.AxisX.Data(IdStt:IdEnd);
prof(:,2:length(mzList)+1) = 0;

switch obj.Format
    case 'profile'
        if ~isempty(AxisY)
            
            for jj = 1:length(mzList)
                ix = findCloser(mzList(jj), AxisY);
                indMZ(jj, 1) = max(1, ix-Dx);
                indMZ(jj, 2) = min(ix+Dx, length(AxisY));
            end
            
            h = waitbar(0,'Calculating profile, please wait');
            for ii = IdStt:IdEnd
                waitbar(ii/length(prof(:,1)))
%                 try
%                     XMS = xpend(obj, obj.ListOfScans{ii});
%                 catch
                    obj.Path2Dat = strrep(obj.Path2Dat,'G:\O meu disco\', 'G:\My Drive\');
                    XMS = xpend(obj, obj.ListOfScans{ii});
%                 end
                
                if isnan(options.setSz)
                    % Modification 21/07/2017
                    % calculate the extracted profile based on the set
                    % limits
                    
                    for jj = 1:length(mzList)
                        prof(ii - IdStt +1,jj+1) = ...
                            max( XMS.Data(indMZ(jj,1):indMZ(jj,2), 2));
                    end
                    
                else
                    % Modification 21/07/2017
                    % calculate the extracted profile using the max
                    % intensity between the set values and a set number of
                    % intervals 
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%TO BE DONE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     [ItMax, IdMax] = max(XMS.Data(indMZStt:indMZEnd, 2));
%                         Sz    = options.setSz;
%                     if ItMax ~=0
%                         IdMax = IdMax(1) + indMZStt - 1;
%                         prof(ii - IdStt +1,2) = sum(XMS.Data(IdMax - Sz:IdMax + Sz, 2));
%                     end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%TO BE DONE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                end
            end
            
            if isnan(options.setSz)
                for jj = 1:length(mzList)
                    InfoTrc{jj}.Title = ['Extracted ion Profiles from ', ...
                        num2str(AxisY(indMZ(jj,1)), obj.AxisY.fo),...
                        ' to ', num2str(AxisY(indMZ(jj,2)), obj.AxisY.fo), ...
                        ' ', obj.AxisY.Unit];
                end
            else
%                 InfoTrc.Title = ['Extracted ion Profiles: max between ', ...
%                     num2str(mzStart, obj.AxisY.fo)                     ,...
%                     ' and ', num2str(mzEnd, obj.AxisY.fo)              , ...
%                     ' ', obj.AxisY.Unit                                ,...
%                     ' +/- ', num2str(Sz), ' interval(s)'];
            end
        else
%             h = waitbar(0,'Calculating profile, please wait');
%             for ii = IdStt:IdEnd
%                 waitbar(ii/length(prof(:,1)))
%                 MS = obj.ListOfScans{ii};
%                 indMZStt = findCloser(mzStart, MS.Data(:,1));
%                 indMZEnd = findCloser(mzEnd,   MS.Data(:,1));
%                 if isnan(options.setSz)
%                     % Modification 21/07/2017
%                     % calculate the extracted profile based on the set
%                     % limits
%                     
%                     try
%                     prof(ii - IdStt +1,2) = sum(MS.Data(indMZStt:indMZEnd, 2));
%                     catch
%                         diap('wtf')
%                     end
%                 else
%                     % Modification 21/07/2017
%                     % calculate the extracted profile using the max
%                     % intensity between the set values and a set number of
%                     % intervals 
%                     
%                     [ItMax, IdMax] = max(MS.Data(indMZStt:indMZEnd, 2));
%                         Sz    = options.setSz;
%                     if ItMax ~=0
%                         IdMax = IdMax(1) + indMZStt - 1;
%                         Id1 = max(1, IdMax - Sz);
%                         Id2 = min(size(MS.Data, 1), IdMax + Sz);
%                         prof(ii - IdStt +1,2) = sum(MS.Data(Id1:Id2, 2));
%                     end
%                     
%                 end
%             end
%             
%             if isnan(options.setSz)
%                 InfoTrc.Title = ['Extracted ion Profiles from ', ...
%                     num2str(mzStart, obj.AxisY.fo)             ,...
%                     ' to ', num2str(mzEnd, obj.AxisY.fo)       , ...
%                     ' ', obj.AxisY.Unit];
%             else
%                 InfoTrc.Title = ['Extracted ion Profiles: max between ', ...
%                     num2str(mzStart, obj.AxisY.fo)                     ,...
%                     ' and ', num2str(mzEnd, obj.AxisY.fo)              , ...
%                     ' ', obj.AxisY.Unit                                ,...
%                     ' +/- ', num2str(Sz), ' interval(s)'];
%             end
        end
        
    case 'centroid'
        
        h = waitbar(0,'Calculating profile, please wait');
        
        mzInt = mzList -  Dx*1e-6*mzList;
        mzInt(:,2) = mzList +  Dx*1e-6*mzList;
            
        for ii = IdStt:IdEnd
            waitbar(ii/length(prof(:,1)))
            
            for jj = 1:length(mzList)
                MS = obj.ListOfScans{ii}.Data;
                InfoTrc.Title = ['Extracted ion Profiles from ', ...
                    num2str(mzInt(jj,1), obj.AxisY.fo),...
                    ' to ', num2str(mzInt(jj,2), obj.AxisY.fo),...
                    ' ', obj.AxisY.Unit];
                
                if isempty(MS)
                    prof(ii - IdStt +1, jj+1) = 0;
                    continue
                end
                
                ind2keep = MS(:,1) >= mzInt(jj,1) & MS(:,1) <= mzInt(jj,2);
                
                if any(ind2keep)
                    prof(ii - IdStt +1, jj+1) = sum(MS(ind2keep, 2));
                else
                    prof(ii - IdStt +1, jj+1) = 0;
                end
            end
        end
end

try close(h); catch, end %#ok<CTCH>

prof(isnan(prof)) = 0;
for jj = 1:length(mzList)
    InfoTrc{jj}.TT     = 'SEP';
    [~, partial]   = decipherLog(obj.Log);
    InfoTrc{jj}.FT     = partial{1};
    InfoTrc{jj}.AxisX  = Axis(obj.AxisX.InfoAxis);
    InfoTrc{jj}.AxisY  = Axis(obj.AxisZ.InfoAxis);
    InfoTrc{jj}.Loc    = 'inTrace';
    InfoTrc{jj}.AdiPrm = {};
    InfoTrc{jj}.P2Fin  = obj.InfoDts.P2F;
    
    IdStt     = find(prof(:,1) <  timeList(jj) - interval, 1, 'last');
    if isempty(IdStt), IdStt = 1; end

    IdEnd     = find(prof(:,1) > timeList(jj) + interval, 1, 'first');
    if isempty(IdEnd), IdEnd = length(prof(:,1)); end
    
    clear newTitle
    newTitle{1} =  InfoTrc{jj}.Title;
    newTitle{2} = sprintf('Migration time %.2f min', timeList(jj));
    InfoTrc{jj}.Title = newTitle;
    s{jj} = Trace(InfoTrc{jj}, prof(IdStt:IdEnd, [1, jj+1]));
end

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
