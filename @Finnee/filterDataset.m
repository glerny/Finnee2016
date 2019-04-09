%% DESCRIPTION
% the FILTERDATASET method is used to filter a dataset using a particular
% method and create a new dataset with the filtered scans. FILTERDATASET is
% used in particular to remove the spikes or the background noise. For more
% information see
% <https://github.com/glerny/Finnee2016/wiki/Basic-operation-with-Dataset>
%
%% INPUT PARAMETERS
% *Compulsory
%   _obj_     : The Finnee object
%   _dts_     : The indices to the dataset to correct
%       (i.e. myFinnee.Datasets{dts}
%   _method_  : The method to implement in the format
%   'methodName:param1:param2:...:paramn'
%       + 'RemoveSpikes:pm1'
%           Allow to remove spikes in each MS spectra where spikes are
%           defined as any series of points delimited by null values whose
%           number of points is lower or equal to pm1.
%       + 'thresholding:pm1:pm2:pm3:pm4'
%           Use to remove value below a set threshold. For each points, will
%           reconstitued a small matrix, centered around the point of interest
%           and of size 2*pm1+1 in the time dimension and 2*pm2+1 in the
%           m/z dimension. The center point intensity will be set to zero
%           if no values within the matrix is higher than pm3. 
%       + 'RemoveNoise:pm1:pm2:pm3:pm4'
%           Use to remove background noise. For each points, will
%           reconstitued a small matrix, centered around the point of interest
%           and of size 2*pm1+1 in the time dimension and 2*pm2+1 in the
%           m/z dimension. The center point intensity will be set to zero
%           if no values within the matrix has a signal to noise higher than
%           pm3. The function only works in the backgroun noise has been
%           estimated, pm4 defined how the noise is calculated in m/z were
%           no data are available. If pm4 start with % (%5 or %25 for
%           example), missing dat are the mean of the percentil, otherwise
%           it is a set value.
%   - varargin: Possible options
%       + 'spikes:pm1'
%          Only actives with 'RemoveNoise:pm1:pm2:pm3'. Will also performe
%          a spikes removal filter after the noise removal filter.
%
%% Copyright
% BSD 3-Clause License
% Copyright 2016-2017 G. Erny (guillaume@fe.up.pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function obj = filterDataset(obj, dts, method, varargin)

%% CORE OF THE FUNCTION
% 1- Initialisation and options
narginchk(3, inf)
dtsIn            = obj.Datasets{dts};
infoDts          = dtsIn.InfoDts;
infoDts.Path2Dat = {};
[options, MtU]   = checkVarargin(infoDts, dts, method, varargin{:});

% Create new dat file
[~, rndStr]         = fileparts(tempname);
allProfiles(:,1)    = dtsIn.AxisX.Data;
allProfiles(:,3)    = 0;
infoDts.ListOfScans = {};
AxisMZ              = dtsIn.AxisY.Data;
AxisMZ(:,4)         = 0;
fln                 = 1;
m                   = length(obj.Datasets)+1;

% 2- Load each scan and run the filter algorithm
switch MtU{1}
    case 'RemoveSpikes'
        %% REMOVESPIKES HERE
        infoDts.Title = 'Spikes removed dataset';
        infoDts.Log   = ['PRF=', num2str(m), ' SPR=', options.method,...
            '|', infoDts.Log];
        h = waitbar(0,'processing scans');
        for ii = 1:size(allProfiles, 1)
            waitbar(ii/size(allProfiles, 1))
            
            infoDts.Path2Dat{fln} = fullfile(obj.Path2Fin, rndStr);
            Scanii  = dtsIn.ListOfScans{ii};
            XMS     = dtsIn.xpend(Scanii.filterTrace(options.method));
            infoScn = Scanii.InfoTrc;
            Log     = decipherLog(infoDts.Log, 1);
            infoScn.FT        = Log{1};
            infoScn.Path2Dat  = infoDts.Path2Dat{fln};
            infoScn.Loc       = 'inFile';
            infoScn.Precision = 'single';
            AxisMZ(:,2)        = AxisMZ(:,2) + XMS.Data(:,2);
            iNZ = XMS.Data(:,2) > 0;
            AxisMZ(iNZ, 3)     = AxisMZ(iNZ, 3) + 1;
            AxisMZ(:,4)        = max([AxisMZ(:,4), XMS.Data(:,2)], [], 2);
            allProfiles(ii,2) = sum(XMS.Data(:,2));
            allProfiles(ii,3) = max(XMS.Data(:,2));
            
            % reduced trailing zero in excess
            provMat      = [XMS.Data(2:end, 2); 0];
            provMat(:,2) = XMS.Data(:, 2);
            provMat(:,3) = [0; XMS.Data(1:end-1, 2)];
            MS           = XMS.Data(sum(provMat, 2) > 0, :);
            
            % recorded each scans
            if isempty(MS)
                infoDts.ListOfScans{ii} = Trace(infoScn);
            else
                infoDts.ListOfScans{ii} = Trace(infoScn, MS);
            end
            
            
            s = dir(infoDts.Path2Dat{fln});
            if s.bytes > obj.Options.MaxFileSize;
                [~, rndStr]           = fileparts(tempname);
                fln                   = fln + 1;
                infoDts.Path2Dat{fln} = fullfile(obj.Path2Fin, rndStr);
                infoScn.Path2Dat      = infoDts.Path2Dat{fln};
            end
        end
        
    case 'threshold'
        %% THRESHOLD HERE
        infoDts.Title = 'Thresholded dataset';
        infoDts.Log   = ['PRF=', num2str(m), ' THR =', options.method,...
            '|', infoDts.Log];
        
        S2NMat  = zeros(size(AxisMZ, 1), 2*MtU{2}+1);
        for ii  = 1:MtU{2}+1
            Scanii  = dtsIn.ListOfScans{ii};
            if ~isempty(Scanii.Data)
                nnMS = Scanii.Data(:,2);
                [~, loc] = ismember(Scanii.Data(:,1), AxisMZ(:,1));
                IZ = loc~= 0;
                S2NMat(loc(IZ), MtU{2}+ii)  =  nnMS(IZ);
            else
                S2NMat(:, MtU{2}+ii)  = 0;
            end
        end
        ZB      = zeros(MtU{3}, 2*MtU{2}+1);
        S2NMat  = [ZB; S2NMat; ZB];
        
        h = waitbar(0,'processing scans');
        for ii = 1:length(dtsIn.ListOfScans)
            waitbar(ii/length(dtsIn.ListOfScans))
			
			TM = [];
			for jj = -MtU{3}:1:MtU{3}
				TM = [TM, circshift(S2NMat(:, MtU{2}+1), [jj, 0])]; %#ok<*AGROW>
            end
			ind2null = max(TM, [], 2) < MtU{4} & max(S2NMat, [], 2) < MtU{4};
            infoDts.Path2Dat{fln} = fullfile(obj.Path2Fin, rndStr);
            MS2Cor    = AxisMZ(:,1);
            Scanii = dtsIn.ListOfScans{ii};
            
            if ~isempty(Scanii.Data)
                [~, loc] = ismember(Scanii.Data(:,1), MS2Cor(:,1));
                IZ = loc~= 0;
                MS2Cor(loc(IZ), 2) = Scanii.Data(IZ,2);
                MS2Cor(ind2null(MtU{3}+1:end-MtU{3}), 2) = 0;
            else
                MS2Cor(:,2) = 0;
            end
            
            % find and remove spikes
            if options.RemSpks
                spkSz  = options.SpkSz;
                MS2Cor = spikesRemoval(MS2Cor, spkSz );
            end
            
            infoScn           = Scanii.InfoTrc;
            [~, partial]      = decipherLog(infoDts.Log);
            
            infoScn.FT        = partial{1};
            infoScn.Path2Dat  = infoDts.Path2Dat{fln};
            infoScn.Loc       = 'inFile';
            infoScn.Precision = 'single';
            AxisMZ(:,2)       = AxisMZ(:,2) + MS2Cor(:,2);
            iNZ = MS2Cor(:,2) > 0;
            AxisMZ(iNZ, 3)    = AxisMZ(iNZ, 3) + 1;
            AxisMZ(:,4)       = max([AxisMZ(:,4), MS2Cor(:,2)], [], 2);
            allProfiles(ii,2) = sum(MS2Cor(:,2));
            allProfiles(ii,3) = max(MS2Cor(:,2));
            
            % reduced trailing zero in excess
            MS2Cor = trailRem(MS2Cor, 2);
            
            % recorded each scans
            if isempty(MS2Cor)
                infoDts.ListOfScans{ii} = Trace(infoScn);
            else
                infoDts.ListOfScans{ii} = Trace(infoScn, MS2Cor);
            end
            
            % Load the  next scan
            if ii+MtU{2}+1 > length(dtsIn.ListOfScans)
                S2A  = AxisMZ(:,1);
                S2A(:, 2) = 0;
            else
                scan2add = dtsIn.ListOfScans{ii+MtU{2}+1};
                S2A  = AxisMZ(:,1);
                S2A(:, 2) = 0;
                if ~isempty(scan2add.Data)
                    nnMS = scan2add.Data(:,2);
                    [~, loc] = ismember(scan2add.Data(:,1), AxisMZ(:,1));
                    IZ = loc~= 0;
                    S2A(loc(IZ), 2) = nnMS(IZ);
                end
            end
            scan2add = [ZB(:,1); S2A(:,2); ZB(:,1)];
            S2NMat   = [S2NMat(:, 2:end), scan2add];

            
            s = dir(infoDts.Path2Dat{fln});
            if isempty(s), continue, end
            if s.bytes > obj.Options.MaxFileSize
                [~, rndStr]           = fileparts(tempname);
                fln                   = fln + 1;
                infoDts.Path2Dat{fln} = fullfile(obj.Path2Fin, rndStr);
                infoScn.Path2Dat      = infoDts.Path2Dat{fln};
            end
        end
        
        case 'remNoise'
        %% REMNOISE HERE
        infoDts.Title = 'Noise corrected dataset';
        infoDts.Log   = ['PRF=', num2str(m), ' NCR=', options.method,...
            '|', infoDts.Log];
        Noise   =  obj.Datasets{dts}.AddInfo.Peak2PeakNoise;
        nd = Noise.Data;
        
        ZerosNoise = MtU{5};
        if ZerosNoise(1) == '%'
            val = max(1, str2double(ZerosNoise(2:end)));
            minNoise = round(mean(nd(nd(:,2) <= ...
                prctile(nd(nd(:,2) > 0, 2), val) & nd(:,2) > 0, 2)))
        else
            val = str2num(ZerosNoise);
            minNoise = val;
        end
        nd(nd(:,2) <= minNoise, 2) = minNoise;
        
        S2NMat  = zeros(size(AxisMZ, 1), 2*MtU{2}+1);
        for ii  = 1:MtU{2}+1
            Scanii  = dtsIn.ListOfScans{ii};
            if ~isempty(Scanii.Data)
                [~, loc] = ismember(Scanii.Data(:,1), nd(:,1));
                nnMS = Scanii.Data(:,2)./nd(loc,2);
                [~, loc] = ismember(Scanii.Data(:,1), AxisMZ(:,1));
                IZ = loc~= 0;
                S2NMat(loc(IZ), MtU{2}+ii)  =  nnMS(IZ);
            else
                S2NMat(:, MtU{2}+ii)  = 0;
            end
        end
        ZB      = zeros(MtU{3}, 2*MtU{2}+1);
        S2NMat  = [ZB; S2NMat; ZB];
        
        h = waitbar(0,'processing scans');
        for ii = 1:length(dtsIn.ListOfScans)
            waitbar(ii/length(dtsIn.ListOfScans))
			
			TM = [];
			for jj = -MtU{3}:1:MtU{3}
				TM = [TM, circshift(S2NMat(:, MtU{2}+1), [jj, 0])]; %#ok<*AGROW>
            end
			ind2null = max(TM, [], 2) < MtU{4} & max(S2NMat, [], 2) < MtU{4};
            infoDts.Path2Dat{fln} = fullfile(obj.Path2Fin, rndStr);
            MS2Cor    = AxisMZ(:,1);
            Scanii = dtsIn.ListOfScans{ii};
            
            if ~isempty(Scanii.Data)
                [~, loc] = ismember(Scanii.Data(:,1), MS2Cor(:,1));
                IZ = loc~= 0;
                MS2Cor(loc(IZ), 2) = Scanii.Data(IZ,2);
                MS2Cor(ind2null(MtU{3}+1:end-MtU{3}), 2) = 0;
                MS2Cor(MS2Cor(:,2) < 0, 2) = 0;
            else
                MS2Cor(:,2) = 0;
            end
            
            
            % find and remove spikes
            if options.RemSpks
                spkSz  = options.SpkSz;
                MS2Cor = spikesRemoval(MS2Cor, spkSz );
            end
            
            infoScn           = Scanii.InfoTrc;
            [~, partial]      = decipherLog(infoDts.Log);
            
            infoScn.FT        = partial{1};
            infoScn.Path2Dat  = infoDts.Path2Dat{fln};
            infoScn.Loc       = 'inFile';
            infoScn.Precision = 'single';
            AxisMZ(:,2)       = AxisMZ(:,2) + MS2Cor(:,2);
            iNZ = MS2Cor(:,2) > 0;
            AxisMZ(iNZ, 3)    = AxisMZ(iNZ, 3) + 1;
            AxisMZ(:,4)       = max([AxisMZ(:,4), MS2Cor(:,2)], [], 2);
            allProfiles(ii,2) = sum(MS2Cor(:,2));
            allProfiles(ii,3) = max(MS2Cor(:,2));
            
            % reduced trailing zero in excess
            MS2Cor = trailRem(MS2Cor, 2);
            
            % recorded each scans
            if isempty(MS2Cor)
                infoDts.ListOfScans{ii} = Trace(infoScn);
            else
                infoDts.ListOfScans{ii} = Trace(infoScn, MS2Cor);
            end
            
            % Load the  next scan
            if ii+MtU{2}+1 > length(dtsIn.ListOfScans)
                S2A  = AxisMZ(:,1);
                S2A(:, 2) = 0;
            else
                scan2add = dtsIn.ListOfScans{ii+MtU{2}+1};
                S2A  = AxisMZ(:,1);
                S2A(:, 2) = 0;
                if ~isempty(scan2add.Data)
                    [~, loc] = ismember(scan2add.Data(:,1), nd(:,1));
                    nnMS = scan2add.Data(:,2)./nd(loc,2);
                    [~, loc] = ismember(scan2add.Data(:,1), AxisMZ(:,1));
                IZ = loc~= 0;
                    S2A(loc(IZ), 2) = nnMS(IZ);
                end
            end
            scan2add = [ZB(:,1); S2A(:,2); ZB(:,1)];
            S2NMat   = [S2NMat(:, 2:end), scan2add];

            
            s = dir(infoDts.Path2Dat{fln});
            if isempty(s), continue, end
            if s.bytes > obj.Options.MaxFileSize
                [~, rndStr]           = fileparts(tempname);
                fln                   = fln + 1;
                infoDts.Path2Dat{fln} = fullfile(obj.Path2Fin, rndStr);
                infoScn.Path2Dat      = infoDts.Path2Dat{fln};
            end
        end
        
    case 'SGolay'
        %% STAVINSKY GOLAY FILTERING HERE
        infoDts.Title = 'SGolay filtered dataset'; %%% One d SGolay through time, what about 2D?
        infoDts.Log   = ['PRF=', num2str(m), ' SGF=', options.method,...
            '|', infoDts.Log];
        infoDts.Path2Dat{fln} = fullfile(obj.Path2Fin, rndStr);
        dt2keep = false(length(dtsIn.AxisX.Data), 1);
        ii = 1:MtU{3}:length(dtsIn.AxisX.Data);
        dt2keep(ii) = true;
        
        SGFMat  = zeros(size(AxisMZ, 1), MtU{2});
        for ii  = 1:(MtU{2}+1)/2
            Scanii  = dtsIn.ListOfScans{ii};
            if ~isempty(Scanii.Data)
                [~, loc] = ismember(Scanii.Data(:,1), AxisMZ(:,1));
                SGFMat(loc, ii+(MtU{2}-1)/2)  =  Scanii.Data(:,2);
            else
                SGFMat(:,ii+(MtU{2}-1)/2)  = 0;
            end
        end
        
        h = waitbar(0,'processing scans');
        for ii = 1:length(dtsIn.ListOfScans)
            waitbar(ii/length(dtsIn.ListOfScans))
            MS2Cor      = AxisMZ(:,1);
            SGD = (sgolayfilt(SGFMat', 1, MtU{2}))';
            MS2Cor      = AxisMZ(:,1);
            MS2Cor(:,2) = SGD(:, (MtU{2}+1)/2);%round(SGD(:, (MtU{2}+1)/2));
            
            % Load the  next scan
                if ii+(MtU{2}-1)/2 <= length(dtsIn.ListOfScans)
                    scan2add = dtsIn.ListOfScans{ii+(MtU{2}-1)/2};
                    if ~isempty(scan2add.Data)
                        [~, loc] = ismember(scan2add.Data(:,1), AxisMZ(:,1));
                        SGFMat(loc, end+1)  =  scan2add.Data(:,2);
                    else
                        SGFMat(:,end+1)  = 0;
                    end
                else
                    SGFMat(:,end+1)  = 0;
                end
                SGFMat   =  SGFMat(:, 2:end);
                
            if dt2keep(ii)
                % find and remove spikes
                
                % MS2Cor(MS2Cor(:,2) < 0, 2) = 0;
                
                
                if options.RemSpks
                    spkSz  = options.SpkSz;
                    MS2Cor = spikesRemoval(MS2Cor, spkSz );
                end
                
                infoScn           = Scanii.InfoTrc;
                [~, partial]      = decipherLog(infoDts.Log);
                
                infoScn.FT        = partial{1};
                infoScn.Path2Dat  = infoDts.Path2Dat{fln};
                infoScn.Loc       = 'inFile';
                infoScn.Precision = 'single';
                AxisMZ(:,2)       = AxisMZ(:,2) + MS2Cor(:,2);
                iNZ = MS2Cor(:,2) > 0;
                AxisMZ(iNZ, 3)    = AxisMZ(iNZ, 3) + 1;
                AxisMZ(:,4)       = max([AxisMZ(:,4), MS2Cor(:,2)], [], 2);
                allProfiles(ii,2) = sum(MS2Cor(:,2));
                allProfiles(ii,3) = max(MS2Cor(:,2));
                
                % reduced trailing zero in excess
                MS2Cor = trailRem(MS2Cor, 2);
                
                % recorded each scans
                if isempty(MS2Cor)
                    infoDts.ListOfScans{ii} = Trace(infoScn);
                else
                    infoDts.ListOfScans{ii} = Trace(infoScn, MS2Cor);
                end
                
                
                s = dir(infoDts.Path2Dat{fln});
                if isempty(s), continue, end
                if s.bytes > obj.Options.MaxFileSize;
                    [~, rndStr]           = fileparts(tempname);
                    fln                   = fln + 1;
                    infoDts.Path2Dat{fln} = fullfile(obj.Path2Fin, rndStr);
                    infoScn.Path2Dat      = infoDts.Path2Dat{fln};
                end
            end
        end
        allProfiles   = allProfiles(dt2keep, :);
        infoDts.ListOfScans = infoDts.ListOfScans(dt2keep);
    otherwise
        error('%s is not a recognised method', MtU{1})
end
try close(h), catch, end

% creating Dataset
% reduced trailing zero in AxisMZ to
AxisMZ        = trailRem(AxisMZ, 2);

infoAxis           = dtsIn.AxisX.InfoAxis;
infoAxis.Loc       = 'inFile';
infoAxis.Precision = 'single';
infoAxis.Path2Dat  = infoDts.Path2Dat{fln};

if isempty(allProfiles)
    infoDts.AxisX  = Axis(infoAxis);
else
    infoDts.AxisX  = Axis(infoAxis, allProfiles(:,1));
end

infoAxis           = dtsIn.AxisY.InfoAxis;
infoAxis.Loc       = 'inFile';
infoAxis.Precision = 'single';
infoAxis.Path2Dat  = infoDts.Path2Dat{fln};
if isempty(AxisMZ)
    infoDts.AxisY  = Axis(infoAxis);
else
    infoDts.AxisY  = Axis(infoAxis, AxisMZ(:,1));
end

infoPrf.AxisX      = Axis(dtsIn.AxisX.InfoAxis);
infoPrf.AxisY      = Axis(dtsIn.AxisZ.InfoAxis);
infoPrf.Loc       = 'inFile';
infoPrf.Precision = 'single';
infoPrf.Path2Dat  = infoDts.Path2Dat{fln};
infoPrf.FT        = infoDts.Log;
infoPrf.TT        = 'SEP';
infoPrf.AdiPrm    = {};
infoPrf.P2Fin     = obj.Path2Fin;

infoPrf.Title     = 'Base Peak Profiles';
if isempty(allProfiles)
    infoDts.BPP   = Trace(infoPrf);
else
    infoDts.BPP   = Trace(infoPrf, [allProfiles(:,1), allProfiles(:,3)]);
end

infoPrf.Title     = 'Total Ion profiles';
if isempty(allProfiles)
    infoDts.TIP   = Trace(infoPrf);
else
    infoDts.TIP   = Trace(infoPrf, [allProfiles(:,1), allProfiles(:,2)]);
end

infoPrf.AxisX      = Axis(dtsIn.AxisY.InfoAxis);
infoPrf.Title     = 'Total Ion Spectrum';
if isempty(AxisMZ)
    infoDts.TIS   = Trace(infoPrf);
else
    infoDts.TIS   = Trace(infoPrf, [AxisMZ(:,1), AxisMZ(:,2)]);
end

infoPrf.Title     = 'Frequency Ion Spectrum';
if isempty(AxisMZ)
    infoDts.FIS   = Trace(infoPrf);
else
    infoDts.FIS   = Trace(infoPrf, [AxisMZ(:,1), AxisMZ(:,3)]);
end

infoPrf.Title     = 'Base Ion Spectrum';
if isempty(AxisMZ)
    infoDts.BIS   = Trace(infoPrf);
else
    infoDts.BIS   = Trace(infoPrf, [AxisMZ(:,1), AxisMZ(:,4)]);
end

infoDts.LAST      = Trace();

infoDts.Option4crt.function = 'filterDataset';
infoDts.Option4crt.Options  = options;
obj.Datasets{end+1}         = Dataset(infoDts);
obj.save;

%% SUB FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CHECKVARARGIN
    function [options, MtU] = checkVarargin(infoDts, dts, method, varargin)
        % CHECKVARARGIN is used to check the input paramters
        % and create the options parameter.
        
        % 1- Verify the method and its parameters
        MtU = strsplit(method, ':');
        switch lower(MtU{1})
            case 'removespikes'
                % Check Dataset format
                if ~strcmp(infoDts.Format, 'profile')
                    error('The original dataset should be in profile mode, dataset %i is in %f mode',...
                        dts, infoDts.Format);
                end
                
                % check parameter for method
                MtU{1} = 'RemoveSpikes';
                if length(MtU) == 1
                    MtU{2} = 2;
                else
                    MtU{2} = str2double(MtU{2});
                end
                options.method = ['RemoveSpikes:', num2str(MtU{2})];
                
            case 'thresholding'
                % Check Dataset format
                if ~strcmp(infoDts.Format, 'profile')
                    error('The original dataset should be in profile mode, dataset %i is in %f mode',...
                        dts, infoDts.Format);
                end
                
                % check parameter for method
                MtU{1} = 'threshold';
                if length(MtU) < 2
                    MtU{2} = 10;
                else
                    MtU{2} = str2double(MtU{2});
                end
                
                if length(MtU) < 3
                    MtU{3} = 3;
                else
                    MtU{3} = str2double(MtU{3});
                end
                
                if length(MtU) < 4
                    MtU{4} = 1000;
                else
                    MtU{4} = str2double(MtU{4});
                end
                
                
                options.method = ['thresholding:', num2str(MtU{2}),':',...
                    num2str(MtU{3}),':', num2str(MtU{4})];
                
            case 'remnoise'
                % Check Dataset format
                if ~strcmp(infoDts.Format, 'profile')
                    error('The original dataset should be in profile mode, dataset %i is in %f mode',...
                        dts, infoDts.Format);
                end
                
                % check parameter for method
                MtU{1} = 'remNoise';
                if length(MtU) < 2
                    MtU{2} = 10;
                else
                    MtU{2} = str2double(MtU{2});
                end
                
                if length(MtU) < 3
                    MtU{3} = 3;
                else
                    MtU{3} = str2double(MtU{3});
                end
                
                if length(MtU) < 4
                    MtU{4} = 10;
                else
                    MtU{4} = str2double(MtU{4});
                end
                
                if length(MtU) < 5
                    MtU{5} = '%5';
                end
                
                options.method = ['RemoveNoise:', num2str(MtU{2}),':',...
                    num2str(MtU{3}),':', num2str(MtU{4}),':', MtU{5}];
                
            case 'sgolay'
                % Check Dataset format
                if ~strcmp(infoDts.Format, 'profile')
                    error('The original dataset should be in profile mode, dataset %i is in %f mode',...
                        dts, infoDts.Format);
                end
                
                % check parameter for method
                MtU{1} = 'SGolay';
                if length(MtU) < 2
                    MtU{2} = 5;
                else
                    MtU{2} = str2double(MtU{2});
                end
                
                if length(MtU) < 3
                    MtU{3} = 1;
                else
                    MtU{3} = str2double(MtU{3});
                end
                
                options.method = ['SGolay:', num2str(MtU{2}),':', ...
                    num2str(MtU{3})];
                
            otherwise
                error('Method not recognised')
        end
        
        % 2- Default parameters
        options.RemSpks = false;
        options.SpkSz   = 1;
        
        % 3- Decipher varargin
        input = @(x) find(strcmpi(varargin,x),1);
        tgtIx = input('spikes');
        if ~isempty(tgtIx)
            spks = varargin{tgtIx +1};
            if spks == 0
                options.RemSpks = false;
            else
                options.RemSpks = true;
                options.SpkSz   =  spks;
            end
        end
    end
end

