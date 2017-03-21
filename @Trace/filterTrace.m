%% DESCRIPTION
%
%% Copyright
% BSD 3-Clause License
% Copyright 2016-2017 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function obj = filterTrace(obj, method, varargin)

narginchk(1, inf)
infoTrc       = obj.InfoTrc;
if nargin == 1, method = 'RemoveSpikes:2'; end
[options, MtU] = checkVarargin(infoTrc, method, varargin{:});
infoTrc.Loc   = 'inTrace';
infoTrc.Title = ['SPR - ', infoTrc.Title];

switch lower(MtU{1})
    case 'removespikes'
        obj = Trace(infoTrc, spikesRemoval(obj.Data, MtU{2}));
    otherwise
        error('%s, is not a recognised method', MtU{1})
end

    function [options, MtU] = checkVarargin(info, method, varargin)
        %%% Check method
        MtU = strsplit(method, ':');
        switch lower(MtU{1})
            case 'removespikes'
                MtU{1} = 'RemoveSpikes';
                if length(MtU) == 1
                    MtU{2} = 2;
                else
                    MtU{2} = str2double(MtU{2});
                end
                options.method = ['RemoveSpikes:', num2str(MtU{2})];
            otherwise
                error('Method not recognised')
        end
        
        if ~strcmp(info.TT, 'PRF')
            error('RemoveSpikes only works with MS scans recorded in profile mode');
        end
        
        %         %%% Decipher varargin
        %         input = @(x) find(strcmpi(varargin,x),1);
        %         tgtIx = input('tLim');
        %         if ~isempty(tgtIx)
        %             tLim         = varargin{tgtIx +1};
        %             options.XLim = [min(tLim) max(tLim)];
        %         end
    end
end

