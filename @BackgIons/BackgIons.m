%% DESCRIPTION
% BackgIons is the class that deals with three-dimensional representations.
% Those can a cut from the dataset after mzAxis alignment to the Matser mz Axis
% others.
%
%% LIST OF THE CLASS'S PROPERTIES
%
%% LIST OF THE CLASS'S METHODS
%
%% Copyright
% BSD 3-Clause License
% Copyright 2019-2020 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef BackgIons
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ROIs         % list of ROI
        ListIons     % BackgroundIons
        path
        n = 1;
        P = {};
        name = 'BckIons.mat'
    end
    
    methods
        function obj = BackgIons(dts, polarity, WdSz)
            
            thres = 10;
            %1. Claculated the ROIs
            switch polarity
                case {'+', 'pos', 'Pos', 'POS'}
                    str = which('BckgIonsPos.mat');
                    if ~isempty(str)
                        BckgIons = load(str, 'BckgIonsPos');
                        BckgIons = BckgIons.BckgIonsPos;
                    else
                        error('BckgIonsPos.mat not found on the path')
                    end
                    
                case {'-', 'neg', 'Neg', 'NEG'}
                    str = which('BckgIonsNeg.mat');
                    if ~isempty(str)
                        BckgIons = load(str, 'BckgIonsNeg');
                        BckgIons = BckgIons.BckgIonsNeg;
                    else
                        error('BckgIonsNeg.mat not found on the path')
                    end
            end
            ROIs = dts.mkMnROI(BckgIons.MonoMass, WdSz, zeros(size(BckgIons.MonoMass)), inf, obj.name);
            
            %2. Validate ROIs as background ions
            SelBckg = {};
            Summary = [];
            ListIons = table;
            for ii = 1:length(ROIs.ROI)
                if nnz(sum(ROIs.ROI{ii}.StoredData, 1))/length(sum(ROIs.ROI{ii}.StoredData, 1))*100 >= 10
                    XY = [ROIs.ROI{ii}.AxisMZ.Data, mean(ROIs.ROI{ii}.StoredData, 2)];
                    Res = LocalMaxima(XY, 3, 0);
                    if size(Res, 1) == 1
                        SelBckg{end+1} = ROIs.ROI{ii};
                        Summary(end+1, 1) = SelBckg{end}.TgtMz;
                        Summary(end,   2) = Res(1);
                        ListIons(end+1, :) = BckgIons(ii, :);
                    end
                end
            end
            
            Summary(:,3) = (Summary(:,2) - Summary(:,1))./Summary(:,2)*1000000;
            TF = isoutlier(Summary(:,3));
            while any(TF)
                Summary(TF,:)  = [];
                SelBckg(TF) = [];
                ListIons(TF, :) = [];
                TF = isoutlier(Summary(:,3));
            end
            
            obj.ROIs = SelBckg;
            obj.ListIons = ListIons;
            obj.path = fullfile(dts.Path2Fin, obj.name);
            myBackIons = obj;
            save(obj.path, 'myBackIons');
        end
        
    end
end

