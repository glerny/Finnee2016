%% DESCRIPTION
% PIP is the class that is contain all information related to a single Pure
% Ion Profile
% 
%% LIST OF THE CLASS'S PROPERTIES
% 
%% LIST OF THE CLASS'S METHODS
% 
%% Copyright
% BSD 3-Clause License
% Copyright 2016-2017 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef PIP
    
    properties
        Data    % 3 columns matrix, m/z|Int|Time
        x       % x axis for this PIP
        IdS
        AxisX   
        AxisY
        AxisZ
        Lbl4FOM
    end
    
    properties (Hidden)
        deconvolved = false
    end
    
    properties (Dependent)
        y       % Intensity for this PIP
        FOM     % Figures of merit
                % max I|time @ max I|M0|M1|M2|M3|mean(m/z)|std(m/z)|AccMass
    end
    
    methods
        function obj  = PIP(DataIn, axe, infoDts) % Creator 
            obj.Data    = DataIn;
            obj.x       = axe(min(DataIn(:,3)):max(DataIn(:,3)));
            obj.IdS     = min(DataIn(:,3));
            obj.AxisX   = infoDts.AxisX;
            obj.AxisY   = infoDts.AxisY;
            obj.AxisZ   = infoDts.AxisZ;
            obj.Lbl4FOM = {'IntMax', 'Tm@IntMax', 'M0', 'M1', 'M2', 'M3', ...
                'MeanMZ', 'stdevMZ', 'AccuMass', 'DeltaInt'} ;
        end
        
        function y   = get.y(obj)
            y = zeros(size(obj.x,1), 1);
            C = unique(obj.Data(:,3));
            
            if length(C) ~= size(obj.Data, 1)
                for jj = 1:length(C)
                    y(C(jj) - min(C) + 1) = ...
                        sum(obj.Data(obj.Data(:,3) == C(jj),2));
                end
            else
                y(obj.Data(:,3) - min(obj.Data(:,3))+1) = obj.Data(:,2);
            end
        end
        
        function FOM = get.FOM(obj)
            [FOM(1), Id2x] = max(obj.y);
            FOM(2) = obj.x(Id2x);
            FOM(3) = trapz(obj.x, obj.y);
            FOM(4) = trapz(obj.x, obj.x.*obj.y)/FOM(3);
            FOM(5) = trapz(obj.x, (obj.x-FOM(4)).^2.*obj.y)/FOM(3);
            FOM(6) = trapz(obj.x, (obj.x-FOM(4)).^3.*obj.y)/FOM(3);
            FOM(7) = mean(obj.Data(:,1));
            FOM(8) = std(obj.Data(:,1));
            FOM(9) = sum(obj.Data(:,1).*obj.Data(:,2))/sum(obj.Data(:,2));
            FOM(10) = FOM(1) - min(obj.y);
        end
        
    end
    
end

