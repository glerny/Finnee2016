classdef Options4bslCor
    
    properties
        Path2Fin     = ''
        Link2Log
        Link2Noise  = ''
        Link2SelPrf = ''
        Link2RemPrf = ''
        TimeAxe
        MzAxe
        FIS
        IndMin
        IndMax
        Precision
        Profiles
        ProfilesTitle = {'Time', 'Noise', 'TIP1', 'BPP1', 'TIP2', 'BPP2'}
        Fmin
        Fmax
        NoiseEstim
        BaselineChoice
    end
    
    methods
        function obj = Options4bslCor(DtsIn)
            if strcmp(DtsIn.Format, 'MS profile')
                
            else
                error('Incorrect Format. The dataset should be made of MS profile scans');
            end
            
            obj.Path2Fin    = DtsIn.Path2Fin;
            obj.Link2Log    = DtsIn.Log;
            obj.TimeAxe     = DtsIn.TimeAxe;
            obj.MzAxe       = DtsIn.MzAxe;
            obj.FIS         = DtsIn.FIS.Data(:,2)/...
                length(DtsIn.TimeAxe.Data)*100;
            obj.Link2SelPrf = tempname;
            obj.Link2RemPrf = tempname;
            obj.Precision   = 'uint32';
            
            figure
            [N,edges] = histcounts(obj.FIS ,1:1:100);
            bar(edges, [0, N]);
            axis([0 100 0 inf])
            
            hPlot = gcf;
            title('Frequency profile - Lower part:noise; higher part : bckg correction');
            
            h = msgbox('Select the low-end value in the frequency profile', 'Noise','modal');
            uiwait(h)
            figure(hPlot);
            [fmin, ~] = ginput(1);
            if isempty(fmin)
                return
            end
            fmin = inputdlg('Confirm the value', 'Noise', 1, {num2str(fmin)});
            
            h = msgbox('Select the high-end value in the frequency profile', 'Correct','modal');
            uiwait(h)
            figure(hPlot);
            [fmax, ~] = ginput(1);
            if isempty(fmax)
                return
            end
            fmax = inputdlg('Confirm the value', 'Correct', 1, {num2str(fmax)});
            
            obj.IndMin = obj.FIS <= str2double(fmin{1});
            obj.IndMax = obj.FIS >= str2double(fmax{1});
            obj.Fmin = str2double(fmin{1});
            obj.Fmax = str2double(fmax{1});
            
            h = waitbar(0, 'Generating matrix of profiles');
            fidWSelPrf = fopen(obj.Link2SelPrf, 'ab');
            fidWRemPrf = fopen(obj.Link2RemPrf , 'ab');
            
            profile(:,1) = DtsIn.TimeAxe.Data;
            data(:,1) = obj.MzAxe.Data;
            for ii = 1:length(profile(:,1))
                waitbar(ii/length(profile(:,1)))
                
                data(:,2) = 0;
                MS = DtsIn.ListOfScans{ii}.Data;
                cst = DtsIn.ListOfScans{ii}.Variables;
                rdg = DtsIn.RndFactor;
                axe = MS(:,1)/(1-cst);
                [tf, loc] = ismember(round(axe*10^rdg)/10^rdg, ...
                    round(data(:,1)*10^rdg)/10^rdg);
                
                indZeros = find(tf == 0);
                if ~isempty(indZeros)
                    [~, locc] = ismember(ceil(axe*10^rdg)/10^rdg, ...
                        ceil(data(:,1)*10^rdg)/10^rdg);
                    loc(indZeros) = locc(indZeros);
                end
                data(loc,2) =  MS(:,2);
                
                
                fwrite(fidWSelPrf, data(obj.IndMax, 2), obj.Precision);
                profile(ii,2) = max( data(obj.IndMin, 2));
                profile(ii,3) = sum( data(obj.IndMax, 2));
                profile(ii,4) = max( data(obj.IndMax, 2));
                fwrite(fidWRemPrf, data(~obj.IndMax, 2), obj.Precision);
                profile(ii,5) = sum(data(~obj.IndMax, 2));
                profile(ii,6) = max( data(~obj.IndMax, 2));
                
                
            end
            close(h)
            fclose(fidWSelPrf);
            fclose(fidWRemPrf);
            obj.Profiles = profile;
            
            figure('Name', 'Selected Profiles for Baseline Corrections')
            subplot(5, 1, 1)
            plot(profile(:,1), profile(:,2));
            obj.NoiseEstim = round(mean(profile(:,2)));
            title(['Estimated noise: ', num2str(obj.NoiseEstim)])
            ylabel([DtsIn.ZLabel, ' / ', DtsIn.ZUnit]);
            
            subplot(5, 1, [2 3])
            plot(profile(:,1), profile(:,4));
            title(['Selected profiles for correction (max plot): ', num2str(sum(obj.IndMax))])
            ylabel([DtsIn.ZLabel, ' / ',  DtsIn.ZUnit]);
            
            subplot(5, 1, [4 5])
            plot(profile(:,1), profile(:,6));
            title(['Non-selected profiles (max plot): ', num2str(sum(~obj.IndMax))])
            xlabel([DtsIn.XLabel, ' / ', DtsIn.XUnit]);
            ylabel([DtsIn.ZLabel, ' / ', DtsIn.ZUnit]);
            
        end
        
        function obj = setBslParameters(obj)
            % Figure and axes
            f = figure('Units', 'normalized', 'Visible', 'off', 'CloseRequestFcn', @done);
            ax1 = axes('Units', 'normalized', 'Position', [0.03 0.6 0.67, 0.344]);
            ax2 = axes('Units', 'normalized', 'Position', [0.03 0.08 0.67, 0.344]);
            edNoise = uicontrol('Style', 'edit',...
                'String', num2str(obj.NoiseEstim),...
                'Units', 'normalized',...
                'Callback', @changeEdit ,...
                'HorizontalAlignment', 'left',...
                'FontWeight', 'bold',...
                'Position', [0.6 0.38 0.02 0.035]);
            
            uicontrol('Style', 'text',...
                'String', 'Noise: ',...
                'Units', 'normalized',...
                'HorizontalAlignment', 'left',...
                'FontWeight', 'bold',...
                'Position', [0.57 0.38 0.02 0.035]);
            
            edWdz = uicontrol('Style', 'edit',...
                'String', '3',...
                'Units', 'normalized',...
                'Callback', @changeEdit ,...
                'HorizontalAlignment', 'left',...
                'FontWeight', 'bold',...
                'Position', [0.6 0.35 0.02 0.035]);
            
            uicontrol('Style', 'text',...
                'String', 'WDZ: ',...
                'Units', 'normalized',...
                'HorizontalAlignment', 'left',...
                'FontWeight', 'bold',...
                'Position', [0.57 0.35 0.02 0.035]);
            
            % List box
            listMethod = uicontrol('Style', 'listbox', ...
                'String', {'arPLS'},...
                'Value', 1,...
                'Units', 'normalized',...
                'Tag', 'listMethod',...
                'Enable', 'off',...
                'Position', [0.8 0.7 0.16 0.25]);
            
            % Edit box
            editPara1 = uicontrol('Style', 'edit', ...
                'String', '1E6',...
                'Units', 'normalized',...
                'Callback', @changeEdit ,...
                'Tag', 'editPara1',...
                'Position', [0.8 0.63 0.05 0.03]);
            
            uicontrol('Style', 'text',...
                'String', 'smoothness (0 - 1E15)',...
                'Units', 'normalized',...
                'HorizontalAlignment', 'left',...
                'FontWeight', 'bold',...
                'Position', [0.86 0.63, 0.1, 0.03]);
            
            % Push buttions
            
            uicontrol('Style', 'pushbutton', ...
                'String', 'Next profile',...
                'Units', 'normalized',...
                'Callback', @newProfile,...
                'Position', [0.79 0.37 0.06 0.044]);
            
            uicontrol('Style', 'pushbutton', ...
                'String', 'Done',...
                'Units', 'normalized',...
                'Callback', @done,...
                'Position', [0.89 0.37 0.06 0.044]);
            
            yy = getProfile(obj);
            parametersBaseline = doBaselineCorrection(yy);
            
            function newProfile(~, ~)
                yy = getProfile(obj);
                parametersBaseline = doBaselineCorrection(yy);
            end
            
            function done(~, ~)
                parametersBaseline.obj = obj;
                parametersBaseline.noise =  str2double(edNoise.String);
                parametersBaseline.wdz = str2double(edWdz.String);
                assignin('base', 'par4bas', parametersBaseline)
                delete(gcf)
            end
            
            function parametersBaseline = changeEdit(~, ~)
                parametersBaseline = doBaselineCorrection(yy);
            end
            
            function currentProfile = getProfile(obj)
                nbrProf = sum(obj.IndMax);
                lengthTm = length(obj.TimeAxe.Data);
                tgtProf = round(rand*(nbrProf-1))+1;
                fidRead = fopen(obj.Link2SelPrf, 'rb');
                fseek(fidRead, (tgtProf-1)*4, 'bof');
                currentProfile(:,1) = obj.TimeAxe.Data;
                currentProfile(:,2) = ...
                    fread(fidRead, lengthTm, obj.Precision, (nbrProf-1)*4);
                fclose(fidRead);
            end
            
            function parametersBaseline = doBaselineCorrection(yy)
                value = get(listMethod, 'Value');
                switch value
                    case 1
                        parameter(1) = str2double(get(editPara1, 'String'));
                        notZeros = find(yy(:,2) ~= 0);
                        [baseline, ~, parameter] = doArPLS(yy(notZeros,2),...
                            parameter(1));
                        set(editPara1, 'String', parameter(1));
                        parametersBaseline.type = 'arPLS';
                        parametersBaseline.parameter(1) = str2double(get(editPara1, 'String'));
                        
                end
                
                
                axes(ax1)
                plot(yy(:,2), 'k'), hold on, plot(notZeros, baseline, 'r'), hold off;
                
                bsl = zeros(length(yy(:,2)), 1);
                bsl(notZeros) = round(baseline);
                filteredyy = yy(:,2) - bsl;
                assignin('base', 'ftyy', filteredyy)
                noise = str2double(edNoise.String);
                wdz = str2double(edWdz.String);
                mat4noise = zeros(length(filteredyy), 2*wdz+1);
                for ii = 1:2*wdz+1
                    id1 = max(wdz+2-ii, 1);
                    id2 = max(ii-wdz, 1);
                    mat4noise(id1:end-id2+1, ii) = filteredyy(id2:end-id1+1);
                end
                w = max(mat4noise, [], 2) < 3*noise;
                filteredyy(w) = 0;
                indNeg = filteredyy <= 0;
                filteredyy(indNeg) = 0;
                
                axes(ax2)
                plot(filteredyy, 'k');
                
            end
        


set(f, 'Visible', 'on');
        end
        
    end
    
end

