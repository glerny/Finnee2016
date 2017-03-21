%% DESCRIPTION
%
%% Copyright 
% BSD 3-Clause License
% Copyright 2016-2017 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function [prmOut, decInt, model, SSE] = doPeakMinimization(x, X, peakFct, prmIn)

if size(X,1) > size(X,2)
    X = X';
end

a=[];
switch peakFct
    case 'Gauss'
        for ii = 1:length(prmIn(:,1))
            a(end+1) = prmIn(ii, 1);
            a(end+1) = prmIn(ii, 2);
        end
        
        [a,FVAL,EXITFLAG] = fminsearch(@(a) gaussModel(a, x, X), a);
        [SSE, model, decInt] = gaussModel(a, x, X);
        X_ = (decInt*model);
        
        for ii = 1:2:length(a)
            prmOut((ii-1)/2 + 1, 1) = a(ii);
            prmOut((ii-1)/2 + 1, 2) = a(ii+1);
        end
        
    case 'EMG'
        for ii = 1:length(prmIn(:,1))
            a(end+1) = prmIn(ii, 1);
            a(end+1) = prmIn(ii, 2);
            a(end+1) = -1*a(end)/2; %distortion fixed for optimisation
        end
        
        % Test with negative a3
        [an,FVAL,EXITFLAG] = fminsearch(@(a) EMGModel(a, x, X), a);
        
        while EXITFLAG ~= 1
            [an,FVAL,EXITFLAG] = fminsearch(@(a) EMGModel(a, x, X), an);
        end
        
        [SSEn, modeln, decIntn] = EMGModel(an, x, X);
        % Test with positive a3 (all peaks shoudl have the same distortion
        % type!)
        
        for ii = 3:3:length(a)
            a(ii) = 1*a(ii-1)/2;
        end
        
        [ap,FVAL,EXITFLAG] = fminsearch(@(a) EMGModel(a, x, X), a);  
        
        while EXITFLAG ~= 1
            [ap,FVAL,EXITFLAG] = fminsearch(@(a) EMGModel(a, x, X), ap);
        end
        [SSEp, modelp, decIntp] = EMGModel(ap, x, X);
        
        if SSEn <= SSEp
            a = an; SSE = SSEn; model = modeln; decInt = decIntn;
        else
            a = ap; SSE = SSEp; model = modelp; decInt = decIntp;
        end
        
        X_ = (decInt*model);
        
        for ii = 1:3:length(a)
            prmOut((ii-1)/3 + 1, 1) = a(ii);
            prmOut((ii-1)/3 + 1, 2) = a(ii+1);
            prmOut((ii-1)/3 + 1, 3) = a(ii+2);
        end
        
        
    case 'PMG2'
        for ii = 1:length(prmIn(:,1))
            a(end+1) = prmIn(ii, 1);
            a(end+1) = prmIn(ii, 2);
            a(end+1) = 0;
            a(end+1) = 0;
        end
        
        [a,FVAL,EXITFLAG] = fminsearch(@(a) PMG2Model(a, x, X), a);
        
         while EXITFLAG ~= 1
             [a,FVAL,EXITFLAG] = fminsearch(@(a) PMG2Model(a, x, X), a);
         end
        
        [SSE, model, decInt] = PMG2Model(a, x, X);
        X_ = (decInt*model);
        
        for ii = 1:4:length(a)
            prmOut((ii-1)/4 + 1, 1) = a(ii);
            prmOut((ii-1)/4 + 1, 2) = a(ii+1);
            prmOut((ii-1)/4 + 1, 3) = a(ii+2);
            prmOut((ii-1)/4 + 1, 4) = a(ii+3);
        end
end



    function [SSE, model, dec] = gaussModel(a, x, X)
        nbrPeaks =  length(a)/2;
        model = zeros( nbrPeaks, length(x));
        for jj = 1:nbrPeaks
            model(jj,:) = gaussPeak(x, a(2*jj-1:2*jj));
            model(jj,:) = model(jj,:)/max(model(jj,:));
        end
        
        dec = X/model;
        SSE = sum(sum((X - (dec*model)).^2));
    end

    function [SSE, model, dec] = EMGModel(a, x, X)
        nbrPeaks =  length(a)/3;
        model = zeros( nbrPeaks, length(x));
        for jj = 1:nbrPeaks
            model(jj,:) = EMGPeak(x, a([3*jj-2 2 3*jj]));
            model(jj,:) = model(jj,:)/max(model(jj,:));
        end
        
        dec = X/model;
        SSE = sum(sum((X - (dec*model)).^2));
    end

    function [SSE, model, decInt] = PMG2Model(a, x, X)
        nbrPeaks =  length(a)/4;
        model = zeros(nbrPeaks, length(x));
        for jj = 1:nbrPeaks
            model(jj,:) = PMG2Peak(x, a([4*jj-3 2 4*jj-1 4*jj]));
            model(jj,:) = model(jj,:)/max(model(jj,:));
        end

        decInt = X/model;
        SSE = sum(sum((X - (decInt*model)).^2));
    end

end

