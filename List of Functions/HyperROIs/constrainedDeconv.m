function [Unconstr, dec_Data]  = constrainedDeconv(Data, Model, mode)

dec_Data = zeros(size(Data, 1), size(Model, 1));
Unconstr = zeros(size(Data, 1), size(Model, 1));
opts = optimset('Display','off');
mCP = 5;

% TODO: Check mode

if ~strcmp(mode, 'full'), mode = 'fast'; end
for ii = 1:size(dec_Data, 1)
    
    if minContinuousPoints(Data(ii, :)', 1, 0) > mCP
        Unconstr(ii, :) = Data(ii, :)/Model;
        
        if strcmp(mode, 'full')
            x0 = Unconstr(ii, :);
            x0(x0 < 0 ) = 0;
            x0(isnan(x0)) = 0;
            
            [x, ~, exitflag] = fminsearch(@myfun, x0, opts);
            
            count = 1;
            while exitflag == 0
                [x, ~, exitflag] = fminsearch(@myfun, x, opts);
                count = count +1;
                if count > 10
                    break
                end
            end
            
            dec_Data(ii, :) = x';
        end
    end
end

    function [f, cModel] = myfun(x)
        
            cModel   =  zeros(size(Model));
        for jj = 1:size(Model, 1)
            if x(jj) < 0
                cModel(jj, :) = inf(1, size(Model, 2));
            else
                cModel(jj, :) = x(jj)*Model(jj, :);
            end
        end
        
        f = sqrt(sum((Data(ii, :) - sum(cModel, 1)).^2));
        if ~isfinite(f)
            f = 5*sqrt(sum((Data(ii, :)).^2));
        end
    end
end

