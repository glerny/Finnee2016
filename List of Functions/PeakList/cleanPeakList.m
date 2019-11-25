function [myPL, graph] = cleanPeakList(myPL, SN)

PIP2Kill = false(length(myPL.LstPIP{1}), 1);
for ii = 1:length(PIP2Kill)
    Data       = myPL.LstPIP{1}{ii}.x;
    Data(:,2)  = myPL.LstPIP{1}{ii}.y;
    Data(:,3)  = sgolayfilt(Data(:,2), 3, 5);
    cFOM       = myPL.LstPIP{1}{ii}.FOM;
    
    
    D4N = Data(Data(:, 2) ~=0, 2);
    cNoise = 2*1.96*std(Data(Data(:,2) < prctile(D4N, 25), 2));
    
    % Local maxima
    wdz = 4;
    yy = Data(:,3);
    yy = [zeros(wdz, 1); Data(:,3); zeros(wdz, 1)];
    A  = zeros(length(yy), 2*wdz+1);
    for jj = -wdz:wdz
        A(:,  jj + wdz + 1) = circshift(yy, jj);
    end
    [~, Id] = max(A, [], 2);
    IM = find(Id == wdz+1) - wdz;
    IM = IM(find(Data(IM, 2) >= SN*cNoise));
    
    % Local minima
    yy = Data(:,3);
    yy = [inf(wdz, 1); Data(:,3); inf(wdz, 1)];
    A  = zeros(length(yy), 2*wdz+1);
    for jj = -wdz:wdz
        A(:,  jj + wdz + 1) = circshift(yy, jj);
    end
    [~, Id] = min(A, [], 2);
    Im = find(Id == wdz+1) - wdz;
    Im = Im(find(Data(Im, 2) < 0.25*max(Data(:,2))));
    
    if isempty(IM)
        PIP2Kill(ii) = true;
        
    else
        if length(Im) < 2
            PIP2Kill(ii) = true;
        else
            im = find(Im >= IM(end), 1, 'first');
            ip = find(Im <= IM(1), 1, 'last');
            if Im(im) - Im(ip) <= 5
                PIP2Kill(ii) = true;
            end
        end
    end
    
end

graph.TIPin  = myPL.TIP{1}.Data(:,1); graph.TIPin(:,2)  = 0;
graph.TIPout = myPL.TIP{1}.Data(:,1); graph.TIPout(:,2) = 0;
graph.BPPin  = myPL.TIP{1}.Data(:,1); graph.BPPin(:,2)  = 0;
graph.BPPout = myPL.TIP{1}.Data(:,1); graph.BPPout(:,2) = 0;

for ii = 1:length(PIP2Kill)
    cData       = myPL.LstPIP{1}{ii}.x;
    cData(:, 2) = myPL.LstPIP{1}{ii}.y;
    ix = find(ismember(graph.TIPin(:,1), cData(:,1)));
    
    if PIP2Kill(ii)
        graph.TIPout(ix, 2) = graph.TIPout(ix, 2) + cData(:, 2);
        graph.BPPout(ix, 2) = max([graph.BPPout(ix, 2), cData(:, 2)], [], 2);
        
    else
        graph.TIPin(ix, 2)  = graph.TIPin(ix, 2) + cData(:, 2);
        graph.BPPin(ix, 2)  = max([graph.BPPin(ix, 2), cData(:, 2)], [], 2);
        
    end
end

myPL.LstPIP{1}(PIP2Kill) = [];
myPL.FOM{1}.Data(PIP2Kill, :) = [];
myPL.FOM{1}.Data(:,1) = 1:length(myPL.FOM{1}.Data(:,1));


end

