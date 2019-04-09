%% Description
% SPIKESREMOVAL removes spikes, that is peaks with spkSz or fewer data
% points for a nx2 dataIn array (1st column axe, 2nd column value).  
%
%% Copyright 
% BSD 3-Clause License
% Copyright 2016-2017 G. Erny (guillaume@fe.up,pt), FEUP, Porto, Portugal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dataOut = spikesRemoval(dataIn, spkSz )

% if ~(spkSz == 1 ||spkSz == 2 ||spkSz == 3)
%     error('spksSz should be 1, 2 or 3')
% end

dataOut = dataIn;

if spkSz >= 10
    findZeros = find(dataOut(:,2) == 0);
    ind2null = findZeros(diff(findZeros) > 10 & diff(findZeros) <= 11);
    dataOut(ind2null+1, 2) = 0;
end

if spkSz >= 9
    findZeros = find(dataOut(:,2) == 0);
    ind2null = findZeros(diff(findZeros) > 9 & diff(findZeros) <= 10);
    dataOut(ind2null+1, 2) = 0;
end

if spkSz >= 8
    findZeros = find(dataOut(:,2) == 0);
    ind2null = findZeros(diff(findZeros) > 8 & diff(findZeros) <= 9);
    dataOut(ind2null+1, 2) = 0;
end

if spkSz >= 7
    findZeros = find(dataOut(:,2) == 0);
    ind2null = findZeros(diff(findZeros) > 7 & diff(findZeros) <= 8);
    dataOut(ind2null+1, 2) = 0;
end

if spkSz >= 6
    findZeros = find(dataOut(:,2) == 0);
    ind2null = findZeros(diff(findZeros) > 6 & diff(findZeros) <= 7);
    dataOut(ind2null+1, 2) = 0;
end

if spkSz >= 5
    findZeros = find(dataOut(:,2) == 0);
    ind2null = findZeros(diff(findZeros) > 5 & diff(findZeros) <= 6);
    dataOut(ind2null+1, 2) = 0;
end

if spkSz >= 4
    findZeros = find(dataOut(:,2) == 0);
    ind2null = findZeros(diff(findZeros) > 4 & diff(findZeros) <= 5);
    dataOut(ind2null+1, 2) = 0;
end

if spkSz >= 3
    findZeros = find(dataOut(:,2) == 0);
    ind2null = findZeros(diff(findZeros) > 3 & diff(findZeros) <= 4);
    dataOut(ind2null+1, 2) = 0;
end

if spkSz >= 2
    findZeros = find(dataOut(:,2) == 0);
    ind2null = findZeros(diff(findZeros) > 2 & diff(findZeros) <= 3);
    dataOut(ind2null+1, 2) = 0;
end

findZeros = find(dataOut(:,2) == 0);
ind2null = findZeros(diff(findZeros) > 1 & diff(findZeros) <=2);
dataOut(ind2null+1, 2) = 0;
end

