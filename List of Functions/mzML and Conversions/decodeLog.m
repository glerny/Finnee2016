function strLog = decodeLog(Log)

A = strsplit(Log, ' ');
strLog.definition = A{1};

if length(A) > 1
    for ii = 2:length(A)
        B = strsplit(A{ii}, '=');
        strLog.(B{1}) = B{2};
    end
end
    

end

