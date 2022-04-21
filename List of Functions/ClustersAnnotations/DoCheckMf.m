function [mf, mi, a, Z] = DoCheckMf(mfString, Data)
MassElec = 0.000548579909;

a = Data;
a.nbElements = zeros(size(a, 1), 1);
a = sortrows(a, 'Priority', 'descend');
Z = 0;

for ii = 1:size(a, 1)
    if any(strfind(mfString, a.Symbol{ii}))
        
        IList = strfind(mfString, a.Symbol{ii});
        for jj = length(IList):-1:1
            is = IList(jj);
            
            if is + a.Letters(ii) > length(mfString)
                ie = length(mfString);
                nbr = 1;
                
            elseif isnan(str2double(mfString(IList(jj)+a.Letters(ii)))) &&...
                    mfString(IList(jj)+a.Letters(ii)) ~= '-'
                ie = IList(jj)+a.Letters(ii)-1;
                nbr = 1;
                
            else
                ie  = IList(jj)+a.Letters(ii);
                while 1
                    if ie+1 > length(mfString), break; end
                    if isnan(str2double(mfString(ie+1))), break; end
                    ie = ie + 1;
                end
                nbr = str2double(mfString(is+a.Letters(ii):ie));
            end
            
            mfString(is:ie) = [];
            a.nbElements(ii) = a.nbElements(ii) + nbr; 
        end
    end
end

if ~isempty(mfString)
    is = strfind(mfString, '(')+1;
    if length(is) ~= 1, error(''); end
    ie = strfind(mfString, ')')-1;
    if length(ie) ~= 1, error(''); end
    Z = str2double(mfString(is:ie));
end

a = sortrows(a, 'Z', 'ascend');
mi =  sum(a.MonoisotopicMass.*a.nbElements);
mf =[];

I2mf = find(a.nbElements ~= 0);
for ii = 1:length(I2mf)
    mf = [mf a.Symbol{I2mf(ii)} num2str(a.nbElements(I2mf(ii)))];
end

mi = mi -Z*MassElec;
if Z > 0
    mf = [mf '(+' num2str(Z) ')'];
elseif Z < 0 
    mf = [mf '(' num2str(Z) ')'];
end

