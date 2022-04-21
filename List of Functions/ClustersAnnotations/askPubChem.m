function [OutputList, Count] = askPubChem(QueryName, QueryValue, options)


eutil = 'http://www.ncbi.nlm.nih.gov/entrez/eutils';
option = weboptions('Timeout',10);
PropIsDouble = {'MolecularWeight', 'MonoisotopicMass', 'XLogP'};

switch lower(QueryName)
    case 'monoisotopicmass'
        if ~isnumeric(QueryValue)
            error('QueryValue shoud be a mass interval: [mi_min mi_max]');
        end
        
        if min(size(QueryValue)) ~= 1 || max(size(QueryValue))~= 2
            error('QueryValue shoud be a mass interval: [mi_min mi_max]');
        end
        
        if QueryValue(1) >= QueryValue(2)
            error('QueryValue shoud be a mass interval: [mi_min mi_max]');
        end
        
        query = sprintf('%.6f:%.6f[MonoisotopicMass]', QueryValue); 
end

esearch = ['/esearch.fcgi?db=pccompound&term=' query '&retmax=0&usehistory=y'];
output  = webread([eutil esearch]);
is = strfind(output, '<Count>');
ie = strfind(output, '</Count>');
Count = str2double(output(is(1)+7:ie(1)-1));
if Count == 0
    OutputList = table();
    return
end

is = strfind(output, '<QueryKey>');
ie = strfind(output, '</QueryKey>');
QueryKey = output(is+10:ie-1);

is = strfind(output, '<WebEnv>');
ie = strfind(output, '</WebEnv>');
WebEnv = output(is+8:ie-1);

cgi = 'https://pubchem.ncbi.nlm.nih.gov/list_gateway/list_gateway.cgi?';
action = sprintf('action=entrez_to_pug&entrez_db=pccompound&entrez_query_key=%s&entrez_webenv=%s', QueryKey, WebEnv);
url = [cgi action];
output  = webread(url);

is = strfind(output, '<Response_pug-listkey>');
ie = strfind(output, '</Response_pug-listkey>');
listkey = output(is+22:ie-1);
pugrest    = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/';
input      = sprintf('compound/listkey/%s/', listkey);
operation  = 'Property/MolecularFormula,MolecularWeight,MonoisotopicMass,XLogP,CanonicalSMILES,InChIKey,IUPACName/';
output     = 'XML';

pug_url    = [pugrest input operation output];

while 1
    try
        OutputList     = webread(pug_url, option);
        break
        
    catch EM
        if ~strcmp(EM.identifier, 'MATLAB:webservices:Timeout')
            error('EM')
        else
            pause(30)
            fprintf('\nMATLAB:webservices:Timeout')
            fprintf('\nIncreasing Timeout (ctrl^C to stop)\n')
            option = weboptions('Timeout',60);
        end
    end
end

% Error as time with CSV as output need to decifer XML
%% GET FIELD FROM operation %%
Fields = strsplit(operation, '/');
BlckField = strsplit(Fields{2}, ',');

%% Filling the table
List  = strsplit(OutputList, '<Properties>');
for ii = 2:size(List, 2)
    for jj = 1:length(BlckField)
        is = strfind(List{ii}, ['<' BlckField{jj} '>']);
        ie = strfind(List{ii}, ['</' BlckField{jj} '>']);
        
        if any(is) && any(ie)
            Valor = List{ii}(is+length(BlckField{jj})+2:ie-1);
            if any(strcmp(PropIsDouble, BlckField{jj}))
                Valor = str2double(Valor);
            end
            newOutput.(BlckField{jj}){ii-1, 1} = Valor;
        else
            if any(strcmp(PropIsDouble, BlckField{jj}))
                newOutput.(BlckField{jj}){ii-1, 1} = NaN;
            else
                newOutput.(BlckField{jj}){ii-1, 1} = '';
            end
            
        end
    end
    
end
OutputList = struct2table(newOutput);
end