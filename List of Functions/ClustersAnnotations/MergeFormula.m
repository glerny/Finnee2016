function rmf = MergeFormula(rMotif, radduct)

rmf = {};
elements = unique([fieldnames(rMotif); fieldnames(radduct)]);

for ii = 1:length(elements)
    nbr = 0;
    
    if isfield(rMotif, elements{ii})
        nbr = nbr + rMotif.(elements{ii});
    end
    
    if isfield(radduct, elements{ii})
        nbr = nbr + radduct.(elements{ii});
    end
    
    if nbr < 0 && ~strcmp(elements{ii}, 'Q')
        rmf = {};
        return
    end
    
    rmf.(elements{ii}) = nbr;
    
end

end

