function [offsets] = fcn_Extract_Offsets(coefficients)
 %% La funzione permette di estrarre le intercette delle rette di calibrazione ottenute dall'impulsazione dei canali
    
    names = fieldnames(coefficients);
    a = names{1};
    b = names{2};
    n_offs = numel(names);
    match = zeros(1, min(numel(a), numel(b)));
    
    for i = 1:min(numel(a), numel(b))
        match(i) = strcmp (a(i), b(i));
    end
    match = logical(match);
    word = a(match);
    found_ch = zeros(1, n_offs);
    
    for i = 1:n_offs
        index = strfind(names{i}, word) + numel(word);
        found_ch(i) = str2double( names{i}(index:end));
    end
    
    offsets = zeros(n_offs, 1);
    
    for i = 1:n_offs
        offsets(found_ch(i)) = coefficients.(strcat(word,num2str(found_ch(i)))).q;
    end

    offsets = round(offsets);
    
    selection = menu('SAVE OFFSET .mat VARIABLE:','Yes','No');
    if selection == 1
        uisave('offsets')
    end

end

