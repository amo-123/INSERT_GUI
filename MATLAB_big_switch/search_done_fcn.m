%% ricerca funzione precedente


function [ output,ws ] = search_done_fcn( string_fcn_type, output,ws, do_it)

%% ricerca
num = numel(output);
for i=num:-1:0
    if(i>0)
        if (strcmp(output{i}.fcn_type, string_fcn_type))
            ws.sequence_num = i;
            break
        end
    end
    %% non trovato
    if (i==0 && do_it==1)
        %% conversione
        selected_fnc = conversione_fcnType2num(string_fcn_type);
        if isnan(selected_fnc)
            %% non esiste un campo così nominato
            disp(['ERROR: fnc_type = "', string_fcn_type, '" can not be found'])
            ws.sequence_num = NaN;
        else
            %% lo eseguo
            disp(['WARNING: fnc_type = "', string_fcn_type, '" not executed yet'])
            disp(['WARNING: it will be executed with sequence number = ', num2str(selected_fnc)])
            big_switch
        end
    end
    if (i==0 && do_it==0)
        ws.sequence_num = NaN;
    end
end

end

