function [ output,ws ] = do_fcn_type( string_fcn_type, output,ws )

selected_fnc = conversione_fcnType2num(string_fcn_type);
if isnan(selected_fnc)
    disp(['ERROR: fnc_type = "', string_fcn_type, '" can not be found'])
else
    big_switch
end

end

