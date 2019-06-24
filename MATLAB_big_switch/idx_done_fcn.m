function [ idx ] = idx_done_fcn( string_fcn_type, output )

%% ricerca
for i=numel(output):-1:1
    if (strcmp(output{i}.fcn_type, string_fcn_type))
        idx = i;
        i=0;
    end
end

end

