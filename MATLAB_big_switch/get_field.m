%% search in 'output' by fieldname and return the value of the searched item

function [ value, output,ws ] = get_field( searched_string, string_fcn_type, output,ws, varargin )

idx_in = ws.idx_fnc;
% string_fcn_type = 'fragment in Frame_node';
% searched_string = 'Frame_node.node_1';  

%% indico nome dell'operazione svolta dal nodo
[ output,ws ] = search_done_fcn( string_fcn_type, output,ws, 1);
ws.idx_fnc = numel(output)+1;
[ output,ws ] = search_done_fcn( string_fcn_type, output,ws, 0);

%% indico nome della struttura che sto cercando
if ( any(strcmp(fieldnamesr(output{ws.sequence_num}), searched_string)) )
    if isempty(varargin)
        value = eval(['output{',num2str(ws.sequence_num),'}.',searched_string]);
    elseif strcmp(varargin{1},'ref')
        value = ['output{',num2str(ws.sequence_num),'}.',searched_string];
    end
else
    disp(['ERROR: fieldname ', searched_string, ' not found'])
end

ws.idx_fnc = idx_in;
end