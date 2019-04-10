function [ val ] = comdat( sg, tag, data )
%Set or get data in common structure CB

    persistent CB

    if ( nargin == 0 )
        if ( isempty(CB) ), val=0; else val=length(fieldnames(CB)); end
        return
    end

    val = [];

    if strcmp( sg, 'set' )

        if     ( nargin == 1 ), CB = [];
        elseif ( nargin == 2 ), CB.(tag) = [];
        elseif ( nargin == 3 ), CB.(tag) = data; 
        end
        
    elseif strcmp( sg, 'get' )

        if ( isempty(CB) ), return, end

        if ( nargin == 1 )
            val = fieldnames(CB); 
        elseif isfield( CB, tag )
            val = CB.(tag); 
        end
        
    end
    
end
