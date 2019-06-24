function [result_fit,Delta]=FcnGaussFitV2(x,y_flt,x_max,y_max,soglia_fit_dx,soglia_fit_sx )

if (isrow(x))
    x = x';
end

if(isrow(y_flt))
    y_flt=y_flt';
end

    indice=0;
    while y_flt(dsearchn(x,x_max)+indice)>soglia_fit_sx*y_max
        indice=indice+1;
    end
    Delta_sup=dsearchn(x,x_max)+indice;

    indice=0;
    while y_flt(dsearchn(x,x_max)-indice)>soglia_fit_dx*y_max
        indice=indice+1;
    end
    Delta_inf=dsearchn(x,x_max)-indice;
    
    fitOptions = fitoptions('gauss1');
    fitOptions.Lower = zeros(1,3);
    fitOptions.Upper = [+Inf +Inf +Inf];

    result_fit=fit(x(Delta_inf:Delta_sup), y_flt(Delta_inf:Delta_sup),'gauss1',fitOptions);
    Delta.inf=Delta_inf;
    Delta.sup=Delta_sup;
    
end

