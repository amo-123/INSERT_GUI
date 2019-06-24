function [ result_fit,Delta] = FcnGaussFitV2(x,y_flt,x_max,y_max,Soglia_fit_dx,Soglia_fit_sx )
%FCNGAUSSFIT Summary of this function goes here
%   Detailed explanation goes here
    indice=0;
    while y_flt(dsearchn(x,x_max)+indice)>Soglia_fit_sx*y_max
        indice=indice+1;
    end
    Delta_sup=dsearchn(x,x_max)+indice;

    indice=0;
    while y_flt(dsearchn(x,x_max)-indice)>Soglia_fit_dx*y_max
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

