function [] = vline_spectrum_plot(bins,E_min,E_max)

hold on
v_plot=1:max(bins);
v_resh_min=E_min.*ones(length(v_plot),1);
v_resh_max=E_max.*ones(length(v_plot),1);
plot(v_resh_min,v_plot,'--c',v_resh_max,v_plot,'--c')
hold off

end

