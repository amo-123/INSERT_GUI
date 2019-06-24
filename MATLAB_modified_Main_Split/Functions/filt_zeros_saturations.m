function [Frame_others]=filt_zeros_saturations(Frame,zero_thr,sat_thr,op_mode)

bins=100;

%find Frame portions corresponding to zeros and saturations
[under_zero_thr, x_under_zero_thr]=hist(sum(Frame(min(Frame,[],2)<=zero_thr,:),2),bins);
[over_sat_thr, x_over_sat_thr]=hist(sum(Frame(max(Frame,[],2)>=sat_thr,:),2),bins);
Frame_others=Frame(and(min(Frame,[],2)>zero_thr,max(Frame,[],2)<sat_thr),:);
[others,x_others]=hist(sum(Frame_others,2),bins);
[all,x_all]=hist(sum(Frame,2),bins);
tot=sum(all);

if strcmp(op_mode,'plot')
    plot(x_all,all,'k',x_others,others,'b',x_under_zero_thr,under_zero_thr,'c',x_over_sat_thr,over_sat_thr,'r','Linewidth',2);
    xlabel('Ssum')
    ylabel('counts')
    legend('original','filtered','zeros','saturations')
elseif strcmp(op_mode,'filt')
    disp([['zeros = ',num2str(round(sum(under_zero_thr)/tot*100,2)),' %'], 32, ['saurations = ',num2str(round(sum(over_sat_thr)/tot*100,2)),' %'], 32, ['num events survived = ',num2str(size(Frame_others,1))]])
end

end

