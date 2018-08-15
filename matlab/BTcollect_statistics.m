function [ stat, viol ] = BTcollect_statistics( bounds_old, bounds_new, current_val, statistics, check_viol  )
%collect statistical information about the tightening.

%compute some intermediate values
nb=numel(bounds_old.min);
ind_tight=bounds_new.min>bounds_old.min | bounds_new.max<bounds_old.max;
range_old=bounds_old.max-bounds_old.min;
range_new=bounds_new.max-bounds_new.min;
drange=(range_old-range_new)./range_old;

%collect general statistics
stat.nb_p=sum(ind_tight)*100/nb;
stat.range_max=max(range_new);
stat.range_mean=mean(range_new);
stat.range_median=median(range_new);
stat.range_min=min(range_new);
if (stat.nb_p>0)
    stat.drange_max_p=100*max(drange(ind_tight));
    stat.drange_mean_p=100*mean(drange(ind_tight));
else
    stat.drange_max_p=0;
    stat.drange_mean_p=0;
end

%collect more detailed statistics if need be
if (statistics==2)
    ind_tight=bounds_new.min>bounds_old.min;
    stat.n_lb_p=sum(ind_tight)*100/nb;
    if (stat.n_lb_p>0)
        stat.d_lb_max=max(bounds_new.min(ind_tight)-bounds_old.min(ind_tight));
        stat.d_lb_mean=mean(bounds_new.min(ind_tight)-bounds_old.min(ind_tight));
    else
        stat.d_lb_max=0;
        stat.d_lb_mean=0;
    end
    ind_tight=bounds_new.max<bounds_old.max;
    stat.n_ub_p=sum(ind_tight)*100/nb;
    if (stat.n_ub_p>0)
        stat.d_ub_max=max(bounds_old.max(ind_tight)-bounds_new.max(ind_tight));
        stat.d_ub_mean=mean(bounds_old.max(ind_tight)-bounds_new.max(ind_tight));
    else
        stat.d_ub_max=0;
        stat.d_ub_mean=0;
    end
end

%record violations if need be
if (check_viol)
    ind_lb=current_val<bounds_new.min;
    ind_ub=current_val>bounds_new.max;
    viol.n_lb=sum(ind_lb);
    viol.n_ub=sum(ind_ub);
    if (viol.n_lb>0)
        viol.lb_max=max(bounds_new.min(ind_lb)-current_val(ind_lb));
        viol.lb_mean=mean(bounds_new.min(ind_lb)-current_val(ind_lb));
    else
        viol.lb_max=0;
        viol.lb_mean=0;
    end
    if (viol.n_ub>0)
        viol.ub_max=max(current_val(ind_ub)-bounds_new.max(ind_ub));
        viol.ub_mean=mean(current_val(ind_ub)-bounds_new.max(ind_ub));
    else
        viol.ub_max=0;
        viol.ub_mean=0;
    end
else
    viol=[];
end
   
end
