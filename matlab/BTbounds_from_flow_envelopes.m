function [ theta, Vdif, stat_iter ] = BTbounds_from_flow_envelopes( bus, branch, opt, Vdif_in, tighten_angle )
%Performs iterative bound tightening and returns bounds on theta, Vdif, or 
%both theta and Vdif for all branches. If squares of voltage magnitudes are
%used, iterative tightening


For Vdif, only one iteration makes sense because they do not affect convex 
%envelopes.
   
%record initial values of output parameters
theta=branch.theta;
Vdif=Vdif_in;
stat_iter=struct('theta',[],'Vdif',[]);

%adjust parameters depending on what type of bounds is tightened
if (opt.


if (tighten_angle)
    max_iter=opt.max_iter;
    eps_tol=opt.max_change;
    bounds=theta;
else
    max_iter=1;
    eps_tol=inf;
    pow=opt.Vdif_type;
    bounds=struct('min',bus.Vmin(branch.ind_bus_F).^pow-bus.Vmax(branch.ind_bus_T).^pow,...
        'max',bus.Vmax(branch.ind_bus_F).^pow-bus.Vmin(branch.ind_bus_T).^pow);
end

%start iterative process of tightening bounds
iter=1;
while iter<=max_iter
    %build all envelopes
    if (tighten_angle)
        Env = BTfinal_flow_envelopes( bus, branch, bounds, opt.Vdif_type-1, Vdif_in );
    else
        Env = BTfinal_flow_envelopes( bus, branch, theta, opt.Vdif_type-1, Vdif_in );
    end
    
    %tighten bounds
    [ bounds_new ] = BTtighten_through_flows( bus, branch, Env, opt, tighten_angle );
    
    %ensure lb<ub
    ind_temp=bounds_new.min>bounds_new.max;
    lb_temp=bounds_new.min(ind_temp);
    bounds_new.min(ind_temp)=bounds_new.max(ind_temp);
    bounds_new.max(ind_temp)=lb_temp;
    
    %compute max changes in bounds compared to previous iteration
    max_bound_change=max(max(bounds.max-bounds_new.max),max(bounds_new.min-bounds.min));
    bounds_temp=bounds;
    
    %record tighter bounds
    bounds=BTupdate_bounds(bounds, bounds_new);
    
    %collect iteration-specific statistics if needed
    if (opt.statistics==2)
        stat_iter=[stat_iter; BTcollect_statistics( bounds_temp, bounds_new, [], 2, 0 )];
    end
    
    %check stopping condition
    if (max_bound_change<eps_tol)
        break;
    else
        iter=iter+1;
    end
end


% postprocess bounds (to store tightest bounds for all parallel branches)
[ bounds.min, bounds.max ] = BTpostprocess_bounds( bounds.min, bounds.max, branch.uniq );
    
end

