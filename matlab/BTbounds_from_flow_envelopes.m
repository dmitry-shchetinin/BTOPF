function [ theta, Vdif, stat_iter ] = BTbounds_from_flow_envelopes( bus, branch, opt, Vdif, tighten_Vdif_only )
%Performs iterative bound tightening and returns bounds on theta, Vdif, or 
%both theta and Vdif for all branches. If squares of voltage magnitudes are
%used, iterative tightening


For Vdif, only one iteration makes sense because they do not affect convex 
%envelopes.
   
%record initial values of output parameters
theta=branch.theta;
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
    Env = BTfinal_flow_envelopes( bus, branch, theta, Vdif, opt.Vdif_type-1 );
    
    %tighten bounds
    [ theta_new, Vdif_new ] = BTtighten_through_flows( bus, branch, Env, opt, tighten_angle );
    
    %ensure lb<ub
    theta_new=ensure_consistent_bounds(theta_new);
    Vdif_new=ensure_consistent_bounds(Vdif_new);
    
    %compute max changes in bounds compared to previous iteration
    max_bound_change=max([theta.max-theta_new.max, theta_new.min-theta.min, ...
        Vdif.max-Vdif_new.max, Vdif_new.min-Vdif.min]);
    theta_temp=theta;
    Vdif_temp=Vdif;
    
    %record tighter bounds
    theta=BTupdate_bounds(theta, theta_new);
    Vdif=BTupdate_bounds(Vdif, Vdif_new);
    
    %collect iteration-specific statistics if needed
    if (opt.statistics==2)
        stat_iter.theta=[stat_iter.theta; BTcollect_statistics( theta_temp, theta_new, [], 2, 0 )];
        stat_iter.Vdif=[stat_iter.theta; BTcollect_statistics( Vdif_temp, Vdif_new, [], 2, 0 )];  %CHECK WHEN TO USE IT!!!!
    end
    
    %check stopping condition
    if (max_bound_change<eps_tol)
        break;
    else
        iter=iter+1;
    end
end

% postprocess bounds (to store tightest bounds for all parallel branches)
[ theta.min, theta.max ] = BTpostprocess_bounds( theta.min, theta.max, branch.uniq );
[ Vdif.min, Vdif.max ] = BTpostprocess_bounds( Vdif.min, Vdif.max, branch.uniq );
    
end


function bounds=ensure_consistent_bounds(bounds)
    ind_temp=bounds.min>bounds.max;
    lb_temp=bounds.min(ind_temp);
    bounds.min(ind_temp)=bounds.max(ind_temp);
    bounds.max(ind_temp)=lb_temp;
end

