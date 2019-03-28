function [ theta, Vdif, stat_iter ] = BTbounds_from_flow_envelopes( bus, branch, Vdif_in, opt, tighten_Vdif_only )
%Performs iterative bound tightening and returns bounds on theta, Vdif, or 
%both theta and Vdif for all branches. If voltage magnitudes are used as
%variables, angles can be tightened iteratively and for Vdif, only one
%iteration is possible. If squares of voltage magnitudes are used as
%variables, simultaneously Vdif and angles can be tightened iteratively, or
%only Vdif or theta can be tightened iteratively. 
   
%record initial values of output parameters
theta=branch.theta;
Vdif=Vdif_in;
stat_iter=struct('theta',[],'Vdif',[]);

%adjust parameters depending on what type of bounds is tightened
if (tighten_Vdif_only && opt.Vdif_type==1)
    max_iter=1;
else
    max_iter=opt.max_iter;
end

%start iterative process of tightening bounds
iter=1;
while iter<=max_iter
    %build all envelopes
    Env = BTfinal_flow_envelopes( bus, branch, theta, Vdif, opt.Vdif_type-1 );
    
    %tighten bounds
    [ theta_new, Vdif_new ] = BTtighten_through_flows( bus, branch, Env, opt, tighten_Vdif_only );

    %check which bound types were tightened
    if (isempty(theta_new))
        theta_new=theta;
    end
    if (isempty(Vdif_new))
        Vdif_new=Vdif;
    end
    
    %compute max changes in bounds compared to previous iteration
    max_bound_change=max([theta.max-theta_new.max, theta_new.min-theta.min, ...
        Vdif.max-Vdif_new.max, Vdif_new.min-Vdif.min]);
    theta_temp=theta;
    Vdif_temp=Vdif;
    
    %record tighter bounds
    theta=BTupdate_bounds(theta, theta_new);
    Vdif=BTupdate_bounds(Vdif, Vdif_new);
    
    %ensure lb<ub
    theta=ensure_consistent_bounds(theta);
    Vdif=ensure_consistent_bounds(Vdif);
    
    %collect iteration-specific statistics if needed
    if (opt.statistics==2)
        stat_iter.theta=[stat_iter.theta; BTcollect_statistics( theta_temp, theta, [], 2, 0 )];
        stat_iter.Vdif=[stat_iter.Vdif; BTcollect_statistics( Vdif_temp, Vdif, [], 2, 0 )];
    end
    
    %check stopping condition
    if (max_bound_change<opt.max_change)
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

