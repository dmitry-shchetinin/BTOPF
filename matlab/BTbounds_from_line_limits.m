function [ theta, Vdif ] = BTbounds_from_line_limits( bus, branch )
%Returns bounds on angle and voltage differences for all lines based on the
%feasible set of thermal limit constraints. If a line does not have a
%thermal limit constraint (modeled as branch.I_max=0), the limits on phase
%andgle differences are +-pi/2, and limits on voltage differences are
%taken from the bounds on voltages
    
%% initialize output structures
L=numel(branch.g);
theta=struct('min',-pi*ones(L,1)/2,'max',pi*ones(L,1)/2);
Vdif=struct('min',bus.Vmin(branch.ind_bus_F)-bus.Vmax(branch.ind_bus_T),...
        'max',bus.Vmax(branch.ind_bus_F)-bus.Vmin(branch.ind_bus_T));
extra=struct('slope1',zeros(L,1),'offset1',zeros(L,1),'slope2',zeros(L,1),'offset2',zeros(L,1));


%iterate over all branches
for i=1:L
    if (branch.I_max(i)>0)
        ratio=branch.ratio(i);

        %% obtain limits on V_i and V_j (after transformer)
        Vbox=struct('V_i_min',bus.Vmin(branch.ind_bus_F(i))/ratio,...
            'V_j_min',bus.Vmin(branch.ind_bus_T(i)),...
            'V_i_max',bus.Vmax(branch.ind_bus_F(i))/ratio,...
            'V_j_max',bus.Vmax(branch.ind_bus_T(i)));

        %% deal with the beginning of the line
        branch_params=BTbranch_parameters(branch.g(i), branch.b(i), branch.bs(i), ratio, branch.I_max(i), 1);
        [ theta_begin, Vdif_begin ] = BTbounds_from_one_line_end(branch_params, Vbox, ratio);

        %% deal with the end of the line
        branch_params=BTbranch_parameters(branch.g(i), branch.b(i), branch.bs(i), ratio, branch.I_max(i), 2);
        [ theta_end, Vdif_end ] = BTbounds_from_one_line_end(branch_params, Vbox, ratio);

        %% record final angle limits
        theta.min(i)=min(theta_begin.min,theta_end.min)+branch.shift(i);
        theta.max(i)=max(theta_begin.max,theta_end.max)+branch.shift(i);
        
        
        %% record voltage limits in form |Vj-alpha*Vi|<beta
        [extra.slope1(i), extra.offset1(i), extra.slope2(i), extra.offset2(i)]=deal(Vdif_begin.slope, ...
            Vdif_begin.offset, Vdif_end.slope, Vdif_end.offset);
        
        %% do postprocessing to obtain |Vi-Vj| constraints
        Vbox.V_i_min=Vbox.V_i_min*ratio;
        Vbox.V_i_max=Vbox.V_i_max*ratio;
        if (extra.slope1(i)>=extra.slope2(i))
            Vdif.max(i)=compute_bound_Vdif(extra.slope1(i), -extra.offset1(i), Vbox, Vdif.max(i));
            Vdif.min(i)=compute_bound_Vdif(extra.slope2(i), extra.offset2(i), Vbox, Vdif.min(i));
        else
            Vdif.max(i)=compute_bound_Vdif(extra.slope2(i), -extra.offset2(i), Vbox, Vdif.max(i));
            Vdif.min(i)=compute_bound_Vdif(extra.slope1(i), extra.offset1(i), Vbox, Vdif.min(i));
         end
    end
end

%% record all obtained voltage constraints
Vdif.extra=extra;

%% postprocess bounds (to store tightest bounds for all parallel branches)
[ theta.min, theta.max ] = BTpostprocess_bounds( theta.min, theta.max, branch.uniq );
[ Vdif.min, Vdif.max ] = BTpostprocess_bounds( Vdif.min, Vdif.max, branch.uniq );
    
end


%% returns bound on |Vi-Vj| from Vj-slope*Vi=offset constraint
function bound_out=compute_bound_Vdif(slope, offset, Vbox, bound_in)
    %check if line intersects Vi=Vi_min edge
    Vj_temp=slope*Vbox.V_i_min+offset;
    if (Vj_temp>=Vbox.V_j_min && Vj_temp<=Vbox.V_j_max)
        bound_out=Vbox.V_i_min-Vj_temp;
        return;
    end
    
    %check if line intersects Vj=Vj_min edge
    Vi_temp=(Vbox.V_j_min-offset)/slope;
    if (Vi_temp>=Vbox.V_i_min && Vi_temp<=Vbox.V_i_max)
        bound_out=Vi_temp-Vbox.V_j_min;
        return;
    end
    
    %if line does not intersect with either, set bound to initial value
    bound_out=bound_in;
end


