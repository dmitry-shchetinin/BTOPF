function [ theta, Vdif ] = BTbounds_from_one_line_end( branch_data, Vbox, ratio )
%returns the bounds on theta_ij and |Vj-alphaVi| computed from the thermal
%limit constraint at either beginning or end of the line
    
    
%% compute parameters of boundary lines
alpha=branch_data(9);
beta_max=sqrt(branch_data(5)/branch_data(2));
beta_min=-beta_max;


%% record bounds on voltage magnitude differences (|Vj-slope*Vi|<offset)
Vdif=struct('slope',alpha/ratio,'offset',beta_max);


%% compute angle bounds at V_j=V_j_min edge
if (alpha*Vbox.V_i_min+beta_min<=Vbox.V_j_min) %if this edge has feasible points
    Vrange=struct('minval',Vbox.V_i_min,'maxval',Vbox.V_i_max,'minisbound',0,'maxisbound',0);
    Vfixed=struct('val',Vbox.V_j_min,'type',0);
    if (alpha*Vbox.V_i_min+beta_max<Vbox.V_j_min)
        Vrange.minval=(Vbox.V_j_min-beta_max)/alpha;
        Vrange.minisbound=1;
    end
    if (alpha*Vbox.V_i_max+beta_min>Vbox.V_j_min)
        Vrange.maxval=(Vbox.V_j_min-beta_min)/alpha;
        Vrange.maxisbound=1;
    end
    theta_Vj_edge=BTangle_extremum( branch_data, Vrange, Vfixed );
else
    theta_Vj_edge.min=inf;
    theta_Vj_edge.max=-inf;
end

%% compute angle bounds at V_i=V_i_min edge
if (alpha*Vbox.V_i_min+beta_max>=Vbox.V_j_min) %if this edge has feasible points
    Vrange=struct('minval',Vbox.V_j_min,'maxval',Vbox.V_j_max,'minisbound',0,'maxisbound',0);
    Vfixed=struct('val',Vbox.V_i_min,'type',1);
    if (alpha*Vbox.V_i_min+beta_min>Vbox.V_j_min)
        Vrange.minval=alpha*Vbox.V_i_min+beta_min;
        Vrange.minisbound=1;
    end
    if (alpha*Vbox.V_i_min+beta_max<Vbox.V_j_max)
        Vrange.maxval=alpha*Vbox.V_i_min+beta_max;
        Vrange.maxisbound=1;
    end
    theta_Vi_edge=BTangle_extremum( branch_data, Vrange, Vfixed );
else
    theta_Vi_edge.min=inf;
    theta_Vi_edge.max=-inf;
end

%% compute bounds on theta_ij for entire feasible set
theta.min=min(theta_Vi_edge.min,theta_Vj_edge.min);
theta.max=max(theta_Vi_edge.max,theta_Vj_edge.max);

end

