function [ theta ] = BTangle_extremum( branch_data, Vrange, Vfixed )
%returns lower and upper bounds on theta_ij for a given edge of ViVj box
    
%compute value of V for which slope is zero at the given edge 
if (Vfixed.type)
    V=sqrt((branch_data(1)*Vfixed.val^2-branch_data(5))/branch_data(2));
else
    V=sqrt((branch_data(2)*Vfixed.val^2-branch_data(5))/branch_data(1));
end

%check if this V value is within feasible range
if (isreal(V) && V>=Vrange.minval && V<=Vrange.maxval)
    %compute theta_ij for this V at lower and upper boundary part
    if (Vfixed.type)
        theta.min=BTmax_angle_for_given_current(Vfixed.val, V, branch_data, -1);
        theta.max=BTmax_angle_for_given_current(Vfixed.val, V, branch_data, 1);
    else
        theta.min=BTmax_angle_for_given_current(V, Vfixed.val, branch_data, -1);
        theta.max=BTmax_angle_for_given_current(V, Vfixed.val, branch_data, 1);
    end
else
    %lower part: compute theta_ij values at range ends and choose lower
    if (Vrange.minisbound)
        theta_Vmin_lower=inf;
        theta_Vmin_upper=-inf;
    else
        if (Vfixed.type)
            theta_Vmin_lower=BTmax_angle_for_given_current(Vfixed.val, Vrange.minval, branch_data, -1);
            theta_Vmin_upper=BTmax_angle_for_given_current(Vfixed.val, Vrange.minval, branch_data, 1);
        else
            theta_Vmin_lower=BTmax_angle_for_given_current(Vrange.minval, Vfixed.val, branch_data, -1);
            theta_Vmin_upper=BTmax_angle_for_given_current(Vrange.minval, Vfixed.val, branch_data, 1);
        end
    end
    if (Vrange.maxisbound)
        theta_Vmax_lower=inf;
        theta_Vmax_upper=-inf;
    else
        if (Vfixed.type)
            theta_Vmax_lower=BTmax_angle_for_given_current(Vfixed.val, Vrange.maxval, branch_data, -1);
            theta_Vmax_upper=BTmax_angle_for_given_current(Vfixed.val, Vrange.maxval, branch_data, 1);
        else
            theta_Vmax_lower=BTmax_angle_for_given_current(Vrange.maxval, Vfixed.val, branch_data, -1);
            theta_Vmax_upper=BTmax_angle_for_given_current(Vrange.maxval, Vfixed.val, branch_data, 1);
        end
    end
    theta.min=min(theta_Vmin_lower,theta_Vmax_lower);
    theta.max=max(theta_Vmin_upper,theta_Vmax_upper);
end

end