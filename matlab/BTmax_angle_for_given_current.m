function [delta_critical] = BTmax_angle_for_given_current( V_i, V_j, branch_data, surface_part )
%returns the value of theta_ij for point (V_i,V_j) on either upper part of
%the boundary surface (surface_part=1) or its lower part (surface_part=-1).
%If |theta_ij|>pi/2, it is fixed at +-pi/2 since this is steady-state limit
    
%set some intermediate variables for ease of computation
t = branch_data(1) * V_i / V_j + branch_data(2) * V_j / V_i - branch_data(5) / (V_i*V_j);
D = branch_data(6) - t*t;

%check the determinant and compute the sine of angle
if (D < 0) %based on hwo we select point, this means that thermal limit is unattainable
    delta_critical = surface_part*pi/2;
else
    %compute the value of the sine of angle
    sin_delta = branch_data(7) * t + surface_part * branch_data(8) * sqrt(D);

    %compute the angle
    delta_critical = asin(sin_delta);
    delta_critical1 = surface_part * pi - delta_critical;

    %check if the equation is satisfied
    w = abs(branch_data(1) * V_i*V_i + branch_data(2) * V_j*V_j + 2.0 * V_i*V_j*(branch_data(3) * sin(delta_critical) - ...
        branch_data(4) * cos(delta_critical)) - branch_data(5));
    w1 = abs(branch_data(1) * V_i*V_i + branch_data(2) * V_j*V_j + 2.0 * V_i*V_j*(branch_data(3) * sin(delta_critical1) - ...
        branch_data(4) * cos(delta_critical1)) - branch_data(5));
    if (w > w1)
        delta_critical = delta_critical1;
    end
    if (delta_critical>pi/2)
        delta_critical = pi/2;
    elseif (delta_critical<-pi/2)
        delta_critical = -pi/2;
    end
end

end

