function branch_data = BTbranch_parameters( g, b, bs, ratio, I_max, flow_side )
%returns some coefficients derived from branch parameters to accelerate
%subsequent computations
    
%take transformer into account
if (flow_side == 1)
    I_max = I_max*ratio;
end
bs = bs/2;

%compute frequenly used values derived from g, b, b_sh
if (flow_side == 1)
    branch_data(1) = g^2 + (b + bs)^2;
    branch_data(2) = g^2 + b^2;
    branch_data(3) = g * bs;
else
    branch_data(1) = g^2 + b^2;
    branch_data(2) = g^2 + (b + bs)^2;
    branch_data(3) = -g * bs;
end
branch_data(4) = g^2 + b^2 + b*bs;
branch_data(5) = I_max^2;
branch_data(6) = 4 * branch_data(1) * branch_data(2);
branch_data(7) = -branch_data(3) / (2 * branch_data(1) * branch_data(2));
branch_data(8) = abs(branch_data(4)) / (2 * branch_data(1) * branch_data(2));
branch_data(9) = sqrt(branch_data(1) / branch_data(2));

end

