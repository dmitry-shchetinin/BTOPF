function [ A, b ] = BTenvelope_for_cos(indices, varbounds, dims, shifts)
%returns sparse matrix and structure with vectors of lower and upper bounds
%for linear envelopes of cosine function. Each envelope is represented as
%y=ax+b, l<b<u

%compute vector of coefficients to be used
coef_temp=(cos(varbounds.lb)-cos(varbounds.ub))./(varbounds.lb-varbounds.ub);

%create sparse matrix
ind_rows=repmat((1:dims.nrows)',2,1);
ind_cols=[indices.theta_F; indices.theta_T];
ind_vals=[coef_temp; -coef_temp];
A=sparse(ind_rows, ind_cols, ind_vals, dims.nrows, dims.ncols);

%record lower bounds
b.l=cos(varbounds.lb)-coef_temp.*(varbounds.lb+shifts);

%work on upper bounds
x_tangent=-asin(coef_temp)+shifts; %tangent point
b.u=cos(x_tangent-shifts)-coef_temp.*x_tangent;
end

