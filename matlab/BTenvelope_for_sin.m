function [ A, b ] = BTenvelope_for_sin( indices, varbounds, dims, shifts )
%returns sparse matrix and structure with vectors of lower and upper bounds
%for linear envelopes of sine function. Each envelope is represented as
%y=ax+b, l<b<u
    
%initialize output structure
b=struct('l',zeros(dims.nrows,1),'u',zeros(dims.nrows,1));
temp=(1:dims.nrows)';

%deal with variables with lb>=0
ind_temp=varbounds.lb>=0; ind_pos=ind_temp;
coef_temp=(sin(varbounds.lb(ind_temp))-sin(varbounds.ub(ind_temp)))./(varbounds.lb(ind_temp)-varbounds.ub(ind_temp));
rows_pos=repmat(temp(ind_temp),2,1);
cols_pos=[indices.theta_F(ind_temp); indices.theta_T(ind_temp)];
vals_pos=[coef_temp; -coef_temp];
b.l(ind_temp)=sin(varbounds.lb(ind_temp))-coef_temp.*(varbounds.lb(ind_temp)+shifts(ind_temp));
x_tangent=acos(coef_temp)+shifts(ind_temp);
b.u(ind_temp)=sin(x_tangent-shifts(ind_temp))-coef_temp.*x_tangent;

%deal with variables with ub<=0
ind_temp=varbounds.ub<=0 & ~ind_pos; ind_neg=ind_temp;
coef_temp=(sin(varbounds.lb(ind_temp))-sin(varbounds.ub(ind_temp)))./(varbounds.lb(ind_temp)-varbounds.ub(ind_temp));
rows_neg=repmat(temp(ind_temp),2,1);
cols_neg=[indices.theta_F(ind_temp); indices.theta_T(ind_temp)];
vals_neg=[coef_temp; -coef_temp];
b.l(ind_temp)=sin(varbounds.lb(ind_temp))-coef_temp.*(varbounds.lb(ind_temp)+shifts(ind_temp));
x_tangent=-acos(coef_temp)+shifts(ind_temp);
b.u(ind_temp)=sin(x_tangent-shifts(ind_temp))-coef_temp.*x_tangent;

%deal with rest
ind_temp=~ind_pos & ~ind_neg;
coef_temp=max(abs(varbounds.lb(ind_temp)),abs(varbounds.ub(ind_temp)))/2;
rows_norm=repmat(temp(ind_temp),2,1);
cols_norm=[indices.theta_F(ind_temp); indices.theta_T(ind_temp)];
vals_norm=[cos(coef_temp); -cos(coef_temp)];
b.l(ind_temp)=-sin(coef_temp)-cos(coef_temp).*(shifts(ind_temp)-coef_temp);
b.u(ind_temp)=sin(coef_temp)-cos(coef_temp).*(shifts(ind_temp)+coef_temp);

%build matrix
A=sparse([rows_pos; rows_neg; rows_norm], ...
    [cols_pos; cols_neg; cols_norm], ...
    [vals_pos; vals_neg; vals_norm], dims.nrows, dims.ncols);
end