function [ A, b ] = BTenvelope_for_xy( indices, varbounds, dims )
%returns sparse matrix and structure with vectors of lower and upper bounds
%for linear envelopes of z=xy function. Each envelope is represented as
%z=a1*x+a2*y+b, l<b<u
    
%matrix
ind_rows=repmat((1:dims.nrows)',2,1);
ind_cols=[indices.x; indices.y];
ind_vals=[varbounds.lb(indices.y)+varbounds.ub(indices.y); varbounds.lb(indices.x)+varbounds.ub(indices.x)]/2;
A=sparse(ind_rows, ind_cols, ind_vals, dims.nrows, dims.ncols);
%bounds
b.l=-((varbounds.lb(indices.x)).*(varbounds.lb(indices.y))+(varbounds.ub(indices.x)).*(varbounds.ub(indices.y)))/2;
b.u=-((varbounds.ub(indices.x)).*(varbounds.lb(indices.y))+(varbounds.lb(indices.x)).*(varbounds.ub(indices.y)))/2;
    
end