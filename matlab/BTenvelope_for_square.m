function [ A, b ] = BTenvelope_for_square( varbounds, nrows )
%returns sparse matrix and structure with vectors of lower and upper bounds
%for linear envelopes of x^2 function. Each envelope is represented as
%y=ax+b, l<b<u
    
%matrix
A=sparse(1:nrows,1:nrows,varbounds.lb+varbounds.ub, nrows, nrows);
%offsets
b.l=-0.25*(varbounds.lb+varbounds.ub).^2;
b.u=-(varbounds.lb).*(varbounds.ub);
    
end