function Env = BTfinal_flow_envelopes( bus, branch, theta, use_Vsquared, Vdif )
%returns matrix and bounds on vector that describe linear envelopes of all
%nonlinear terms that enter power flow equations.
    
%unpack some data
[N, L, Vmin, Vmax, ind_bus_F, ind_bus_T, shift]=deal(bus.nbus, branch.nbranch, ...
    bus.Vmin, bus.Vmax, branch.ind_bus_F, branch.ind_bus_T, branch.shift);

%ensure there is distance between angle bounds (for stable computation)
ind_small_angles=(theta.max-theta.min)<1e-7;
theta.min(ind_small_angles)=theta.min(ind_small_angles)-5e-8;
theta.max(ind_small_angles)=theta.max(ind_small_angles)+5e-8;

%envelopes for V^2
varbounds=struct('lb',Vmin,'ub',Vmax);
if (use_Vsquared)
    b1.l=[]; b1.u=[];
else
    [A1, b1]= BTenvelope_for_square(varbounds, N );
end

%envelopes for Vi*Vj
indices=struct('x', ind_bus_F, 'y', ind_bus_T);
dims=struct('nrows',L,'ncols',N);
if (use_Vsquared)
    [ A2, b2 ] = BTenvelope_for_sqrt_xy(indices, varbounds, Vdif, dims);
else
    [ A2, b2 ] = BTenvelope_for_xy(indices, varbounds, dims);
end

%envelopes for sin(theta_ij)
indices=struct('theta_F', ind_bus_F, 'theta_T', ind_bus_T);
varbounds=struct('lb',theta.min-shift,'ub',theta.max-shift);
[ A3, b3 ] = BTenvelope_for_sin(indices, varbounds, dims, shift);

%envelopes for cos(theta_ij)
[ A4, b4 ] = BTenvelope_for_cos(indices, varbounds, dims, shift);

%envelopes for (ViVj)*(sin(theta_ij))
indices=struct('x', (1:L)', 'y', (L+1:2*L)');
dims=struct('nrows',L,'ncols',2*L);
varbounds=struct('lb',[Vmin(ind_bus_F).*Vmin(ind_bus_T); sin(theta.min-shift)],...
    'ub',[Vmax(ind_bus_F).*Vmax(ind_bus_T); sin(theta.max-shift)]);
[A5, b5] = BTenvelope_for_xy(indices, varbounds, dims);

%envelopes for (ViVj)*(cos(theta_ij))
%max of cosine is either 1 or value at delta_min or delta_max if delta_max<0 or delta_min>0
cos_max=ones(L,1);
ind_temp=(theta.min-shift)>0 | (theta.max-shift)<0;
cos_max(ind_temp)=max(cos(theta.min(ind_temp)-shift(ind_temp)),cos(theta.max(ind_temp)-shift(ind_temp)));
varbounds=struct('lb',[Vmin(ind_bus_F).*Vmin(ind_bus_T); cos(theta.min-shift)],...
    'ub',[Vmax(ind_bus_F).*Vmax(ind_bus_T); cos_max]);
[A6, b6] = BTenvelope_for_xy(indices, varbounds, dims);

%construct final envelope matrix and bounds on vector b
if (use_Vsquared)
    A=[speye(N), sparse(N,N), sparse(N,5*L); ...
        A5*[A2, sparse(L,N); sparse(L,N), A3], A5, sparse(L,L), speye(L), sparse(L,L); ...
        A6*[A2, sparse(L,N); sparse(L,N), A4], A6(:,1:L), sparse(L,L), A6(:,L+1:end), sparse(L,L), speye(L)];
else
    A=[A1, sparse(N,N), speye(N), sparse(N,5*L); ...
        A5*[A2, sparse(L,N); sparse(L,N), A3], sparse(L,N), A5, sparse(L,L), speye(L), sparse(L,L); ...
        A6*[A2, sparse(L,N); sparse(L,N), A4], sparse(L,N), A6(:,1:L), sparse(L,L), A6(:,L+1:end), sparse(L,L), speye(L)];
end

%compute matrix A that includes information from envelopes and admittance matrix
A=bus.A0*A;

%delete row corresponding to load bus with small injection and column corresponding to angle of this bus (make it reference and assume it's zero)
A(N+bus.ind_ref,:)=[];
A(:,N+bus.ind_ref)=[];

%record outputs
Env.A=A;
Env.bl=[b1.l; b2.l; b3.l; b4.l; b5.l; b6.l];
Env.bu=[b1.u; b2.u; b3.u; b4.u; b5.u; b6.u];
end

