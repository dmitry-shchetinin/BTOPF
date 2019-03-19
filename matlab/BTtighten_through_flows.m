function [ bounds ] = BTtighten_through_flows( bus, branch, Env, opt, tighten_angle )
%Performs one iteration of bound tightening and returns bounds on either 
%theta (if tighten_angle=1) or Vdif (if tighten_angle=0) for all branches.
    
%unpack data
[ind_bus_F, ind_bus_T, L, N, opt_prob_flow, use_mex, ind_ref, load_PQ, ind_gen, linsolver]=...
    deal(branch.ind_bus_F, branch.ind_bus_T, branch.nbranch, bus.nbus, opt.opt_prob_flow, ...
    opt.use_mex, bus.ind_ref, bus.load_PQ, bus.ind_gen_PQ, opt.linsolver);
    
%extract part of A related to slack variables (envelopes' wiggle room)
A_rhs=-transpose(Env.A(:,2*N:end)); %we will only need transpose

%factorize transpose of part of A related to V and theta
if (linsolver==1)
    [LL, UU, p, q] = lu(transpose(Env.A(:,1:2*N-1)),'vector');    %Alhs(q,p) = L_lhs*U_lhs
    LU_struct=[];
else
    LU_struct=klu(transpose(Env.A(:,1:2*N-1)));
    LL=[]; UU=[]; p=[]; q=[];
end

%record bounds on vector of variables
x_min=[bus.Pmin(bus.ind_genbus); bus.Qmin(bus.ind_genbus); Env.bl];
x_max=[bus.Pmax(bus.ind_genbus); bus.Qmax(bus.ind_genbus); Env.bu];
ind_temp=find(bus.ind_genbus==bus.ind_ref);
if (~isempty(ind_temp)) %remove Qgen of a bus, whose reactive power balance was removed
    x_min(N+ind_temp)=[];
    x_max(N+ind_temp)=[];
end

%record bounds on voltage magnitude at the beginning of each branch
Vmin_busF=bus.Vmin(ind_bus_F);
Vmax_busF=bus.Vmax(ind_bus_F);


%initialize arrays for bounds
bound_min=-inf(L,1);
bound_max=inf(L,1);

%loop over all branches
parfor i=1:L
    i_begin=ind_bus_F(i);
    i_end=ind_bus_T(i);
    
    %compute rhs vector for obtaining required row of inverse
    ind_1=0; ind_2=0;
    if (tighten_angle)
        %get row of inverse corresponding to theta_ij
        if (i_begin<ind_ref)
            ind_1=N+i_begin;
        elseif (i_begin>ind_ref)
            ind_1=N+i_begin-1;
        end
        if (i_end<ind_ref)
            ind_2=N+i_end;
        elseif (i_end>ind_ref)
            ind_2=N+i_end-1;
        end
    else
        %get row of inverse corresponding to Vi-Vj
        ind_1=i_begin;
        ind_2=i_end;
    end
    e_vector=zeros(2*N-1,1);
    if (linsolver==1)
        e_vector(p==ind_1)=1;
        e_vector(p==ind_2)=-1;
        row_bound=UU\(LL\e_vector);
        row_bound(q)=row_bound;
    else
        if (ind_1~=0)
            e_vector(ind_1)=1;
        end
        if (ind_2~=0)
            e_vector(ind_2)=-1;
        end
        row_bound=klu(LU_struct,'\',e_vector);
    end
    
    %get contribution from loads to objective
    objconst=load_PQ*row_bound;
    %get vector of coefficients in the objective function
    c=[row_bound(ind_gen); A_rhs*row_bound];
    
    %compute lower and upper bounds
    if (opt_prob_flow==1) 
        %solve by inspection
        temp=[c.*x_min, c.*x_max];
        bound_min(i)=objconst+sum(min(temp,[],2));
        bound_max(i)=objconst+sum(max(temp,[],2));
    else
        %solve LP with one constraint
        e_vector=zeros(2*N-1,1);
        if (linsolver==1)
            e_vector(p==i_begin)=1;
            row_Vi=UU\(LL\e_vector);
            row_Vi(q)=row_Vi;
        else
           e_vector(i_begin)=1;
           row_Vi=klu(LU_struct,'\',e_vector);
        end
        %get contribution from loads to objective
        b_const=load_PQ*row_Vi;
        c=c';
        %get vector of coefficients in the constraint
        a=[row_Vi(ind_gen); A_rhs*row_Vi]';
        if (use_mex)
            bound_min(i)=objconst-greedy_LP_solution_mex( -c, a, x_min, x_max, Vmin_busF(i)-b_const, Vmax_busF(i)-b_const);
            bound_max(i)=objconst+greedy_LP_solution_mex(  c, a, x_min, x_max, Vmin_busF(i)-b_const, Vmax_busF(i)-b_const);
        else
            bound_min(i)=objconst-greedy_LP_solution( -c, a, x_min, x_max, Vmin_busF(i)-b_const, Vmax_busF(i)-b_const);
            bound_max(i)=objconst+greedy_LP_solution(  c, a, x_min, x_max, Vmin_busF(i)-b_const, Vmax_busF(i)-b_const);
        end
    end
end

%put computed bounds into output structure
bounds=struct('min',bound_min,'max',bound_max);
end


