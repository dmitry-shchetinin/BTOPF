function [ theta, Vdif ] = BTtighten_through_flows( bus, branch, Env, opt, tighten_Vdif_only )
%Performs one iteration of bound tightening and returns bounds on either 
%theta (if tighten_angle=1) or Vdif (if tighten_angle=0) for all branches.
    
%unpack data
[ind_bus_F, ind_bus_T, L, N, opt_prob_flow, use_mex, ind_ref, load_PQ, ind_gen, linsolver, pow]=...
    deal(branch.ind_bus_F, branch.ind_bus_T, branch.nbranch, bus.nbus, opt.opt_prob_flow, ...
    opt.use_mex, bus.ind_ref, bus.load_PQ, bus.ind_gen_PQ, opt.linsolver, opt.Vdif_type);
nrows=2*N-1;

%record information on what should be tightened
tighten_angle=~tighten_Vdif_only;
tighten_Vdif=(tighten_Vdif_only | (opt.bounds~=1 && opt.Vdif_type==2));
    
%extract part of A related to slack variables (envelopes' wiggle room)
A_rhs=-transpose(Env.A(:,2*N:end)); %we will only need transpose

%factorize transpose of part of A related to V and theta
if (linsolver==1)
    [LU_struct.L, LU_struct.U, LU_struct.p, LU_struct.q] = ...
        lu(transpose(Env.A(:,1:nrows)),'vector');    %Alhs(q,p) = L_lhs*U_lhs
else
    LU_struct=klu(transpose(Env.A(:,1:nrows)));
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
Vmin_busF=(bus.Vmin(ind_bus_F)).^pow;
Vmax_busF=(bus.Vmax(ind_bus_F)).^pow;

%initialize arrays for bounds
if (tighten_angle)
    theta_min=-inf(L,1);
    theta_max=inf(L,1);
end
if (tighten_Vdif)
    Vdif_min=-inf(L,1);
    Vdif_max=inf(L,1);
end

%loop over all branches
parfor i=1:L
    i_begin=ind_bus_F(i);
    i_end=ind_bus_T(i);
    
    %compute row of inverse for constraint in LP if need be
    if (opt_prob_flow==2)
        row_Vi=compute_row_of_inverse(LU_struct, i_begin, 0, nrows, linsolver);
    end
    
    %tighten angle bounds
    if (tighten_angle)
        %compute rhs vector for obtaining required row of inverse
        ind_1=0; ind_2=0;
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
        row=compute_row_of_inverse(LU_struct, ind_1, ind_2, nrows, linsolver);
        
        %get contribution from loads to objective
        objconst=load_PQ*row;
        %get vector of coefficients in the objective function
        c=[row(ind_gen); A_rhs*row];
        
        %compute lower and upper bounds
        if (opt_prob_flow==1)
            %solve by inspection
            temp=[c.*x_min, c.*x_max];
            theta_min(i)=objconst+sum(min(temp,[],2));
            theta_max(i)=objconst+sum(max(temp,[],2));
        else
            %solve LP with one constraint
            %get contribution from loads to objective
            b_const=load_PQ*row_Vi;
            c=c';
            %get vector of coefficients in the constraint
            a=[row_Vi(ind_gen); A_rhs*row_Vi]';
            if (use_mex)
                theta_min(i)=objconst-greedy_LP_solution_mex( -c, a, x_min, x_max, Vmin_busF(i)-b_const, Vmax_busF(i)-b_const);
                theta_max(i)=objconst+greedy_LP_solution_mex(  c, a, x_min, x_max, Vmin_busF(i)-b_const, Vmax_busF(i)-b_const);
            else
                theta_min(i)=objconst-greedy_LP_solution( -c, a, x_min, x_max, Vmin_busF(i)-b_const, Vmax_busF(i)-b_const);
                theta_max(i)=objconst+greedy_LP_solution(  c, a, x_min, x_max, Vmin_busF(i)-b_const, Vmax_busF(i)-b_const);
            end
        end
    end
    
    %tighten Vdif bounds
    if (tighten_Vdif)
        row=compute_row_of_inverse(LU_struct, i_begin, i_end, nrows, linsolver);
        
        %get contribution from loads to objective
        objconst=load_PQ*row;
        %get vector of coefficients in the objective function
        c=[row(ind_gen); A_rhs*row];
        
        %compute lower and upper bounds
        if (opt_prob_flow==1)
            %solve by inspection
            temp=[c.*x_min, c.*x_max];
            Vdif_min(i)=objconst+sum(min(temp,[],2));
            Vdif_max(i)=objconst+sum(max(temp,[],2));
        else
            %solve LP with one constraint
            %get contribution from loads to objective
            b_const=load_PQ*row_Vi;
            c=c';
            %get vector of coefficients in the constraint
            a=[row_Vi(ind_gen); A_rhs*row_Vi]';
            if (use_mex)
                Vdif_min(i)=objconst-greedy_LP_solution_mex( -c, a, x_min, x_max, Vmin_busF(i)-b_const, Vmax_busF(i)-b_const);
                Vdif_max(i)=objconst+greedy_LP_solution_mex(  c, a, x_min, x_max, Vmin_busF(i)-b_const, Vmax_busF(i)-b_const);
            else
                Vdif_min(i)=objconst-greedy_LP_solution( -c, a, x_min, x_max, Vmin_busF(i)-b_const, Vmax_busF(i)-b_const);
                Vdif_max(i)=objconst+greedy_LP_solution(  c, a, x_min, x_max, Vmin_busF(i)-b_const, Vmax_busF(i)-b_const);
            end
        end
    end
end

%put computed bounds into output structures
if (tighten_angle)
    theta=struct('min',theta_min,'max',theta_max);
else
    theta=[];
end
if (tighten_Vdif)
    Vdif=struct('min',Vdif_min,'max',Vdif_max);
else
    Vdif=[];
end
end


%compute row (or difference between two rows) of matrix inverse
function row=compute_row_of_inverse(LUstruct, ind1, ind2, nrows, linsolver)
    e_vector=zeros(nrows,1);
    if (linsolver==1)
        e_vector(LUstruct.p==ind1)=1;
        e_vector(LUstruct.p==ind2)=-1;
        row=LUstruct.U\(LUstruct.L\e_vector);
        row(LUstruct.q)=row;
    else
        if (ind1~=0)
            e_vector(ind1)=1;
        end
        if (ind2~=0)
            e_vector(ind2)=-1;
        end
        row=klu(LUstruct,'\',e_vector);
    end
end


