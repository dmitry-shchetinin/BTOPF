function [ theta ] = BTbounds_from_inj_envelopes( bus, branch, opt )
%returns bounds on angle differences obtained by BT based on convex
%envelopes of bus current injections
    
%unpack data
[N, L, L_uniq, LU_solver, ordering, adj_buses, adj_buses_all, adj_branches, ...
    n_adj_buses, adj_signs, opt_prob_inj, uniq]=deal(bus.nbus, branch.nbranch, ...
    branch.uniq.L_uniq, opt.linsolver, bus.ordering, bus.adj_buses, ...
    bus.adj_buses_all, bus.adj_branches, bus.n_adj_buses, bus.adj_signs, ...
    opt.opt_prob_inj, branch.uniq);
LP_opt=struct('max_iter',opt.max_iter,'eps_tol',opt.max_change,'use_mex',opt.use_mex);

%reorder bus-related vectors
Pmin=bus.Pmin(ordering);
Pmax=bus.Pmax(ordering);
Qmin=bus.Qmin(ordering);
Qmax=bus.Qmax(ordering);
Vmin=bus.Vmin(ordering);
Vmax=bus.Vmax(ordering);
Ybus=bus.Ybus(ordering,ordering);

%get bounds on angles for unique connections
theta_min=branch.theta_min(uniq.full_to_uniq);
theta_max=branch.theta_max(uniq.full_to_uniq);

%for each bus, record some info about adjacent buses/branches to enable parfor loop
adj_Vmin=cell(N,1);
adj_Vmax=cell(N,1);
adj_theta_min=cell(N,1);
adj_theta_max=cell(N,1);
for i=1:N
    adj_Vmin{i}=Vmin(adj_buses{i});
    adj_Vmax{i}=Vmax(adj_buses{i});
    adj_theta_min{i}=theta_min(adj_branches{i});
    adj_theta_max{i}=theta_max(adj_branches{i});
end

%compute values of bus shunts that will replace part of the load
g_shunt=-(Pmin+Pmax)./(2*(Vmin).*(Vmax));
b_shunt=(Qmin+Qmax)./(2*(Vmin).*(Vmax));

%compute bounds on magnitudes of bus current injections
I_bound = BTinj_current_bound_circle( Vmin.^2, Vmax.^2, Pmin, Pmax, Qmin, Qmax, g_shunt, b_shunt );

%update admittance matrix
Ybus=transpose(Ybus+sparse(1:N, 1:N, complex(g_shunt,b_shunt), N, N)); %we'll only need the transpose

%initizalize arrays
ind_all=(1:N)';
V_imag_min=cell(N,1);
V_imag_max=cell(N,1);

%loop over all nodes
for i=1:N
    %get the list of buses that bus i is connected to
    N_adj=n_adj_buses(i);
    buses=adj_buses{i};
    buses_all=adj_buses_all{i};
    
    %check if this bus has any connections that have not been dealt with
    if (N_adj==0)
        continue;
    end
    
    %update i-th diagonal element of the matrix
    Y_temp=Ybus-sparse(i, i, complex(g_shunt(i),b_shunt(i)), N, N);
    
    %delete i-th column from the matrix (i-th row from its transpose) and record it in a separate vector
    Ycol=transpose(Y_temp(i,:));
    Y_temp(i,:)=[];
    
    %delete the row (column of transpose) corresponding to any bus not connected to bus i
    if (median(buses_all)<N/2)
        i_start=N; i_end=1; i_step=-1;
    else
        i_start=1; i_end=N; i_step=1;
    end
    for counter=i_start:i_step:i_end
        row_ind=counter;
        if (all(buses_all~=row_ind)==1 && row_ind~=i)
            break;
        end
    end
    Y_temp(:,row_ind)=[];
    
    %factorize matrix Y_temp into LU
    if (LU_solver==1)
        [LL,UU,p,q] = lu(Y_temp,'vector');    %Y_temp(q,p) = LL*UU
    else
        LU_struct=klu(Y_temp);
    end
    
    %initialize vectors for current injections limits in all nodes
    [I_real_min, I_imag_min, I_real_max, I_imag_max]=deal(-I_bound, -I_bound, ...
        I_bound, I_bound);
    
    %update bounds on I for nodes adjacent with node i
    temp=-[real(Ycol(buses_all))*Vmin(i), real(Ycol(buses_all))*Vmax(i)]; 
    I_real_min(buses_all)=I_real_min(buses_all)+min(temp,[],2);
    I_real_max(buses_all)=I_real_max(buses_all)+max(temp,[],2);
    temp=-[imag(Ycol(buses_all))*Vmin(i), imag(Ycol(buses_all))*Vmax(i)]; 
    I_imag_min(buses_all)=I_imag_min(buses_all)+min(temp,[],2);
    I_imag_max(buses_all)=I_imag_max(buses_all)+max(temp,[],2);
    
    %update bounds on I for node i
    [I_real_min(i),I_real_max(i)]=BTinj_current_bound_box(Pmin(i),Pmax(i),Vmin(i),Vmax(i),real(Ycol(i)));
    [I_imag_min(i),I_imag_max(i)]=BTinj_current_bound_box(-Qmax(i),-Qmin(i),Vmin(i),Vmax(i),imag(Ycol(i)));
    
    %record indices of nonzero entries of dI
    ind_nonzeros=I_real_min~=0 | I_real_max~=0 | I_imag_min~=0 | I_imag_max~=0;
    
    %loop over those adjacent buses, the connections with which were untouched before
    V_imag_min{i}=-inf(N_adj,1);
    V_imag_max{i}=inf(N_adj,1);
    for counter=1:N_adj
        %obtain index of the row of inverse of Y_temp corresponding to Vj
        inv_row_ind=buses(counter)-(buses(counter)>i);
        
        %get row of inverse of Y_temp corresponding to Vj
        e_vector=zeros(N-1,1);
        if (LU_solver==1)
            e_vector(p==inv_row_ind)=1;
            a_vector=UU\(LL\e_vector);
            a_vector(q)=a_vector;
        else
            e_vector(inv_row_ind)=1;
            a_vector=klu(LU_struct,'\',e_vector);
        end
        
        %insert zero element corresponding to the deleted entry of delta_I
        a_vector=[a_vector(1:row_ind-1); 0; a_vector(row_ind:end)];
        a_real=real(a_vector); 
        a_imag=imag(a_vector);
        
        %initialize values of the bounds
        V_min_SOCP=-inf; V_max_SOCP=inf;
        V_min_LP=-inf; V_max_LP=inf;
        
        %solve SOCP if need be
        if (opt_prob_inj~=2)
            ind_connected=[buses_all;i];
            %compute contributions from all nodes not connected to node i
            ind_temp=~ismember(ind_all, ind_connected) & ind_nonzeros;
            obj_from_cones=sum(abs(a_vector(ind_temp)).*(I_bound(ind_temp)));
            %compute contributions from adjacent nodes and node i
            temp=[[a_imag(ind_connected).*I_real_min(ind_connected), a_imag(ind_connected).*I_real_max(ind_connected)]; ...
                [a_real(ind_connected).*I_imag_min(ind_connected), a_real(ind_connected).*I_imag_max(ind_connected)]];
            V_min_SOCP=sum(min(temp,[],2))-obj_from_cones;
            V_max_SOCP=sum(max(temp,[],2))+obj_from_cones;
        end

        %solve LP with one constraint if need be
        if (opt_prob_inj~=1)
            %solve LP: max/min Vj_real s.t.  0<Vj_imag<Vj_max
            if (adj_signs{i}(counter)==1) %line is from i to j
                theta_temp_min=-adj_theta_max{i}(counter);
                theta_temp_max=-adj_theta_min{i}(counter);
            else
                theta_temp_min=adj_theta_min{i}(counter);
                theta_temp_max=adj_theta_max{i}(counter);
            end
            
            LP_prob=struct('c',[a_imag(ind_nonzeros); a_real(ind_nonzeros)]',...
                'a',[a_real(ind_nonzeros); -a_imag(ind_nonzeros)]',...
                'xl',[I_real_min(ind_nonzeros); I_imag_min(ind_nonzeros)],...
                'xu',[I_real_max(ind_nonzeros); I_imag_max(ind_nonzeros)],...
                'bu',adj_Vmax{i}(counter));
            
            %for max(Vj_imag)
            LP_prob.bl=adj_Vmin{i}(counter)*cos(theta_temp_max);
            V_max_LP=BTiterative_LP_for_inj_env(LP_prob,adj_Vmin{i}(counter),LP_opt);
            
            %for min(Vj_imag)
            LP_prob.c=-LP_prob.c;
            LP_prob.bl=adj_Vmin{i}(counter)*cos(theta_temp_min);
            V_min_LP=-BTiterative_LP_for_inj_env(LP_prob,adj_Vmin{i}(counter),LP_opt);
        end
        
        %record values of the bounds
        V_imag_min{i}(counter)=max(V_min_SOCP, V_min_LP);
        V_imag_max{i}(counter)=min(V_max_SOCP, V_max_LP);
    end
end


%compute bounds on angles
theta_min=zeros(L_uniq,1);
theta_max=zeros(L_uniq,1);
for i=1:N
    ind_temp=adj_signs{i}==1;
    theta_min(adj_branches{i}(ind_temp))=-asin(V_imag_max{i}(ind_temp)./adj_Vmin{i}(ind_temp));
    theta_max(adj_branches{i}(ind_temp))=-asin(V_imag_min{i}(ind_temp)./adj_Vmin{i}(ind_temp));
    ind_temp=adj_signs{i}==-1;
    theta_min(adj_branches{i}(ind_temp))=asin(V_imag_min{i}(ind_temp)./adj_Vmin{i}(ind_temp));
    theta_max(adj_branches{i}(ind_temp))=asin(V_imag_max{i}(ind_temp)./adj_Vmin{i}(ind_temp));
end

%limit unrealistic bounds by +-pi/2
theta_min(imag(theta_min)~=0)=-pi/2;
theta_max(imag(theta_max)~=0)=pi/2;

%make sure theta_max>theta_min
ind_temp=theta_min>theta_max;
temp=theta_min(ind_temp);
theta_min(ind_temp)=theta_max(ind_temp);
theta_max(ind_temp)=temp;

%map obtained bounds into list of all branches
theta=struct('min', zeros(L,1), 'max', zeros(L,1));
ind_temp=uniq.signs==1;
theta.min(ind_temp)=theta_min(uniq.ind_original(ind_temp));
theta.max(ind_temp)=theta_max(uniq.ind_original(ind_temp));
ind_temp=uniq.signs==-1;
theta.min(ind_temp)=-theta_max(uniq.ind_original(ind_temp));
theta.max(ind_temp)=-theta_min(uniq.ind_original(ind_temp));
%no postprocessing is needed

end
