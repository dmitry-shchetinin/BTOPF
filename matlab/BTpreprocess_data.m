function [bus, branch]=BTpreprocess_data(mpc, opt)
%get data required by the algorithms to run 

%% remove isolated and out-of-service gens, buses, and branches
buses_out=mpc.bus(mpc.bus(:,2)==4,1);
if isempty(buses_out)
    buses_out=-inf;
end
branches_out=mpc.branch(:,11)==0 | ismember(buses_out,mpc.branch(:,1)) | ...
    ismember(buses_out,mpc.branch(:,2));
mpc.gen(ismember(buses_out,mpc.gen(:,1)) | mpc.gen(:,8)<=0,:) = [];   
mpc.bus(mpc.bus(:,2)==4,:)=[];
mpc.branch(branches_out,:)=[];


%% apply consecutive bus numbering
bus_number=mpc.bus(:,1);
N=numel(bus_number);
mpc.bus(:,1)=(1:N);
e2i=sparse(bus_number,ones(N,1),1:N,max(bus_number), 1); 
mpc.branch(:,1)=full(e2i(mpc.branch(:,1))); 
mpc.branch(:,2)=full(e2i(mpc.branch(:,2))); 
[~,mpc.gen(:,1)]=ismember(mpc.gen(:,1),bus_number);


%% record data required by all methods
%for buses
bus=struct('Vmin',mpc.bus(:,13),'Vmax',mpc.bus(:,12));
bus.nbus=size(mpc.bus,1);

%for branches
y_temp=1./complex(mpc.branch(:,3),mpc.branch(:,4));
branch=struct('g',real(y_temp),'b',imag(y_temp),'ratio',mpc.branch(:,9),...
    'shift',mpc.branch(:,10)*pi/180,'ind_bus_F',mpc.branch(:,1), ...
    'ind_bus_T',mpc.branch(:,2), 'in_service', find(~branches_out),...
    'theta_min',mpc.branch(:,12)*pi/180,'theta_max',mpc.branch(:,13)*pi/180);
branch.ratio(branch.ratio==0)=1;
branch.theta_min(branch.theta_min<-pi/2)=-pi/2; %fix big angles at steady-state limit
branch.theta_max(branch.theta_max>pi/2)=pi/2; %fix big angles at steady-state limit
branch.nbranch=size(mpc.branch,1);
branch.Vdif=mpc.bus(branch.ind_bus_F,8)-mpc.bus(branch.ind_bus_T,8);
branch.theta=(mpc.bus(branch.ind_bus_F,9)-mpc.bus(branch.ind_bus_T,9))*pi/180;

%get unique connections between buses
[~,uniq.full_to_uniq,uniq.ind_original]=unique([(branch.ind_bus_F).*(branch.ind_bus_T), ...
    branch.ind_bus_F+branch.ind_bus_T], 'rows','stable');
uniq.ind_F=branch.ind_bus_F(uniq.full_to_uniq);
uniq.ind_T=branch.ind_bus_T(uniq.full_to_uniq);
uniq.L_uniq=numel(uniq.ind_F);

%for each index in the list of unique connection, find indices of all
%corresponding elements in the original list of branches
[~,uniq.ordering]=sort(uniq.ind_original);
[uniq.n_entries,~] = histcounts(uniq.ind_original,'BinMethod','integers');
uniq.signs=ones(branch.nbranch,1); %if 1, branch is from i to j, -1 otherwise
counter=0;
for i=1:uniq.L_uniq
    ind_temp=uniq.ordering(counter+1:counter+uniq.n_entries(i));
    uniq.signs(ind_temp(branch.ind_bus_F(ind_temp)~=uniq.ind_F(i)))=-1;
    counter=counter+uniq.n_entries(i);
end
branch.uniq=uniq;



%% record data for tightening based on line flow constraints
if (opt.method==1 || opt.method==2)
    branch.bs=mpc.branch(:,5);
    branch.I_max=mpc.branch(:,6)/mpc.baseMVA;
end


%% record common data for tightening based on envelopes
if (opt.method~=2)
    [bus.Ybus, ~, ~] = makeYbus(mpc.baseMVA, mpc.bus, mpc.branch);
    mpc = cleanup_Qgen(mpc); %attempt to remove unrealistic bounds on Qgen
    [~, bus.ind_genbus]=ismember(unique(mpc.gen(:,1)),mpc.bus(:,1)); %indexes of nodes that have generators in the node list
    bus.ngen=size(mpc.gen(:,1),1);
    %compute min and max bus injections
    Cg=sparse(mpc.gen(:,1), (1:bus.ngen)', ones(bus.ngen, 1), bus.nbus, bus.ngen);
    bus.Pmin=(Cg*mpc.gen(:,10)-mpc.bus(:,3))/mpc.baseMVA;
    bus.Pmax=(Cg*mpc.gen(:,9)-mpc.bus(:,3))/mpc.baseMVA;
    bus.Qmin=(Cg*mpc.gen(:,5)-mpc.bus(:,4))/mpc.baseMVA;
    bus.Qmax=(Cg*mpc.gen(:,4)-mpc.bus(:,4))/mpc.baseMVA;
end


%% record data for tightening based on injections envelopes
if (opt.method==1 || opt.method==4)
    [ bus, branch ] = BTdata_for_inj_env( bus, branch );
end


%% record data for tightening based on flow envelopes
if (opt.method==1 || opt.method==3)
    N=bus.nbus;
    L=branch.nbranch;
    %construct admittance-related matrix to linearize power flow equations
    Y_temp.G_f=sparse(branch.ind_bus_F, 1:L, (branch.g)./(branch.ratio), N, L);
    Y_temp.G_t=sparse(branch.ind_bus_T, 1:L, (branch.g)./(branch.ratio), N, L);
    Y_temp.B_f=sparse(branch.ind_bus_F, 1:L, (branch.b)./(branch.ratio), N, L);
    Y_temp.B_t=sparse(branch.ind_bus_T, 1:L, (branch.b)./(branch.ratio), N, L);
    bus.A0=[diag(diag(real(bus.Ybus))), -(Y_temp.B_f-Y_temp.B_t), -(Y_temp.G_f+Y_temp.G_t); ...
        -diag(diag(imag(bus.Ybus))), -(Y_temp.G_f-Y_temp.G_t), Y_temp.B_f+Y_temp.B_t];
    
    %choose reference bus
    bus.ind_ref = BTreference_bus( mpc );
    
    %record indices of power outputs of generators in the list of injections (without reference bus)
    bus.ind_gen_PQ=[bus.ind_genbus; N+[bus.ind_genbus(bus.ind_genbus<bus.ind_ref); bus.ind_genbus(bus.ind_genbus>bus.ind_ref)-1]];
    
    %record load injections (for nodes with generators set them to zero)
    bus.load_PQ=-transpose([mpc.bus(:,3); mpc.bus(:,4)])/mpc.baseMVA;
    bus.load_PQ(N+bus.ind_ref)=[];
    bus.load_PQ(bus.ind_gen_PQ)=0;
end


%% order structure fields by name for convenience
bus=orderfields(bus);
branch=orderfields(branch);

end


%% adjust unrealistic bounds on Q for generators
function [ mpc ] = cleanup_Qgen( mpc )
M=length(mpc.gen(:,1));
for i=1:M
    Q_min=mpc.gen(i,5)/mpc.baseMVA;
    Q_max=mpc.gen(i,4)/mpc.baseMVA;
    P_max=max(abs(mpc.gen(i,9)),abs(mpc.gen(i,8)))/mpc.baseMVA;
    if abs(Q_max)>10 && abs(Q_max)>P_max
        mpc.gen(i,4)=max([P_max*mpc.baseMVA,5*abs(mpc.gen(i,3)),mpc.baseMVA]);
    end
    if abs(Q_min)>10 && abs(Q_min)>P_max
        mpc.gen(i,5)=-max([P_max*mpc.baseMVA,5*abs(mpc.gen(i,3)),mpc.baseMVA]);
    end
end
end





%% TESTING:record angles and voltages
% dif_V=0.0001;
% dif_theta=0.0001;
% dif_PQ=0.0001;
% branch.theta_min=(mpc.bus(branch.ind_bus_F,9)-mpc.bus(branch.ind_bus_T,9))*pi/180-dif_theta;
% branch.theta_max=(mpc.bus(branch.ind_bus_F,9)-mpc.bus(branch.ind_bus_T,9))*pi/180+dif_theta;
% bus.Vmin=mpc.bus(:,8)-dif_V;
% bus.Vmax=mpc.bus(:,8)+dif_V;
%  bus.theta=mpc.bus(:,9)*pi/180;
%  bus.V=mpc.bus(:,8);
% bus.Pmin=(Cg*mpc.gen(:,2)-mpc.bus(:,3))/mpc.baseMVA-dif_PQ;
% bus.Pmax=(Cg*mpc.gen(:,2)-mpc.bus(:,3))/mpc.baseMVA+dif_PQ;
% bus.Qmin=(Cg*mpc.gen(:,3)-mpc.bus(:,4))/mpc.baseMVA-dif_PQ;
% bus.Qmax=(Cg*mpc.gen(:,3)-mpc.bus(:,4))/mpc.baseMVA+dif_PQ;