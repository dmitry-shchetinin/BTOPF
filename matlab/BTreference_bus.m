function [ ind_ref ] = BTreference_bus( mpc )
%returns the index of the bus whose angle is set as reference.
%(preferably choose terminal node with smallest nonzero load)
    
N=size(mpc.bus,1);
L=size(mpc.branch,1);

load_temp=sqrt((mpc.bus(:,3)).^2+(mpc.bus(:,4)).^2);
ind_ref=find(load_temp==min(load_temp(load_temp>0)),1); %smallest nonzero load
%consider only load buses
load_temp1=load_temp(load_temp>0 & mpc.bus(:,2)==1);
if (~isempty(load_temp1))
    ind_ref=find(load_temp==min(load_temp1),1); %smallest nonzero load among load buses
end
%consider only load buses with one connection
Adj=sparse(mpc.branch(:,1), mpc.branch(:,2), ones(L,1),N,N);
Adj=Adj+Adj';
Adj(Adj>1)=1;
number_connected_nodes=Adj*ones(N,1);
load_temp1=load_temp(load_temp>0 & mpc.bus(:,2)==1 & number_connected_nodes==1);
if (~isempty(load_temp1))
    ind_ref=find(load_temp==min(load_temp1),1); %smallest nonzero load among load buses
end
end

