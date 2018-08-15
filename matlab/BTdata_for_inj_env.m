function [ bus, branch ] = BTdata_for_inj_env( bus, branch )
%returns data required for running BT based on convex envelopes of bus
%current injections
    
%unpack data
[N, L_uniq, uniq]=deal(bus.nbus, branch.uniq.L_uniq, branch.uniq);

%compute number of adjacent buses for each bus
Adj=sparse(uniq.ind_F, uniq.ind_T, ones(L_uniq,1),N,N); %half of adjacency matrix
Adj=Adj+Adj'; %full adjacency matrix (buses are not connected to themselves, i.e. main diagonal is all zeros)
n_adj_buses=Adj*ones(N,1); %compute number of adjacent buses for each node
%sort the buses in the descending order of number of connections
[n_adj_buses,ordering]=sort(n_adj_buses,'descend');

%reorder bus numbers
bus_number=(1:N)';
bus_number=bus_number(ordering);

%apply consecutive bus numbering
e2i=sparse(bus_number,ones(N,1),1:N, N, 1);
busF=full(e2i(branch.ind_bus_F));
busT=full(e2i(branch.ind_bus_T));

%get unique connections between reordered buses
[~,ind_temp,~]=unique([busF.*busT, busF+busT], 'rows','stable');
busF=busF(ind_temp);
busT=busT(ind_temp);

%obtain adjacent buses for each bus in the list of reordered buses
Adj=sparse(busF, busT, ones(L_uniq,1), N, N);
Adj=Adj+Adj';
[temp,~,~]=find(Adj);
adj_buses=cell(N,1);
counter=0;
for i=1:N
    adj_buses{i}=temp(counter+1:counter+n_adj_buses(i));
    counter=counter+n_adj_buses(i);
end
bus.adj_buses_all=adj_buses;

%obtain adjacent branches for each bus
Adj=sparse(busF, busT, 1:L_uniq, N, N); %half of adjacency matrix
Adj=Adj+Adj'; %full adjacency matrix (buses are not connected to themselves, i.e. main diagonal is all zeros)
[~,~,temp]=find(Adj);
adj_branches=cell(N,1);
counter=0;
for i=1:N
    adj_branches{i}=temp(counter+1:counter+n_adj_buses(i));
    counter=counter+n_adj_buses(i);
end

%make sure each connection is only touched once
temp=zeros(L_uniq,1);
for i=1:N
    %check which connections were already handled
    ind_repeated=temp(adj_branches{i})==1;
    
    %delete information about these connections
    adj_buses{i}(ind_repeated)=[];
    adj_branches{i}(ind_repeated)=[];
    n_adj_buses(i)=n_adj_buses(i)-sum(ind_repeated);
    
    %record connections handled by this bus
    temp(adj_branches{i})=1;
end

%for each bus, record whether it is From bus or To bus for adjacent branches
adj_signs=cell(N,1);
for i=1:N
    adj_signs{i}=ones(n_adj_buses(i),1);
    adj_signs{i}(busT(adj_branches{i})==i)=-1;
end

%store results
bus.adj_buses=adj_buses;
bus.adj_branches=adj_branches;
bus.adj_signs=adj_signs;
bus.n_adj_buses=n_adj_buses;
bus.ordering=ordering;
end

