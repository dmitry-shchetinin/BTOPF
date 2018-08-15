function [ lb, ub ] = BTpostprocess_bounds( lb, ub, uniq )
%when there is more than one branch between two buses, finds the tightest
%bounds from distributes them to all parallel branches

%number of unique connections
L_uniq=uniq.L_uniq;

%iterate over unique connections
counter=0;
for i=1:L_uniq
    %get indices of parallel branches for particular connection
    ind_temp=uniq.ordering(counter+1:counter+uniq.n_entries(i));
    
    %get bounds and connection direction
    lb_temp=lb(ind_temp);
    ub_temp=ub(ind_temp);
    signs=uniq.signs(ind_temp);
    negative_signs=signs==-1;
    
    %for opposite direction of connection, swap lower and upper bounds    
    b_temp=-ub_temp(negative_signs);
    ub_temp(negative_signs)=-lb_temp(negative_signs);
    lb_temp(negative_signs)=b_temp;
    
    %find tightest lower and upper bounds
    minbound=max(lb_temp);
    maxbound=min(ub_temp);
    
    %record tightest bounds for all parallel branches
    lb(ind_temp(~negative_signs))=minbound;
    lb(ind_temp(negative_signs))=-maxbound;
    ub(ind_temp(~negative_signs))=maxbound;
    ub(ind_temp(negative_signs))=-minbound;
    
    %update counter
    counter=counter+uniq.n_entries(i);
end

end

