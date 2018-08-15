function [obj, info] = greedy_LP_solution( c, a, xl, xu, bl, bu )
%this function solves an LP with one linear constraint by a greedy algorithm:
% max c'x
% s.t. xl<= x <=xu
%      bl<=a'x<=bu
%output info shows: 0 - constraint not binding, 1 - lower bound of
%constraint is binding, 2 - upper bound of constraint is binding

eps=1e-13; %define threshold for very small values

%change variables such that xl=0
xu=xu-xl;

temp=a*xl;
bl=bl-temp;
bu=bu-temp;
obj=c*xl;

%handle the cases with zero or very small a
ind_zero=abs(a)<=eps; num_small_a=sum(ind_zero);
if (num_small_a>0)
    %add elements with positive c to the objective
    obj=obj+c(c>0 & ind_zero)*xu(c>0 & ind_zero);
    if (num_small_a==numel(xl))
        return; %if constraint is essentially non-existent, exit
    end
end

%get indices of elements in x s.t. c>0, a>0
ind_1=c>0 & a>eps;

%get indices of elements in x s.t. c>0, a<0
ind_2=c>0 & a<-eps;

%get indices of elements in x s.t. c>0, a>0 or c>0, a<0
ind_12=ind_1 | ind_2;

%get the maximum possible value of the objective and corressponding value of constraint
temp=xu(ind_12);
sum_a=a(ind_12)*temp;
sum_c=c(ind_12)*temp;
%check if entries with positive c exist
if isempty(sum_a)
    sum_c=0;
    sum_a=0;
end
    
%check if a'x constraint is binding at the optimum solution
if (sum_a>=bl && sum_a<=bu)
    obj=obj+sum_c;
    info=0;
    return;
elseif sum_a>bu % the value is too big
    %we need to either remove some x with a>0,c>0 or add some x with a<0,c<=0 until a'x=bu
    ind_change=ind_1 | (c<=0 & a<-eps);
    %compute the violation of the constraint
    b_viol=sum_a-bu;
    info=2;
else %the value is too small
    %we need to either remove some x with a<0,c>0 or add some x with a>0,c<=0 until a'x=bl
    ind_change=ind_2 | (c<=0 & a>eps);
    b_viol=bl-sum_a;
    info=1;
end

%extract the indexes
ind_change=find(ind_change);

%check if the problem is infeasible
if (isempty(ind_change))
    obj=-inf;
    return;
end

%sort the indices of elements to change in the increasing order of |a/c|
[~,ind_reorder]=sort(abs(c(ind_change)./a(ind_change)));
ind_change=ind_change(ind_reorder);

%get number of elements that can be potentially changed
n=numel(ind_change);

%add and/or remove elements by a greedy strategy until the constraint is satisfied
obj=obj+sum_c;
for i=1:n
    index=ind_change(i);
    %check if we need to remove the entire x entry
    a_temp=abs(a(index)*xu(index));
    if (b_viol>=a_temp)
        b_viol=b_viol-a_temp;
        obj=obj-abs(c(index)*xu(index));
    else
        obj=obj-abs(c(index)*b_viol/a(index));
        b_viol=0;
        break;
    end
end

%check if the problem is feasible
if (abs(b_viol)>1e-6)
    obj=-inf;
end

end