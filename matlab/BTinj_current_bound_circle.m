function I_bound = BTinj_current_bound_circle( VVmin, VVmax, Pmin, Pmax, Qmin, Qmax, g_shunt, b_shunt )
%returns bound on magnitude of bus current injections. This is temporary 
%solution as for generator buses it uses samplings  
    
N=numel(VVmin);
Npoints=20; %number of samples point along each dimension (V, P, Q)
    
I_bound=zeros(N,1);
for i=1:N
    if (Pmin(i)~=Pmax(i) || Qmin(i)~=Qmax(i)) %if injections are not fixed
        %create vectors of sample points along each dimension
        Pall=create_datapoints(Pmin(i),Pmax(i),Npoints);
        Qall=create_datapoints(Qmin(i),Qmax(i),Npoints);
        Vall=create_datapoints(VVmin(i),VVmax(i),Npoints);

        %create one long vector of datapoints for all three variables
        Pall=repmat(Pall,Npoints^2,1);
        Qall=reshape(repmat(Qall,Npoints,Npoints)', Npoints^3, 1);
        Vall=reshape(repmat(Vall,1,Npoints^2)', Npoints^3, 1);

        %compute maximum current among all sample points
        I_bound(i)=sqrt(max((Pall.^2+Qall.^2)./Vall+...
            Vall*(g_shunt(i)^2+b_shunt(i)^2)+2*(Pall*g_shunt(i)-Qall*b_shunt(i))));
    else %if both injections are fixed
        I_bound(i)=(VVmax(i)-VVmin(i))/(VVmin(i)*VVmax(i))*sqrt(Pmin(i)^2+Qmin(i)^2);
    end
end
    
    
end


function points=create_datapoints(lb,ub,N)
    if (lb==ub)
        points=lb*ones(N,1);
    else
        points=(lb:(ub-lb)/(N-1):ub+1e-10)';
    end
end
