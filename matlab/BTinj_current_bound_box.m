function [Imin,Imax] = BTinj_current_bound_box(Smin,Smax,Vmin,Vmax,alpha)
%returns min and max values of injection current considering shunt. This is
%an estimate because the values are computed through sampling

Npoints=20;

%create vectors of sample points along each dimension
S=create_datapoints(Smin,Smax,Npoints);
V=create_datapoints(Vmin,Vmax,Npoints);

%create one long vector of datapoints for all two variables
S_all=repmat(S,Npoints,1);
V_all=reshape(repmat(V,1,Npoints)', Npoints^2, 1);

%compute I_max at all points
f=S_all./V_all-alpha*V_all;
Imax=max(f);
Imin=min(f);
end


function points=create_datapoints(lb,ub,N)
    if (lb==ub)
        points=lb*ones(N,1);
    else
        points=(lb:(ub-lb)/(N-1):ub+1e-10)';
    end
end

