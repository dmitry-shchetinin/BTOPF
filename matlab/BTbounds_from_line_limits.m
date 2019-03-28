function [ theta, Vdif ] = BTbounds_from_line_limits( bus, branch, opt )
%Returns bounds on angle and voltage differences for all lines based on the
%feasible set of thermal limit constraints. If a line does not have a
%thermal limit constraint (modeled as branch.I_max=0), the limits on phase
%andgle differences are +-pi/2, and limits on voltage differences are
%taken from the bounds on voltages
    
%% initialize output structures
L=numel(branch.g);
theta=struct('min',-pi*ones(L,1)/2,'max',pi*ones(L,1)/2);
pow=opt.Vdif_type;
Vdif=struct('min',bus.Vmin(branch.ind_bus_F).^pow-bus.Vmax(branch.ind_bus_T).^pow,...
    'max',bus.Vmax(branch.ind_bus_F).^pow-bus.Vmin(branch.ind_bus_T).^pow);
extra=struct('slope1',zeros(L,1),'offset1',zeros(L,1),'slope2',zeros(L,1),'offset2',zeros(L,1));


%iterate over all branches
for i=1:L
    if (branch.I_max(i)>0)
        ratio=branch.ratio(i);

        %% obtain limits on V_i and V_j (after transformer)
        Vbox=struct('V_i_min',bus.Vmin(branch.ind_bus_F(i))/ratio,...
            'V_j_min',bus.Vmin(branch.ind_bus_T(i)),...
            'V_i_max',bus.Vmax(branch.ind_bus_F(i))/ratio,...
            'V_j_max',bus.Vmax(branch.ind_bus_T(i)));

        %% deal with the beginning of the line
        branch_params=BTbranch_parameters(branch.g(i), branch.b(i), branch.bs(i), ratio, branch.I_max(i), 1);
        [ theta_begin, Vdif_begin ] = BTbounds_from_one_line_end(branch_params, Vbox, ratio);

        %% deal with the end of the line
        branch_params=BTbranch_parameters(branch.g(i), branch.b(i), branch.bs(i), ratio, branch.I_max(i), 2);
        [ theta_end, Vdif_end ] = BTbounds_from_one_line_end(branch_params, Vbox, ratio);

        %% record final angle limits
        theta.min(i)=min(theta_begin.min,theta_end.min)+branch.shift(i);
        theta.max(i)=max(theta_begin.max,theta_end.max)+branch.shift(i);
        
        
        %% record voltage limits in form |Vj-alpha*Vi|<beta
        [extra.slope1(i), extra.offset1(i), extra.slope2(i), extra.offset2(i)]=deal(Vdif_begin.slope, ...
            Vdif_begin.offset, Vdif_end.slope, Vdif_end.offset);
        
        %% do postprocessing to obtain |Vi-Vj| constraints
        Vbox.V_i_min=Vbox.V_i_min*ratio;
        Vbox.V_i_max=Vbox.V_i_max*ratio;
        [Vdif.min(i), Vdif.max(i)]=compute_bound_Vdif(Vbox, Vdif_begin, Vdif_end);
    end
end

%% record all obtained voltage constraints
Vdif.extra=extra;

%% obtain differences of squares of voltage magnitudes if need be
if (pow==2)
    ind_Imax=branch.I_max>0;
    ind_temp=Vdif.min<0;
    Vdif.min(ind_temp & ind_Imax)=Vdif.min(ind_temp & ind_Imax).*(bus.Vmax(branch.ind_bus_F(ind_temp & ind_Imax))+bus.Vmax(branch.ind_bus_T(ind_temp & ind_Imax)));
    Vdif.min(~ind_temp & ind_Imax)=Vdif.min(~ind_temp & ind_Imax).*(bus.Vmin(branch.ind_bus_F(~ind_temp & ind_Imax))+bus.Vmin(branch.ind_bus_T(~ind_temp & ind_Imax)));
    ind_temp=Vdif.max>=0;
    Vdif.max(ind_temp & ind_Imax)=Vdif.max(ind_temp & ind_Imax).*(bus.Vmax(branch.ind_bus_F(ind_temp & ind_Imax))+bus.Vmax(branch.ind_bus_T(ind_temp & ind_Imax)));
    Vdif.max(~ind_temp & ind_Imax)=Vdif.max(~ind_temp & ind_Imax).*(bus.Vmin(branch.ind_bus_F(~ind_temp & ind_Imax))+bus.Vmin(branch.ind_bus_T(~ind_temp & ind_Imax)));
end

%% postprocess bounds (to store tightest bounds for all parallel branches)
[ theta.min, theta.max ] = BTpostprocess_bounds( theta.min, theta.max, branch.uniq );
[ Vdif.min, Vdif.max ] = BTpostprocess_bounds( Vdif.min, Vdif.max, branch.uniq );   
end


%% returns bound on |Vi-Vj| from Vj-slope*Vi=offset constraint
function [Vdif_min, Vdif_max]=compute_bound_Vdif(Vbox, Vdif_begin, Vdif_end)
Vdifs_upper=zeros(4,1);
Vdifs_lower=zeros(4,1);
%Vdif resulting from intersections of lines at the beginning of the branch
[Vdifs_upper(1),Vdifs_upper(2)]=compute_Vdif(Vbox, Vdif_begin.slope, -Vdif_begin.offset);
[Vdifs_lower(1),Vdifs_lower(2)]=compute_Vdif(Vbox, Vdif_begin.slope, Vdif_begin.offset);
%Vdif resulting from intersections of lines at the end of the branch
[Vdifs_upper(3),Vdifs_upper(4)]=compute_Vdif(Vbox, Vdif_end.slope, -Vdif_end.offset);
[Vdifs_lower(3),Vdifs_lower(4)]=compute_Vdif(Vbox, Vdif_end.slope, Vdif_end.offset);
%compute maximum Vdif
Vdifs_upper(Vdifs_upper==inf)=[];
if (~isempty(Vdifs_upper))
    Vdif_max=max(Vdifs_upper);
else
    Vdif_max=Vbox.V_i_max-Vbox.V_j_min;
end
%compute minimum Vdif
Vdifs_lower(Vdifs_lower==inf)=[];
if (~isempty(Vdifs_lower))
    Vdif_min=min(Vdifs_lower);
else
    Vdif_min=Vbox.V_i_min-Vbox.V_j_max;
end
end


%% returns Vdif for given line at two intersection points with Vbox
function [Vdif_begin, Vdif_end]=compute_Vdif(Vbox, slope, offset)
%check if this line has intersections with the box
if ((slope*Vbox.V_i_min+offset>Vbox.V_j_max) || (slope*Vbox.V_i_max+offset<Vbox.V_j_min))
    Vdif_begin=inf;
    Vdif_end=inf;
    return;
end
%record Vdif at the beginning of the line
if (slope*Vbox.V_i_min+offset>=Vbox.V_j_min)
    Vi=Vbox.V_i_min;
    Vj=slope*Vi+offset;
else
    Vj=Vbox.V_j_min;
    Vi=(Vj-offset)/slope;
end
Vdif_begin=Vi-Vj;
%record Vdif at the end of the line
if (slope*Vbox.V_i_max+offset<Vbox.V_j_max)
    Vi=Vbox.V_i_max;
    Vj=slope*Vi+offset;
else
    Vj=Vbox.V_j_max;
    Vi=(Vj-offset)/slope;
end
Vdif_end=Vi-Vj;
end

