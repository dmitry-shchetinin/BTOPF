function [ mpc, info ] = tighten_bounds( mpc, options)
%% FUNCTION DESCRIPTION
%tighten bounds on angles and/or voltage magnitude differences
%INPUTS
%  mpc     - structure containing system data in Matpower format
%  options - structure containing user defined options. This argument is 
%            optional. If it is not provided or is empty, all fields are 
%            initialized to their defaults. If some fields are not provided 
%            by the user, they will also be initizalized to their defaults.
%            The fields are (default values are given in brackets):
%    .bounds         (1)   type of bounds to tighten for each branch:
%                           1 - phase angle difference
%                           2 - voltage magnitude difference. Tightening
%                               only Vdif is not recommended unless tight 
%                               bounds of angles are available. Therefore, 
%                               if Vdif is needed, bounds=3 is recommended. 
%                           3 - both
%    .method         (1)   method to be used for tightening:
%                           1 - all applicable methods
%                           2 - based on line flow constraints
%                           3 - based on envelopes of power flows
%                           4 - based on envelopes of bus injections (only 
%                               applicable for tightening angle bounds)
%    .Vdif_type      (1)   type of voltage difference bounds to tighten:
%                           1 - difference of voltage magnitudes
%                           2 - difference of squares of voltage magnitudes
%    .linsolver      (1)   linear solver to use:
%                           1 - default Matlab LU
%                           2 - KLU from SuiteSparse (installed separately)
%    .use_mex        (0)   usage of mex files to accelerate computations:
%                           0 - do not use mex files
%                           1 - use mex files (must compile them locally or 
%                               download prebuilt binaries into path)
%    .max_iter       (5)   number of iterations for tightening methods that
%                          are based on convex envelopes
%    .max_change   (1e-3)  stopping condition for iterative tightening
%                          methods: maximum change in bounds on two 
%                          consecutive iterations
%    .opt_prob_flow  (1)   type of optimization problem solved for bound 
%                          tightening based on envelopes of flows:
%                           1 - LP with only variable bounds (faster but 
%                               can produce looser bounds)
%                           2 - LP with variable bounds and one constraint 
%                               (slower but can produce tighter bounds)
%    .opt_prob_inj   (1)   type of optimization problem solved for bound 
%                          tightening based on envelopes of injections:
%                           1 - SOCP with simple cone constraints (fastest)
%                           2 - LP with bounds on variables and one linear
%                               constraint (slower)
%                           3 - SOCP+LP (slowest but yields tighter bounds)
%    .statistics     (1)   statistical information collected:
%                           0 - nothing
%                           1 - only general info
%                           2 - detailed info for each method
%%OUTPUTS
%  mpc     - structure containing system data in Matpower format. If bounds
%            on angle differences are tightened, the updated bounds are
%            stored in the corresponding columns of mpc.branch matrix. If 
%            bounds on voltage magnitude differences were tightened, extra 
%            field is added to structure:
%    .Vdif         - it is itself a structure with the following fields:
%       .min                lower bound on Vi-Vj for all branches
%       .max                upper bound on Vi-Vj for all branches
%       .extra (optional)   stores parameters of constraints of the form
%                           |Vj-slope*Vi|<offset, which are obtained from
%                           the feasible set of thermal limit constraint
%            In addition, information on voltage magnitudes differences is
%            added to columns of matrix mpc.branch. Column indices and the
%            corresponding information stored in them is as follows:
%       22     - lower bound on Vi-Vj for all branches
%       23     - upper bound on Vi-Vj for all branches
%       24     - value of slope in |Vj-slope*Vi|<offset for thermal limit
%                constraint at the beginning of the line. If thermal limit 
%                value is not given, slope=0.
%       25     - value of offset in |Vj-slope*Vi|<offset for thermal limit
%                constraint at the beginning of the line. If thermal limit 
%                value is not given, offset=0.
%       26     - value of slope in |Vj-slope*Vi|<offset for thermal limit
%                constraint at the end of the line. If thermal limit value 
%                is not given, slope=0.
%       27     - value of offset in |Vj-slope*Vi|<offset for thermal limit
%                constraint at the end of the line. If thermal limit value 
%                is not given, offset=0.
%  info    - optional output containing statistical information about the
%            tightening results. If options.statistics==1, the following
%            fields in structure info are present:
%    .time            - total elapsed time (excluding preprocessing). if
%                       options.statistics==2, it is structure with fields
%                       containing runtimes of each employed method as well 
%                       as the total time
%    .stat_Vdif_all   - structure containing general statistics on the
%                       bound tightening for Vdif. The fields are:
%       .nb_p             - percentage of branches, for which at least one 
%                           bound was tightened
%       .range_max        - max distance between the bounds (in p.u.)
%       .range_mean       - mean distance between the bounds (in p.u.)
%       .range_median     - median distance between the bounds (in p.u.)
%       .range_min        - min distance between the bounds (in p.u.)
%       .drange_max_p     - max change in distance between the bounds, in
%                           percent (relative to the initial range)
%       .drange_mean_p    - mean change in distance between the bounds, in
%                           percent (relative to the initial range)
%    .stat_theta_all  - structure containing general statistics on the
%                       bound tightening for theta. Same fields as in
%                       'stat_Vdif_all' but distance between bounds is
%                       measured in radian
%          - if options.statistics==2, extra information is collected and
%            returned to the user. On top of fields 'stat_theta_all' and 
%            'stat_Vdif_all' similar structures are added for results of 
%            each of employed methods. These structures have the following
%            extra-fields:
%       .n_lb_p           - percent of branches with tightened lower bounds
%       .d_lb_max         - max tightening of lower bound
%       .d_lb_mean        - mean tightening of lower bound
%       .n_ub_p           - percent of branches with tightened upper bounds
%       .d_ub_max         - max tightening of upper bound
%       .d_ub_mean        - mean tightening of upper bound
%          - if options.statistics==2, another structure field 'viol' is 
%            added to 'info' for each bound tightening method and overall
%            result. It contains violations of the bounds by operating
%            point contained in input mpc structure. The fields are:
%    .viol
%       .n_lb             - number of violated lower bounds
%       .n_ub             - number of violated upper bounds
%       .lb_max           - maximum violation of lower bound
%       .lb_mean          - mean violation of lower bound
%       .ub_max           - maximum violation of upper bound
%       .ub_mean          - mean violation of upper bound


%% check input arguments
if nargin < 1
    error('You must provide at least case structure in matpower format.')
elseif nargin < 2
    options=[];
else
    if (~isempty(options) && ~isstruct(options))
        error('Second input argument must be either a stucture or empty.')
    end
end
info=[];

%% check user-provided options and set missing ones to default values
options=BTcheck_options(options);


%% preprocess input data
[bus, branch]=BTpreprocess_data(mpc, options);


%% record initial values of bounds
theta_initial=branch.theta;
pow=options.Vdif_type;
Vdif_initial=struct('min',bus.Vmin(branch.ind_bus_F).^pow-bus.Vmax(branch.ind_bus_T).^pow,...
    'max',bus.Vmax(branch.ind_bus_F).^pow-bus.Vmin(branch.ind_bus_T).^pow);
Vdif=Vdif_initial;
time=struct('flow_env',0,'inj_env',0,'thermal',0,'Vdif',0,'total',0);


%% tighten bounds based on feasible set of line flow constraints
if (options.method==1 || options.method==2)
    tic;
    [ theta_thermal, Vdif_thermal ] = BTbounds_from_line_limits( bus, branch, options );
    time.thermal=toc;
    %update bounds
    branch.theta=BTupdate_bounds(branch.theta, theta_thermal);
    if (options.bounds~=1)
        %record voltage difference info only if required
        Vdif=BTupdate_bounds(Vdif, Vdif_thermal);
        Vdif.extra=Vdif_thermal.extra;
    end
    %collect detailed statistics if need be
    if (options.statistics==2)
        [ info.stat_theta_thermal, info.viol_theta_thermal ] = BTcollect_statistics( theta_initial, ...
            branch.theta, branch.theta.val, 2, 1);
       [ info.stat_Vdif_thermal, info.viol_Vdif_thermal ] = BTcollect_statistics( Vdif_initial, ...
           Vdif, branch.Vdif, 2, 1);
    end
end


%% tighten bounds based on convex envelopes of bus injections
if ((options.method==1 || options.method==4) && options.bounds~=2)
    temp_bound=branch.theta;
    tic;
    theta_inj = BTbounds_from_inj_envelopes( bus, branch, options );
    time.inj_env=toc;
    %update bounds
    branch.theta=BTupdate_bounds(branch.theta, theta_inj);
    %collect detailed statistics if need be
    if (options.statistics==2)
        [ info.stat_theta_inj, info.viol_theta_inj ] = BTcollect_statistics( temp_bound, ...
            branch.theta, branch.theta.val, 2, 1);
    end
end



%% tighten angle bounds based on convex envelopes of power flows
if ((options.method==1 || options.method==3) && options.bounds~=2)
    temp_bound=branch.theta;
    tic;
    [ theta_flow, stat_iter ] = BTbounds_from_flow_envelopes( bus, branch, options, Vdif, 1 );
    time.flow_env=toc;
    %update bounds
    branch.theta=BTupdate_bounds(branch.theta, theta_flow);
    %collect detailed statistics if need be
    if (options.statistics==2)
        [ info.stat_theta_flow, info.viol_theta_flow ] = BTcollect_statistics( temp_bound, ...
            branch.theta, branch.theta.val, 2, 1);
        info.stat_theta_flow_detailed=stat_iter;
    end
end



%% tighten Vdif bounds based on convex envelopes of power flows
if ((options.method==1 || options.method==3) && options.bounds~=1)
    temp_bound=Vdif;
    tic;
    [ Vdif_flow, ~ ] = BTbounds_from_flow_envelopes( bus, branch, options, Vdif, 0 );
    time.Vdif=toc;
    %update bounds
    Vdif=BTupdate_bounds(Vdif, Vdif_flow);
    %collect detailed statistics if need be
    if (options.statistics==2)
        [ info.stat_Vdif_flow, info.viol_Vdif_flow ] = BTcollect_statistics( temp_bound, ...
            Vdif, branch.Vdif, 2, 1);
    end
end




%% record relevant statistics
time.total=time.flow_env+time.inj_env+time.thermal+time.Vdif;
if (options.statistics==1)
    [ info.stat_theta_all, ~ ] = BTcollect_statistics( theta_initial, ...
        branch.theta, [], 1, 0);
    [ info.stat_Vdif_all, ~ ] = BTcollect_statistics( Vdif_initial, Vdif, [], 1, 0);
    info.time=time.total;
elseif (options.statistics==2)
    [ info.stat_theta_all, info.viol_theta_all ] = BTcollect_statistics( theta_initial, ...
        branch.theta, branch.theta.val, 2, 1);
    [ info.stat_Vdif_all, info.viol_Vdif_all ] = BTcollect_statistics( Vdif_initial, ...
        Vdif, branch.Vdif, 2, 1);
    info.time=time;
end



%% save results
%angle differences
mpc.branch(branch.in_service,12)=branch.theta.min*180/pi;
mpc.branch(branch.in_service,13)=branch.theta.max*180/pi;
%voltage differences
Vbounds=struct('min',-inf(size(mpc.branch,1),1),'max',inf(size(mpc.branch,1),1));
Vbounds.max(branch.in_service)=Vdif.max;
Vbounds.min(branch.in_service)=Vdif.min;
if (isfield(Vdif,'extra'))
    extra=struct('slope1',ones(size(mpc.branch,1),1),'offset1',...
        inf(size(mpc.branch,1),1),'slope2',ones(size(mpc.branch,1),1),...
        'offset2',inf(size(mpc.branch,1),1));
    [extra.slope1(branch.in_service),extra.offset1(branch.in_service),...
        extra.slope2(branch.in_service),extra.offset2(branch.in_service)]=...
        deal(Vdif.extra.slope1,Vdif.extra.offset1,Vdif.extra.slope2,Vdif.extra.offset2);
    Vbounds.extra=extra;
end
mpc.Vdif=Vbounds;
%add voltage differences to mpc.branch
temp=zeros(size(mpc.branch,1),27);
temp(:,1:size(mpc.branch,2))=mpc.branch;
temp(:,22:23)=[Vbounds.min, Vbounds.max];
if (isfield(Vdif,'extra'))
    temp(:,24:27)=[extra.slope1,extra.offset1,extra.slope2,extra.offset2];
end
mpc.branch=temp;

%order fields for convenience
if (~isempty(info))
    info=orderfields(info);
end

end





