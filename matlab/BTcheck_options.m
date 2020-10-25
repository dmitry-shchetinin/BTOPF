function [ opt ] = BTcheck_options( opt )
%check user provided options and set not provided options to their defaults
%note: only relevant options are checked, i.e. if a particular option is
%not applicable to a chosen method, its value is not checked.

%% types of bounds to tighten
opt=add_field(opt,'bounds', 1);
if (opt.bounds~=1 && opt.bounds~=2 && opt.bounds~=3)
    warning('Value of field "bounds" must be 1, 2 or 3. Changing it to default value of 1.')
    opt.bounds=1;
end

%% tightening method to use
opt=add_field(opt,'method', 1);
if (opt.method~=1 && opt.method~=2 && opt.method~=3 && opt.method~=4)
    warning('Value of field "method" must be 1, 2, 3, or 4. Changing it to default value of 1.')
    opt.method=1;
elseif (opt.method==4 && opt.bounds==2)
    warning('Value of field "method" cannot be 3 when "bounds"=2. Changing it to default value of 1.')
    opt.method=1;
end

%% type of voltage differences to tighten
opt=add_field(opt,'Vdif_type', 1);
if (opt.Vdif_type~=1 && opt.Vdif_type~=2)
    warning('Value of field "Vdif_type" must be 1 or 2. Changing it to default value of 1.')
    opt.Vdif_type=1;
end

%% linear solver to use
opt=add_field(opt,'linsolver', 1);
if (opt.linsolver~=1 && opt.linsolver~=2)
    warning('Value of field "linsolver" must be 1 or 2. Changing it to default value of 1.')
    opt.linsolver=1;
elseif (opt.linsolver==2 && exist('klu','file')~=3)
    warning('KLU function was not found. Changing "linsolver" to default value of 1.')
    opt.linsolver=1;
end

%% usage of mex functions to speed up computations
opt=add_field(opt,'use_mex', 0);
if (opt.use_mex~=0 && opt.use_mex~=1)
    warning('Value of field "use_mex" must be 0 or 1. Changing it to default value of 0.')
    opt.use_mex=1;
elseif (opt.use_mex==1 && exist('greedy_LP_solution_mex','file')~=3)
    warning('Required mex files were not found. Changing "use_mex" to default value of 1.')
    opt.use_mex=1;
end

%% number of iterations for system-wide tightening methods
opt=add_field(opt,'max_iter', 5);
if (opt.max_iter<1 && opt.method~=2)
    warning('Value of field "max_iter" must be at least 1. Changing it to default value of 5.')
    opt.max_iter=1;
elseif (opt.max_iter>10)
    warning('Value of field "max_iter" seems too high. Consider changing it to at most 10.')
end

%% stopping criterion: max change in bounds on two consecutive iterations
opt=add_field(opt,'max_change', 1e-3);
if (opt.max_change<=0 && opt.method~=2 && opt.bounds~=2)
    warning('Value of field "max_change" must be positive. Changing it to default value of 1e-3.')
    opt.max_change=1e-3;
end

%% type of optimization problems solved when tightening is based on envelopes of power flows
opt=add_field(opt,'opt_prob_flow', 1);
if (opt.opt_prob_flow~=1 && opt.opt_prob_flow~=2 && opt.method~=2 && opt.method~=4)
    warning('Value of field "opt_prob_flow" must be 1 or 2. Changing it to default value of 1.')
    opt.opt_prob_flow=1;
end

%% type of optimization problems solved when tightening is based on envelopes of bus injections
opt=add_field(opt,'opt_prob_inj', 1);
if (opt.opt_prob_inj~=1 && opt.opt_prob_inj~=2 && opt.opt_prob_inj~=3 && ...
        opt.bounds~=2 && opt.method~=2 && opt.method~=3)
    warning('Value of field "opt_prob_inj" must be 1, 2, or 3. Changing it to default value of 1.')
    opt.opt_prob_inj=1;
end

%% what information to record
opt=add_field(opt,'statistics', 1);
if (opt.statistics~=0 && opt.statistics~=1 && opt.statistics~=2)
    warning('Value of field "statistics" must be 0, 1, or 2. Changing it to default value of 1.')
    opt.statistics=1;
end

%% order fields by their names
opt=orderfields(opt);
end


function struct=add_field(struct, fieldname, fieldvalue)
    if ~isfield(struct,fieldname)
        struct.(fieldname)=fieldvalue;
    end
end