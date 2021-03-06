%example showing how to run bound tightening algorithm
clear;

%load test case
casename='pglib_opf_case2736sp_k'; %case300, case1354pegase, case2383wp, case3120sp
mpc=loadcase(casename);

%set some options to non-default values
options.bounds=1; %tighten both voltage and angle differences
%options.opt_prob_flow=2; %use constraint in LP for BT based on flow envelopes
%options.opt_prob_inj=3; %use LP+SOCP for BT based on injection envelopes

%run bound tightening (refer to function description for info on inputs/outputs)
[ mpc, info ] = tighten_bounds( mpc, options);