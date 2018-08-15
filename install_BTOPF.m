%install bound tightening package
clear;

%% change directory to folder containing this script
cd(fileparts(which(mfilename)));


%% add necessary folders to path
BTroot = pwd;
add_to_path([BTroot '/matlab']);
add_to_path([BTroot '/src']);
add_to_path([BTroot '/tests']);
cd(BTroot)


%% check is prebuilt mex files exist
ext = ['.',mexext]; %get extension of mex files for current platform
mexfiles={'greedy_LP_solution_mex'}; %list of mex files that BTOPF can use
fprintf ('\nChecking for mex functions......\n');
cd([BTroot '/src'])
for i=1:length(mexfiles)
    filename=mexfiles{i};
    
    %check if mex file already exists
    if (exist([filename,ext],'file')==3)
        mex_exists=1;
    else
        mex_exists=0;
    end
   
    if (exist(filename,'file')==3)
        fprintf('\nMex function "%s" already exists.\n',mexfiles{i});
        UserInput = input('Hit enter to recompile it (or "n" to use existing function): ', 's');
        if (isempty (UserInput))
            compile_mex_file(filename, ext, [BTroot '/matlab'], 1);
        end
    else
        compile_mex_file(filename, ext, [BTroot '/matlab'], 0);
    end
    
end
cd(BTroot)


%% check if MATPOWER function "makeYbus" is in path (BTOPF only uses this function from MATPOWER)
if (exist('makeYbus','file')~=2)
    fprintf('\n');
    warning('Matpower function "makeYbus" was not found in Matlab path. BTOPF cannot run without it.');
    fprintf('\n');
end


% %todos:
% check installation of matpower -in tests
% check if KLU is installed - in tests
% check if it is possible to compile c files locally
% compile c files if possible, otherwise use prebuild binaries (depending on operating system)
% offer to run the test


%% remind user to save path
fprintf ('\nBTOPF installation is complete. To use BTOPF in future sessions,\n');
fprintf ('save the following folder in your Matlab path:\n');
fprintf('"%s"\n',[BTroot '\matlab']);
fprintf('by using "savepath" command. Refer to documentation if you cannot\n');
fprintf('permanently save this path because of file permissions.\n\n');


%% run tests if need be
UserInput=input('Hit enter to run tests (or "n" to quit): ', 's');
if (isempty(UserInput))
    BT_test_installation;
end


%% functions for installation
function add_to_path(newpath)
    cd (newpath) ;
    addpath (newpath) ;
end

function compile_mex_file(filename, ext, destination, mex_exists)
    fprintf('\nCompiling mex function "%s".........\n',filename);
    try
        mex([filename, '.cpp']);
        try
            movefile([filename,ext],destination);
        catch
            fprintf('Matlab was unable to move created mex function to folder "%s"\nPlease move it there manually.\n',destination);
        end
    catch ME
        fprintf('Function "%s" was not compiled because of the following error: \n',filename);
        fprintf('  "%s"\n',ME.identifier);
        if (mex_exists)
            fprintf('Existing mex function will be used.\n');
        else
            fprintf('Matlab functions will be used instead in the algorithm.\n');
        end
    end
end
