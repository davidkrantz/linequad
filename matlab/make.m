% make.m – Cross-platform MEX compilation script for MATLAB

% Detect platform
is_mac = ismac;
is_linux = isunix && ~ismac;

% Compiler and flags
openmp = true;     % Default to true, will be disabled on macOS
verbflag = false;  % Set to true to enable verbose compilation

my_dir = pwd();
c_dir = [pwd '/../c'];
build_dir = [pwd '/../build'];

% -Ofast makes code 2x faster than -O3, but causes the results of the unstable parts of
% the computation to differ between MEX and MATLAB (i.e. weights will differ, though
% quadrature results will not)
if is_mac
    cc = 'clang';
    cflags = ['-std=c99 -fPIC -I' c_dir];
    coptimflags = '-O3';
    ldoptimflags = '-O3';
    openmp = false;  % Clang on macOS doesn't support OpenMP by default
elseif is_linux
    cc = 'gcc';
    cflags = ['-std=c99 -fPIC -march=native -fopt-info-vec -I' c_dir];
    coptimflags = '-Ofast';
    ldoptimflags = '-Ofast';
else
    error('Unsupported platform');
end

% Add OpenMP flags if applicable
if openmp
    coptimflags = [coptimflags ' -fopenmp '];
    ldoptimflags = [ldoptimflags ' -fopenmp '];
end

ldflags = ' -lm ';

% Compiler string for MEX
CC = ['CC=''' cc ''''];
LD = '';  % Empty, not used here
CFLAGS = [' CFLAGS=''' cflags ''''];
LDFLAGS = [' LDFLAGS="\$LDFLAGS ' ldflags '" '];
OPTIMFLAGS = [' COPTIMFLAGS=''' coptimflags '''' ' LDOPTIMFLAGS=''' ldoptimflags ''''];
VERBOSE = '';

if verbflag
    VERBOSE = ' -DVERBOSE ';
end

% Build MEX command string
mex_string = ['mex ' CC LD CFLAGS LDFLAGS OPTIMFLAGS VERBOSE ' -outdir bin/'];

% Compile C object file
c_files = [' ' c_dir '/linequad.c '];
obj_file = [' ' build_dir '/linequad.o'];
c_call = [cc ' ' cflags ' ' coptimflags ' -c ' c_files '-o' obj_file ' ' ldflags];
disp(c_call)
assert(system(c_call) == 0)

% Compile MEX files
eval([mex_string obj_file ' mex/rootfinder_initial_guess_mex.c -output rootfinder_initial_guess_mex'])
eval([mex_string obj_file ' mex/rootfinder_mex.c -output rootfinder_mex'])
eval([mex_string obj_file ' mex/rsqrt_pow_weights_mex.c -output rsqrt_pow_weights_mex'])
eval([mex_string obj_file ' mex/bclag_interp_matrix_mex.c -output bclag_interp_matrix_mex'])
