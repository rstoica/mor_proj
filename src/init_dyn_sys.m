function [ sys ] = init_dyn_sys( m, n, p, name )
% init_dyn_sys Initializes a continuous dynamic system
%   
%   INPUTS: m = number of inputs that the system takes
%           n = number of states that the system has
%           p = number of outputs the system gives
%           name = string as name of the system
%
%   OUTPUTS: sys = MATLAB structure containing A, B, C, D empty matrices
%   initialized with zeros based on the system dimensionality
    sys.m = m;   % number of inputs
    sys.n = n;   % number of states
    sys.p = p;   % number of outputs
    sys.name = name;
    sys.A = zeros(n,n);
    sys.B = zeros(n,m);
    sys.C = zeros(p,n);
    sys.D = zeros(p,m);
end

