function [ SV, W ] = sigma_sys( sys, w )
% sigma_sys Wrapper for sigma MATLAB function from dependency toolbox
% Dynamical Systems. This is a light weight implementation, allowing just
% the sys parameter to be passed. This yields the Bode plot of the SV of 
% the system. For more details type "help sigma"

    sys_dummy = ss(sys.A,sys.B,sys.C,sys.D);
    switch nargin
        case 2
            [SV, W] = sigma(sys_dummy,w);
        case 1
            [SV, W] = sigma(sys_dummy);
        otherwise
            assert(0);
    end
end

