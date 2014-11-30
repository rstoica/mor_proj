function [ Y, T, X, YSD ] = step_sys( sys , t )
% step_sys Wrapper for the step function from MATLAB dynamical systems
% toolbox which needs to be installed as dependency. For more details, type
% "help step"

    sys_dummy = ss(sys.A, sys.B, sys.C, sys.D);
    switch nargin
        case 2
            [Y, T, X, YSD] = step(sys_dummy, t);
        case 1
            [Y, T, X, YSD] = step(sys_dummy);
        otherwise
            assert(0);
    end
end

