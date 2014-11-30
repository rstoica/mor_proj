function [ sys ] = anal_sys( sys )
% anal_sys  Analyzes a given system by checks of stability, controllability
% and observability.
%   
%   INPUTS: sys = system to be analyzed
%   
%   OUTPUTS: sys = system analyzed (essentailly the same as before), but
%   with added struct properties poles, R, O, which represent the poles of
%   the system, and respectively the controllability and the observability
%   matrices.

    % compute poles - just the eigenvalues of A
    if ~(isfield(sys,'poles'))
        [~,sys.poles] = eig(sys.A);
        sys.poles = diag(sys.poles);
    end
    if (sum(abs(real(sys.poles))>=10^-12)==length(sys.poles))
        fprintf('System "%s" is stable\n',sys.name);
    else
        fprintf('System "%s" is not stable\n',sys.name);
    end

    % determine whether is controllable
    R = sys.B;
    next = R;
    for i=1:(sys.n-1)
        next = sys.A*next; 
        R = [R,next];
    end
    sys.R = R;
    [~,SR,~] = svd(R);
    sR = diag(SR);
    if (sum(abs(sR)>=10^-12) ~= min(size(R)))
        fprintf('System "%s" is not controllable\n',sys.name);
    else
        fprintf('System "%s" is controllable\n',sys.name);
    end

    % determine whether is observable
    O = sys.C;
    next = O;
    for i = 1:(sys.n-1)
        next = next*sys.A;
        O = [O;next];
    end
    sys.O = O;
    [~, SO, ~] = svd(O);
    sO = diag(SO);
    if (sum(abs(sO)>=10^-12) ~= min(size(O)))
        fprintf('System "%s" is not observable\n',sys.name);
    else
        fprintf('System "%s" is observable\n',sys.name);
    end
end
