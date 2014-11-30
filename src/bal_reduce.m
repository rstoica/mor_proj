function [ sys_red, res ] = bal_reduce( sys, tol )
% bal_reduce  Employs a balanced transformation for stable dynamical 
% system using infinite Grammians for order reduction. If system is not
% stable it asserts.
%   
%   INPUTS: sys = original dynamic system 
%           tol = tolerance considered at the level of the Hankel singular
%           values which determines how many of them are considered in the
%           reduce system based on their relative magnitude towards the
%           highest one of them
%
%   OUTPUTS: sys_reduced = balanced truncated dynamical system
%            res = reduction results returned into a comprehensive
%            MATLAB struct, see parameters commented inline

    % first determine whether the system is stable or not
    if isfield(sys,'poles')
        sA = sys.poles;
    else
        [~, SA] = eig(sys.A);
        sA = diag(SA);  % eig. values of A
    end
    % system is not stable under the tolerance of 10^-12
    assert(sum(abs(real(sA))>=-10^-12)==length(sA)); 

    % now that system is known to be stable, implement balanced tranformation,
    % similar to balreal in MATLAB

    % Step 1: Compute the infinite reachability grammian P using the Lyapunov
    % equation
    res.Plyap = lyap(sys.A,sys.B*sys.B');
    sys.P = res.Plyap;

    % Step 2: Compute the infinite observability grammian Q using the Lyapunov
    % equation
    res.Qlyap = lyap(sys.A',sys.C'*sys.C);
    sys.Q = res.Qlyap;
    
    % Step 3: Compute the orthogonal transformation matrix of
    % Qlyap^0.5*Plyap*Q^0.5 and the eigenvalues of the same matrix, which are
    % the squares of the Hankel sv.
    [V,H2] = eig(res.Qlyap^0.5*res.Plyap*res.Qlyap^0.5);
    % sort the eigenvalues descendingly
    [eV,i] = sort(diag(H2),'descend');
    % interchange the columns of V as needed
    V = V(:,i');
    % get the singular values
    res.hsv = sqrt(eV);
    % normalize them with respect to the highest one
    res.hsv_norm = res.hsv / res.hsv(1);
    % compute the balanced transformation matrix and its inverse
    res.T = H2^-0.25*V'*res.Qlyap^0.5;  % balanced transformation matrix
    res.Ti = inv(res.T);                    % inverse of balanced transformation matrix

    % Step 4: Compute the diagonal, balanced Grammians P = Q
    res.P = res.T*res.Plyap*res.T';
    res.Q = res.Ti'*res.Qlyap*res.Ti;

    % check results under a fixed threshold 10^-8
    assert(sum(sum(res.P-res.Q))<=10^-8);

    % now find minimial number of dimensions which has the property that its
    % last normalized Hankel singular value is above the tolerance threshold 
    % given
    trim_dim = find(res.hsv_norm<=tol);
    if (isempty(trim_dim))
        res.k = length(res.hsv_norm);   % no hsv_norm under the threshold
    elseif (trim_dim(1)-1>=1)
        res.k = trim_dim(1)-1;          % take the first order that is above threshold
    else
        res.k = 1;
    end

    % reduce new system to dimension k
    sys_red = init_dyn_sys(sys.m,res.k,sys.p,'Bal. Red.');
    % first transform A, B, C, D and then truncate to k
    sys_red.A = res.T*sys.A*res.Ti;
    sys_red.B = res.T*sys.B;
    sys_red.C = sys.C*res.Ti;
    sys_red.D = sys.D;
    sys_red.A = sys_red.A(1:res.k,1:res.k);
    sys_red.B = sys_red.B(1:res.k,:);
    sys_red.C = sys_red.C(:,1:res.k);
    sys_red = anal_sys(sys_red);
end