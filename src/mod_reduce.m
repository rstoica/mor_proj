function [ sys_red, res ] = mod_reduce( sys, tol )
% mod_reduce Performs modal reduction on a given dynamical system by
% truncating the entire system based on the eigenvalues of the A matrix
% from the state space representation which represent the frequency modes
% of the system
%   
%   INPUTS: sys = system to truncate
%           tol = tolerance threshold for truncation used for the
%           normalized eigenvalues of A with respect to the biggest one
%   OUTPUTS: sys_red = reduced system
%            res = truncationg results such as truncation order 

    % Step1: perform EVD of A matrix
    [V,poles] = eig(sys.A);

    % Step2: sort eigenvalues ascendingly
    [poles, idx] = sort(diag(poles),'ascend');
    % arrange transform matrix according to sorting
    V = V(:,idx);

    % Step3: compute modal coefficients
    res.CT = sys.C*V;
    res.BT = V^-1*sys.B;
    coeff = zeros(sys.n,1);
    for i=1:sys.n
        coeff(i) = norm(res.CT(:,i)*res.BT(i,:))/abs(real(poles(i)));
    end

    % Step4: determine k truncation order
    % sort coefficients descendingly and normalize them wrt to the biggest
    [res.coeff,res.idx] = sort(coeff,'descend');
    % arrange V and poles according to coeff order and recompute CT, BT
    V = V(:,res.idx);
    poles = poles(res.idx);
    res.CT = sys.C*V;
    res.BT = V^-1*sys.B;
    res.coeff_norm = res.coeff / res.coeff(1);
    trim_dim = find(res.coeff_norm<=tol);
    if (isempty(trim_dim))
        res.k = length(res.coeff_norm);   % no hsv_norm under the threshold
    elseif (trim_dim(1)-1>=1)
        res.k = trim_dim(1)-1;          % take the first order that is above threshold
    else
        res.k = 1;
    end
    sys_red = init_dyn_sys(sys.m,res.k,sys.p,'Mod. Red.');
    sys_red.A = diag(poles(1:res.k));
    sys_red.B = res.BT(1:res.k,:);
    sys_red.C = res.CT(:,1:res.k);
    sys_red.D = sys.D;
    sys_red = anal_sys(sys_red);    
end

