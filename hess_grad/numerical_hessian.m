function H = numerical_hessian(y, mu, params)
%NUMERICAL_HESSIAN Hessian via gradient differences (central):
%   H(:,j) = (g(y + h e_j) - g(y - h e_j)) / (2h)
%
% Optional:
%   params.eps_hess : base FD step (default 1e-5)

    if isfield(params,'eps_hess') && ~isempty(params.eps_hess)
        eps_step = params.eps_hess;
    else
        eps_step = 1e-5;
    end

    n = length(y);
    H = zeros(n,n);

    for j=1:n
        h = eps_step * max(1, abs(y(j)));

        y_plus  = y; y_plus(j)  = y_plus(j)  + h;
        y_minus = y; y_minus(j) = y_minus(j) - h;

        g_plus  = numerical_gradient(y_plus,  mu, params);
        g_minus = numerical_gradient(y_minus, mu, params);

        H(:,j) = (g_plus - g_minus) / (2*h);
    end

    % Symmetrize to reduce numerical asymmetry
    H = 0.5 * (H + H');
end
