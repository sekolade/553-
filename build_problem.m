function prob = build_problem(f_handle, g_ineq_cell, h_eq_cell, opts)
%BUILD_PROBLEM Create a generic constrained problem definition.
%
% Inputs:
%   f_handle    : @(y, params) -> scalar objective f(y)
%   g_ineq_cell : cell array of @(y, params) -> scalar, g_i(y) <= 0
%   h_eq_cell   : cell array of @(y, params) -> scalar, h_j(y) = 0
%   opts (optional):
%     opts.mode         : 'none' | 'penalty' | 'barrier' | 'mixed'
%     opts.Gamma        : SPD matrix (p x p) for inequality penalty
%     opts.W            : SPD matrix (q x q) for equality penalty
%     opts.barrier_type : 'inv' | 'log'

% Output:
%   prob struct that can be attached into "params.prob"

    if nargin < 2 || isempty(g_ineq_cell), g_ineq_cell = {}; end
    if nargin < 3 || isempty(h_eq_cell),   h_eq_cell   = {}; end
    if nargin < 4, opts = struct(); end

    % Defaults
    if ~isfield(opts,'mode') || isempty(opts.mode), opts.mode = 'penalty'; end
    if ~isfield(opts,'barrier_type') || isempty(opts.barrier_type), opts.barrier_type='inv'; end

    p = numel(g_ineq_cell);
    q = numel(h_eq_cell);

    if ~isfield(opts,'Gamma') || isempty(opts.Gamma), opts.Gamma = eye(p); end
    if ~isfield(opts,'W')     || isempty(opts.W),     opts.W     = eye(q); end



    prob.f = f_handle;
    prob.g_ineq = g_ineq_cell;
    prob.h_eq   = h_eq_cell;
    prob.opts   = opts;
end
