%   Compute weighted non-negative parallel factors with multiplicative updates.
%
%   P = OrthParafac(X,k) computes an estimate of the best rank-k PARAFAC
%   model of a tensor X with nonnegative constraints on the factors.
%   The code requires tensor toolbox version 2.6.
%   The tensor toolbox is available at: %https://www.sandia.gov/~tgkolda/TensorToolbox/index-2.6.html
%   P = OrthParafac(X,k,OPTS) specify options:
%   OPTS.tol: Tolerance on the difference in fit {1.0e-4}
%   OPTS.maxiters: Maximum number of iterations {50}
%   OPTS.dimorder: Order to loop through dimensions {1:ndims(A)}
%   OPTS.init: Initial guess [{'random'}|'nvecs'|cell array]
%   OPTS.printitn: Print fit every n iterations {1}
%   OPTS.wegiht_W: Weight tensor with size of X. {ones(size(X))}
%   OPTS.orthogonal: a vector (length = number of modes of X) with the strength of orthogonality constraint. 0 for no orthogonality constraint {false(1,length(size(X))}.
%   OPTS.fix_mode: a logical vector (length = number of modes of X) with indication if the factor matrix for the mode need to be fixed {false(1,length(size(X)))}.
%
%   [P,U0] = OrthParafac(...) also returns the initial guess.
%
function [P, fit_info, Uinit] = OrthParafac(X, k, opts)
switch nargin
    case 2
        opts = struct;
    case 1
        error('Number of factors (k) is not given.');
end


% Number of modes of X.
N = ndims(X);


fit_info = struct();
fits = [];
ortho = cell(N,1);
for n = 1:N
    ortho{n} = [];
end


% Set options and parameters
fitchangetol = setparam(opts,'tol',1e-4);
maxiters = setparam(opts,'maxiters',500);
dimorder = setparam(opts,'dimorder',1:N);
init = setparam(opts,'init','random');
printitn = setparam(opts,'printitn',1);
orthogonal = setparam(opts, 'orthogonal', false(1,length(size(X))));
epsilon = 1e-12;  % Small number to protect against round-off error
WeightW = setparam(opts, 'weight_W', ones(size(X)));
opt_arrange = setparam(opts, 'arrange', 'true');
fix_mode = setparam(opts, 'fix_mode', false(1,N));


%% Error checking 
% Error checking on maxiters
if maxiters < 0
    error('Maximum iteration (OPTS.maxiter) must be positive');
end


% Error checking on dimorder
if ~isequal(1:N,sort(dimorder))
    error('OPTS.dimorder must include all elements from 1 to ndims(X)');
end


% Error checking weight
if sum(WeightW(:) < 0) > 0
    error('weight (OPTS.weight_W) should be non-negative');
end


% error checking orthogonal constraint
if sum(orthogonal(:) < 0) > 0 || sum(orthogonal(:) > 1) > 0
    error('orthogonal constraint has to between 0 and 1. Revise OPTS.orthogonal')
end


% Error checking fix_mode
if ~islogical(fix_mode) || length(fix_mode) ~= N
    error('OPT.fix_mode has to be a length N logical vector')
end






%% Initialization of factors
if iscell(init)
    for n = dimorder(1:end)
        if ~all(size(Unit{n}) == [size(X,n), k])
            error(['Initial factors (OPTS.init{', int2str(n),'}) has incorrect size']);
        end
    end
    Uinit = init;
elseif isstr(init)
   switch init
      case 'random'
        Uinit = cell(N,1);
        for n = dimorder(1:end)
            Uinit{n} = rand(size(X,n),k) + 0.1;
        end
      case 'nvecs'
        Uinit = cell(N,1);
        for n = dimorder(1:end)
            e = min(k,size(X,n)-2);
            Uinit{n} = abs(nvecs(X,n,e));
            if (e < k)
              Uinit{n} = [Uinit{n} rand(size(X,n),k-e)]; 
            end
        end
     otherwise
        error('Unsupported initialization method specified');
    end
else
    error('Unsupported type of initialization (check OPTS.init)')
end


U = Uinit;
fit = 0;




%% Multiplicative update rules
for iter = 1:maxiters


    fitold = fit;
    % Update across the modes
    for n = dimorder(1:end)
        
        if fix_mode(n) == false


            inds = setdiff(1:N, n);
            Y = khatrirao(U{inds(2)},U{inds(1)});
            Y = ((U{n}*Y').*double(tenmat(tensor(WeightW),n,inds)))*Y;
            Y = (1-orthogonal(n))*Y + ...
                orthogonal(n)*U{n}*mttkrp(X.*WeightW,U,n)'*U{n};
        
            % Initialize factor matrix
            Unew = U{n};


            % Update factors
            tmp = mttkrp(X.*WeightW,U,n);
            Unew = Unew .* (tmp + epsilon);
            Unew = Unew ./ (Y + epsilon);
            U{n} = Unew;
        end
    end


    P = ktensor(U);
    fit = 1-sum(sum(sum((double(X)-double(tensor(P))).^2)))/sum(double(X(:)).^2);
    fits = [fits, fit];
    fitchange = abs(fitold - fit);


    % calculate the orthogonality
    for n=1:N
        ortho{n} = [ortho{n}, norm(U{n}'*U{n} -  eye(size(U{n},2)))];
    end
    
    if mod(iter,printitn)==0
      fprintf(' Iter %2d: fit = %e fitdelta = %7.1e\n', iter, fit, fitchange);
    end


    % Check for convergence
    if (iter > 1) && (fitchange < fitchangetol)
        break;
    end


end




% return the outcome
fit_info.orthogonal = ortho;
fit_info.fits = fits;


if opt_arrange == true
    P = arrange(P);
end


if printitn>0
  fit = 1-sum(sum(sum((double(X)-double(tensor(P))).^2)))/sum(double(X(:)).^2);
  fprintf(' Final fit = %e \n', fit);
end


return;




end




function x = setparam(opts,name,default)
if isfield(opts,name);
    x = opts.(name);
else
    x = default;
end


end