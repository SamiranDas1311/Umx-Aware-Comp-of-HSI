function D =  CovarianceSimMetric(X,Y,n)
switch n
    case 1
        q = 2;
        D = airm(X,Y); % Affine Invariant Riemennian Metric
    case 2
        D = lerm(X,Y); % Log Euclidean Riemennian Metric
    case 3
        D = jkld(X,Y); % Jeffrey's KL Divergence
    case 4
        D = CholDis(X,Y); % Cholesky Distance
    case 5
        D = jbld(X,Y); % Jensen Bregman Logdet divergence
end   

function D = airm(X,Y)
LogTerm=X.^-0.5*Y*X^0.5;
%D=logm(LogTerm);
D=norm(LogTerm);
end

function D = lerm(X,Y)
LogTerm=log(X)-log(Y);
D=norm(LogTerm);
end

function D = jkld(X,Y)
Term=pinv(X)*Y+pinv(Y)*X-2*eye(size(X,1));
D=trace(Term);
end

function D = CholDis(X,Y)
L1=chol(X);
L2=chol(Y);
D=norm(L1-L2);
end

function D = jbld(X,Y)
D=log(det((X+Y)/2))-0.5*log(det(X*Y));
end

end