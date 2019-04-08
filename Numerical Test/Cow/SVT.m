function D_approx = SVT(D)
[n, ~] = size(D);
Omega = logical(D);
m = sum(sum(Omega));
%tau = 5 * n;
tau = 5 * n;

delta = 1.2 * (n^2/m);
%delta = 10;

MaxIter = 1000;
tol = 1e-4;
l = 5;
k0 = ceil(tau/(delta*normest(D, 1e-2)));

Y = k0*delta*D;
r = 0;

for k = 1:MaxIter
    s = min([r + l, n]);
    
    flag = 0;
    while ~flag
        [U, S, V] = svds(Y, s);
        flag = (S(s,s) <= tau) || ( s == n );
        s = min(s + l, n);
    end
    sigma = diag(S); 
    r = sum(sigma > tau);
    
    X = U(:, 1:r) * diag(sigma(1:r)-tau)*V(:,1:r)';
    %X = X - diag(diag(X));
    if(norm(Omega.*X - D, 'fro')/norm(D, 'fro') <= tol)
        fprintf('convergent');
        break;
    else
        fprintf('%dth iteration: %f\n',k,norm(Omega.*X - D, 'fro')/norm(D, 'fro'));
    end

    Y = Omega.*(Y + delta * (D - X));
    
%     if(mod(k,30) == 0)
%         delta = delta/1.1;
%         %tau = tau * 1.2;
%     end
    
end

D_approx = X;

end



