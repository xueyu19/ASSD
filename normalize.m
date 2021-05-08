function [X,DE,y] = normalize(X,y)
[~,p] = size(X);
DE = zeros(p,1);
if nargin == 1
    y = [];
    for k =1:p
        Xk = X(:,k);
        dk = 1/norm(Xk);
        X(:,k) = Xk*dk;
        DE(k) = dk;
    end
else
    y = y - mean(y);
    for k =1:p
        Xk = X(:,k);
        Xk = Xk - mean(Xk);
        dk = 1/norm(Xk);
        X(:,k) = Xk*dk;
        DE(k) = dk;
    end
end
DE = diag(DE);