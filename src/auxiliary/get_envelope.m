function W = get_envelope(S, M, u)

r = size(M, 1);
U = S * S';

[w, D] = eig(U);
[DM ind] = max(diag(D));
W = w(:, ind);

if u > 1
    
    for i = 2 : u
        
        Eu = grams(M * W);
        QEu = eye(r) - Eu * inv(Eu' * Eu) * Eu';
        [w, D] = eigs(QEu * U * QEu, 1);
        W = [W w];
        
    end
    
end


