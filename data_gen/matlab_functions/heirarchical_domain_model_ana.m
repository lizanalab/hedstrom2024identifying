%HEIRARCHICAL DOMAIN MODEL
function [B, twom] = heirarchical_domain_model(A, gamma, a)
    N = size(A,1); % A is square  
    
    if N < 2
        error('N must be at least 2');
    end

    if nargin < 2
        gamma = 1;
    end
    
    if nargin < 3
        a = 0.413;
    end
    
    b = 1/2 - a;  % To keep the matrix normalized

    % Determine the value of 2^n
    two_exp = floor(log2(N));

    ana_diags = zeros(N, 1);
    for d = 0:N-1
        ana_diags(d+1) = hdm_ana(d, two_exp, a);
    end
    matrix = full(spdiags(repmat(ana_diags', N, 1), 0:N-1, N, N));
    matrix = (matrix + triu(matrix, 1)');

    k = full(sum(A));
    twom = sum(k);

    Ds = matrix;

    kD = full(sum(Ds));
    sumD = sum(kD);

    B = @(i) A(:,i) - gamma*Ds(:,i)/(sumD/twom);
end