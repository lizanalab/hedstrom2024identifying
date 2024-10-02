%HEIRARCHICAL DOMAIN MODEL
function [B, twom] = heirarchical_domain_model(A, gamma, a, truncate)

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
    
    if nargin < 4
        truncate = true;
    end
    
    b = 1/2 - a;  % To keep the matrix normalized

    % Determine the value of two_exp
    two_exp = ceil(log2(N)) + 1;
    % two_exp = round(log2(N));
    matrix = zeros(2^two_exp, 2^two_exp);

    % Define the basis matrix
    basis = [a, b; b, a];
    matrix(1:2, 1:2) = basis;

    for i = 1:(two_exp - 1)
        matrix(1:2^i, 1:2^i) = matrix(1:2^i, 1:2^i) * a;
        matrix((2^i+1):(2^(i+1)), (2^i+1):(2^(i+1))) = matrix(1:2^i, 1:2^i);

        matrix(1:2^i, (2^i+1):(2^(i+1))) = b / (2^(2*i));
        matrix((2^i+1):(2^(i+1)), 1:2^i) = b / (2^(2*i));
    end

    if truncate
        %matrix = matrix(1:N, 1:N);
        matrix = imresize(matrix, [N N], "bilinear");
        matrix = (matrix + matrix')/2;
    end

    k = full(sum(A));
    twom = sum(k);

    Ds = matrix;

    kD = full(sum(Ds));
    sumD = sum(kD);

    B = @(i) A(:,i) - gamma*Ds(:,i)/(sumD/twom);
end