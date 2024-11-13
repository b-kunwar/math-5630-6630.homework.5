function Lk_prime = Lk_derivative(x, k, a, h)
    % x  - point at which we are evaluating the derivative
    % k  - index of the Lagrange basis polynomial
    % a  - vector of distinct coefficients [a0, a1, ..., an]
    % h  - step size

    n = length(a) - 1;  % Number of nodes - 1
    x_k = x + a(k) * h; % Compute x_k

    % Initialize L'_k(x) as zero
    Lk_prime = 0;

    % Loop over m != k
    for m = 1:n+1
        if m == k
            continue;  % Skip m = k
        end

        x_m = x + a(m) * h;  % Compute x_m
        product_term = 1;

        % Compute the product over j != k, j != m
        for j = 1:n+1
            if j ~= k && j ~= m
                x_j = x + a(j) * h;
                product_term = product_term * (x - x_j) / (x_k - x_j);
            end
        end

        % Add the term for this m to L'_k(x)
        Lk_prime = Lk_prime + (1 / (x_k - x_m)) * product_term;
    end
end
