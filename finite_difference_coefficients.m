function coefficients = finite_difference_coefficients(x, a, h)
    % finite_difference_coefficients - Compute coefficients for finite difference approximation
    %
    % Inputs:
    %   x  - The point at which we are approximating the derivative
    %   a  - Vector of distinct coefficients [a0, a1, ..., an]
    %   h  - Step size
    %
    % Outputs:
    %   coefficients - Vector of coefficients [c0, c1, ..., cn]

    n = length(a) - 1;            % Number of coefficients - 1
    coefficients = zeros(1, n+1); % Initialize vector to store c_k

    % Loop over each k to compute c_k = h * L'_k(x)
    for k = 1:n+1
        % Compute L'_k(x) for the current k
        x_k = x + a(k) * h;       % x_k = x + a_k * h
        Lk_prime = 0;             % Initialize L'_k(x)

        % Loop over m != k to calculate L'_k(x)
        for m = 1:n+1
            if m == k
                continue;  % Skip m = k
            end

            x_m = x + a(m) * h;   % Compute x_m
            product_term = 1;

            % Compute the product term for j != k and j != m
            for j = 1:n+1
                if j ~= k && j ~= m
                    x_j = x + a(j) * h;
                    product_term = product_term * (x - x_j) / (x_k - x_j);
                end
            end

            % Add to L'_k(x) for this m
            Lk_prime = Lk_prime + (1 / (x_k - x_m)) * product_term;
        end

        % Calculate c_k = h * L'_k(x)
        coefficients(k) = h * Lk_prime;
        coefficients = coefficients(:);
    end
end
