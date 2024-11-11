function c = finite_difference_coefficients(a)
    % Compute finite difference coefficients for approximating the first derivative
    % Inputs:
    %   a : Vector of distinct numbers [a0, a1, ..., an]
    % Output:
    %   c : Vector of coefficients [c0, c1, ..., cn] for f'(x) ≈ (1/h) * (c0*f(x + a0*h) + ... + cn*f(x + an*h))
    
    n = length(a);  % Number of points
    c = zeros(n, 1);  % Initialize coefficients vector
    
    % Loop over each coefficient c_k
    for k = 1:n
        % Initialize the sum for the derivative of Lagrange polynomial L_k'(x)
        L_prime_k = 0;
        
        % Sum over all m != k
        for m = 1:n
            if m ~= k
                term = 1 / (a(k) - a(m));  % Term 1 / (x_k - x_m)
                
                % Product over j != k and j != m
                prod = 1;
                for j = 1:n
                    if j ~= k && j ~= m
                        prod = prod * (a(k) - a(j)) / (a(k) - a(j));
                    end
                end
                
                term = term * prod;
                L_prime_k = L_prime_k + term;
            end
        end
        
        % Store the result as coefficient c_k
        c(k) = L_prime_k;
    end
end
