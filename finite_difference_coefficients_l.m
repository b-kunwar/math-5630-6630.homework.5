function coefficients = finite_difference_coefficients_l(a, h, l)
    % finite_difference_coefficients_l - Compute coefficients for the l-th derivative approximation
    %
    % Inputs:
    %   a  - Vector of distinct coefficients [a0, a1, ..., an]
    %   h  - Step size
    %   l  - The order of the derivative to approximate (1 for first derivative, 2 for second, etc.)
    %
    % Outputs:
    %   coefficients - Vector of coefficients [c0, c1, ..., cn]

    n = length(a) - 1;  % Number of points - 1
    A = zeros(n+1, n+1);
    b = zeros(n+1, 1);

    % Set up the system of equations
    for j = 1:n+1
        for k = 0:n
            A(j, k+1) = a(k+1)^ (j-1);
        end
        % Right-hand side: factorial(l) for row corresponding to j = l+1, otherwise 0
        if j == l + 1
            b(j) = factorial(l);
        end
    end

    % Solve the linear system to find the coefficients
    % coefficients = (A \ b).' / h^l;  % Divide by h^l as per the formula
    coefficients = (A \ b);
    coefficients = coefficients(:);
    
end
