function A_beta_approx = compute_A_beta(beta, m)
    % Function to approximate A(beta) using Richardson extrapolation
    % for the series: A(beta) = sum_{i=0}^{inf} (-1)^i / (2i + 1)^beta.
    % 
    % Inputs:
    %   beta : The exponent beta in the series (0 < beta <= 1)
    %   m    : Number of terms to use in Richardson extrapolation (should be <= 15)
    %
    % Output:
    %   A_beta_approx : Best approximation of A(beta)
    
    % Step 1: Define the powers vector based on beta
    alpha = beta + (0:(m-2));  % powers vector [beta, beta+1, ..., beta + m - 2]

    % Step 2: Compute the data values a_k = f(2^-k, beta) for k = 1, 2, ..., m
    a = zeros(1, m);
    for k = 1:m
        h = 2^(-k);
        a(k) = f_h_beta(h, beta);  % Calculate f(h, beta)
    end

    % Step 3: Use Richardson extrapolation to get the best approximation of A(beta)
    A_beta_approx = richardson_extrapolation(a, alpha);
end

function result = f_h_beta(h, beta)
    % Compute the truncated series sum f(h, beta) = sum_{i=0}^{1/h} (-1)^i / (2i + 1)^beta
    N = round(1/h);  % Number of terms in the series
    result = sum(((-1).^(0:N)) ./ ((2 * (0:N) + 1).^beta));
end
