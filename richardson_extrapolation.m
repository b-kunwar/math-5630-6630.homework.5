function f0_approx = richardson_extrapolation(a, alpha)

m = length(a);

A = zeros(m,m);
b = zeros(m,1);
b(1)= 1;
A(1,:) = 1;

for i = 1: m-1
    for j = 1:m
        A(i+1,j) = (2^(-j))^alpha(i);
        % form A matrix
    end
end

lambda = A\b;
f0_approx = sum(lambda.*a(:));
end
