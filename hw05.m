% Author: Bikash kunwar / bzk006u7@auburn.ed
% Date: 2024-09-01
% Assignment Name: hw05

classdef hw05
    methods (Static)
        function ret = p1(data, powers)
            % Implementation of Richardson Extrapolation
            % Assume the expansion of f(h) = f(0) + c_1 h^{alpha_1} + c_2 h^{alpha_2} + ... + c_n h^{alpha_n} + ...
            %
            %:param: data: a vector of values f(2^(-i) h), i=1,2,...,n
            %:param: powers: a vector of powers (alpha_1, ..., alpha_{n-1})
            %
            %:return: the extrapolated value of f(0) using Richardson Extrapolation

            % Write your code here.
            % ret = inf;
            ret = richardson_extrapolation(data, powers);
        end

        function R = p2(beta)
            % Compute the value of the series 
            %   sum_{k=0}^(\infty) ((-1)^k /(2k + 1)^{beta})

            %:param: beta: a real number on (0, 1].
            %:return: R: the value of the series
            
            % Write your code here.
            m = 10;
             R = compute_A_beta(beta, m);
        end


        function coefs = p3(shifts)
            % Compute the coefficients of the finite difference scheme for f'(x)
            % using the formula
            % f'(x) \approx \frac{1}{h} (c_0 f(x_0) + c_1 f(x_1) + ... + c_n f(x_n)) + O(h^n)

            %:param: shifts: a vector of shifts (a_0, a_1, ..., a_n), the nodes are x_i = x + a_i h
            %:return: coefs: a vector of coefficients (c_0, c_1, ..., c_n)

            m = length(shifts);
            coefs = zeros(m, 1);

            % write your code here. 

            x = 0; h = 0.1; 
            coefs =  finite_difference_coefficients(x, shifts, h);

        end


        function coefs = p4(shifts, l)
            % For 6630 only.
            % Compute the coefficients of the finite difference scheme for f^{(l)}(x)
            % using the formula
            % f^{(l)}(x) \approx \frac{1}{h^l} (c_0 f(x_0) + c_1 f(x_1) + ... + c_n f(x_n)) + O(h^{n + 1 - l})

            %:param: shifts: a vector of shifts (a_0, a_1, ..., a_n), the nodes are x_i = x + a_i h
            %:param: l: the order of the derivative
            %:return: coefs: a vector of coefficients (c_0, c_1, ..., c_n)

            m = length(shifts);
            coefs = zeros(m, 1);
            h = 0.2;
            % write your code here. 
           coefs = finite_difference_coefficients_l(shifts,h,l);

        end
    end
end