clc; clear;
N = 10^8;
[sum, error] = find_sum(N)

% % Uncomment to plot the absolute error 
% errors = zeros(1, 8);
% for i = 1:8
%     [~, error] = find_sum(10^i);
%     errors(i) = error;
% end
% 
% figure;
% semilogy(1:8, errors);
% title("Absolute Error = |ln(2) - computed value|");
% xlabel("log_1_0(N)"); ylabel("Error");
% ylim([10^-10 1]);

function [sum, error] = find_sum(n)
    P1 = 1:n/4;
    P2 = n/4+1:n/2;
    P3 = n/2+1:3*n/4;
    P4 = 3*n/4+1:n;

    sign = 1;

    terms1 = zeros(1, length(P1));
    for i = 1:length(P1)
        terms1(i) = sign*(1/P1(i));
        sign = sign * -1;
    end
    sum1 = kahan(terms1);
    
    terms2 = zeros(1, length(P2));
    for i = 1:length(P2)
        terms2(i) = sign*(1/P2(i));
        sign = sign * -1;
    end
    sum2 = kahan(terms2);
    
    terms3 = zeros(1, length(P3));
    for i = 1:length(P3)
        terms3(i) = sign*(1/P3(i));
        sign = sign * -1;
    end
    sum3 = kahan(terms3);
    
    terms4 = zeros(1, length(P4));
    for i = 1:length(P4)
        terms4(i) = sign*(1/P4(i));
        sign = sign * -1;
    end
    sum4 = kahan(terms4);

    sum5 = sum1 + sum3;
    sum6 = sum2 + sum4;
    
    sum = sum5 + sum6;

    actual = log(2);
    error = abs(actual - sum);
end

function sum = kahan(numbers)
    sum = single(0);
    c = single(0);
    for i = 1:length(numbers)
        y = numbers(i) - c;
        t = sum + y;
        c = (t - sum) - y;
        sum = t;
    end
end
