%{
    use KF to solve problem
%}

clc
clear
close all

N = 100;
n = 0:N-1;
s = 7 * n + rand(1,N);
figure;plot(s,'ro');hold on;

y = zeros(1,N);
y(1) = s(1);
Q = 0.0001;
R = 0.25;
y_r = 0;
for i = 2:N
    y_p = 7*y(i-1);
    yp_r = 7^2 * y_r + Q; 
    k = yp_r/(yp_r + R);
    y(i) = y_p + k*( s(i) - y_p);
    y_r = (1-k)*yp_r;
end

plot(y,'b.');hold off;legend('measure','estimate');