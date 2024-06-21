clc
clear 
close all

% x(k) = sin(3*x(k-1))
% y(k) = x(k)^2
% 注意似然概率是多峰分布，具有强烈的非线性，当y=4时，不知道x=2还是-2

% 生成真实信号和观测数据
t = 0.01:0.01:1;
n = length(t);
x = zeros(1,n);
y = zeros(1,n);
x(1) = 0.1;
y(1) = 0.1^2;
for i = 2:n
   x(i) = sin(3*x(i-1));
   y(i) = x(i)^2 + normrnd(0,0.7);
end
figure;plot(t,x,t,y);
legend('x','y');

% EKF
Xplus = zeros(1,n);
% 设置初值
Pplus = 0.1;
Xplus(1) = 0.1;
Q = 0.0001;
R = 1;
for i = 2:n
    % 预测步
    A = 3 *cos(3*Xplus(i-1));
    Xminus = sin(3 * Xplus(i-1));
    Pminus = A * Pplus * A + Q;
    % 更新步
    C = 2 * Xminus;
    K = Pminus * C / (C * Pminus * C + R);
    Xplus(i) = Xminus + K*(y(i)- Xminus^2);
    Pplus = (eye(1)-K*C)*Pminus;
end
figure;
plot(t,x,'r',t,Xplus,'b');
legend('x','Xplus')
