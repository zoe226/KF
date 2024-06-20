clc
clear
close all

% PF simulate
% predict x(i) = sin(x(i-1)) + 5 * x(i-1)/(x(i-1)^2+1)+ Q
% measure y(i) = x(i)^2 + R

% generate signal
t = 0.01:0.01:1;
x = zeros(1,100);
y = zeros(1,100);
% initial value
x(1) = 0.1;
y(1) = 0.1^3;
% generate real data and measure data
for i = 2:100
    x(i) = sin(x(i-1)) + 5 * x(i-1)/(x(i-1)^2+1);
%     y(i) = x(i)^2 + normrnd(0,1);
    y(i) = x(i)^3;
end
figure;
plot(t,x,t,y,'LineWidth',2);
legend('x','y');

% PF start
% 设粒子集合
n = 100;
xold = zeros(1,n);
xnew = zeros(1,n);
xplus = zeros(1,n); % 用于存放滤波值，就是每一次后验概率的期望
w = zeros(1,n);
% 设置x0(i),可以直接在正态分布采样，也可以如果对初值有自行，让所有粒子都相同
for i = 1:n
    xold(i) = 0.1;
    w(i) = 1/n;
end
% PF核心代码
for i = 2:100
    % 预测步：由x0推出x1
    for j = 1:n
        xold(j) = sin(xold(j)) + 5*xold(j)/(xold(j)^2 + 1) + normrnd(0,1);
    end
    % 更新步
    for j = 1:n
        % w(j) = w(j)*fR(...)
        % fR = (2*pi*R)^(0.5)*exp(-((y(i)-xold(j)^2)^2/(2*R)))
%         w(j) = exp(-((y(i)-xold(j)^2)^2/(2*0.001))); 
        w(j) = exp(-((y(i)-xold(j)^3)^2/(2*0.001))); 
    end
    % 归一化
    w = w/sum(w);
    % 由于w/sum(w)与k*w/sum(k*w)结果一模一样，且(2*pi*R)^(0.5)为常数，
    % 而w(j)如果每次都做重采样，每一次w(j)都会被设置为1/n，也是常数
    % 所以可以将它们去掉
    
    % 重采样
    % 当粒子退化到一定程度才会重采样，N<1/sum(wj^2);
    % 若不是每次都重采样，那么应该将w(j)加上去
    % 生成c
    c = zeros(1,n);
    c(1) = w(1);
    for j = 2:n
        c(j) = c(j-1) + w(j);
    end
    % 生成随机数，轮盘
    % 首先重采样n个粒子，粒子数与之前相同
    for j = 1:n
        a = unifrnd(0,1); % 均匀分布
        for k = 1:n
            if(a < c(k))
                xnew(j) = xold(k);
                break;
            end
        end
    end
    % 重采样完毕
    % 把新的粒子赋给xold，为下一步递推做准备
    xold = xnew;
    % 权重都设为1/n
    for j = 1:n
        w(j) = 1/n;
    end
    % 把每一步的后验概率期望赋值给xplus
    xplus(i) = sum(xnew)/n;
end

figure;
plot(t,x,t,xplus,'LineWidth',2);
legend('x','xplus');

% 效果很差，因为y的概率分布是一个多峰分布，粒子滤波处理非线性问题
% 如果问题本身的性质就是强烈的非线性，比如多峰分布这种，粒子滤波并不能化腐朽为神奇
% 粒子滤波的计算速度是大硬伤
% 也可以尝试没有重采样的滤波代码是什么样的，自己尝试写写，自己调参


