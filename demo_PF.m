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
% �����Ӽ���
n = 100;
xold = zeros(1,n);
xnew = zeros(1,n);
xplus = zeros(1,n); % ���ڴ���˲�ֵ������ÿһ�κ�����ʵ�����
w = zeros(1,n);
% ����x0(i),����ֱ������̬�ֲ�������Ҳ��������Գ�ֵ�����У����������Ӷ���ͬ
for i = 1:n
    xold(i) = 0.1;
    w(i) = 1/n;
end
% PF���Ĵ���
for i = 2:100
    % Ԥ�ⲽ����x0�Ƴ�x1
    for j = 1:n
        xold(j) = sin(xold(j)) + 5*xold(j)/(xold(j)^2 + 1) + normrnd(0,1);
    end
    % ���²�
    for j = 1:n
        % w(j) = w(j)*fR(...)
        % fR = (2*pi*R)^(0.5)*exp(-((y(i)-xold(j)^2)^2/(2*R)))
%         w(j) = exp(-((y(i)-xold(j)^2)^2/(2*0.001))); 
        w(j) = exp(-((y(i)-xold(j)^3)^2/(2*0.001))); 
    end
    % ��һ��
    w = w/sum(w);
    % ����w/sum(w)��k*w/sum(k*w)���һģһ������(2*pi*R)^(0.5)Ϊ������
    % ��w(j)���ÿ�ζ����ز�����ÿһ��w(j)���ᱻ����Ϊ1/n��Ҳ�ǳ���
    % ���Կ��Խ�����ȥ��
    
    % �ز���
    % �������˻���һ���̶ȲŻ��ز�����N<1/sum(wj^2);
    % ������ÿ�ζ��ز�������ôӦ�ý�w(j)����ȥ
    % ����c
    c = zeros(1,n);
    c(1) = w(1);
    for j = 2:n
        c(j) = c(j-1) + w(j);
    end
    % ���������������
    % �����ز���n�����ӣ���������֮ǰ��ͬ
    for j = 1:n
        a = unifrnd(0,1); % ���ȷֲ�
        for k = 1:n
            if(a < c(k))
                xnew(j) = xold(k);
                break;
            end
        end
    end
    % �ز������
    % ���µ����Ӹ���xold��Ϊ��һ��������׼��
    xold = xnew;
    % Ȩ�ض���Ϊ1/n
    for j = 1:n
        w(j) = 1/n;
    end
    % ��ÿһ���ĺ������������ֵ��xplus
    xplus(i) = sum(xnew)/n;
end

figure;
plot(t,x,t,xplus,'LineWidth',2);
legend('x','xplus');

% Ч���ܲ��Ϊy�ĸ��ʷֲ���һ�����ֲ��������˲��������������
% ������Ȿ������ʾ���ǿ�ҵķ����ԣ�������ֲ����֣������˲������ܻ�����Ϊ����
% �����˲��ļ����ٶ��Ǵ�Ӳ��
% Ҳ���Գ���û���ز������˲�������ʲô���ģ��Լ�����дд���Լ�����


