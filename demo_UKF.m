clc
clear
close all

% UKF�㷨
% �Գ�ѧ�߲�̫�Ѻã������ࡢ�����
% plus��+ ,minus��-, hat������
t = 0.01:0.01:1;
% ������ʵ��x��۲�z
x = zeros(2,100);
z = zeros(2,100);
x(1,1) = 0.2;
x(2,1) = 0.2;
% ���дһ�������Ե�״̬���̺͹۲ⷽ�̣�xk= f(xk-1) zk = h(xk) + R
for i = 2:100
    x(1,i) = sin(x(1,i-1)) + 2 * cos(x(2,i-1));
    x(2,i) = 3 * cos(x(1,i-1)) + sin(x(2,i-1));
%     z(1,i) = x(1,i) + x(2,i)^3 + normrnd(0,1);
%     z(2,i) = x(1,i)^3 + x(2,i) + normrnd(0,1);    
    z(1,i) = x(1,i) - x(2,i) + normrnd(0,1);
    z(2,i) = x(1,i) + x(2,i) + normrnd(0,1);   
end
% figure;plot(t,z(1,:),t,z(2,:));
% figure;plot(t,x(1,:),t,x(2,:));

% ���ֵ
X = zeros(2,100);
X(1,1) = 0.1;
X(2,1) = 0.2;
Pplus = eye(2);
Q = 0.1*eye(2);
R = 0.01*eye(2);
% ��w(i)
n = 2; % ����X��ά��
w = zeros(n,2*n + 1);
lamda = 1;
for i = 1:2*n + 1
    w(i) = 1/(2*(n + lamda));
end
w(1) = lamda/(n+lamda);

% UKF start
for i = 2:100
    % ��sigma
    xsigma = zeros(n,2*n+1);
    L = chol(Pplus);
    xsigma(:,1) = X(:,i-1);
    for j = 1:n
        xsigma(:,j+1) = xsigma(:,1) + sqrt(n+lamda) * L(:,j);
        xsigma(:,j+1+n) = xsigma(:,1) - sqrt(n+lamda) * L(:,j);
    end
    % Ԥ�ⲽ ����xsigmaminus
    xsigmaminus = zeros(n,2*n+1);
    for j = 1:2*n+1
        % xk = f(xk-1)
        xsigmaminus(1,j) = sin(xsigma(1,j)) + 2*cos(xsigma(2,j));
        xsigmaminus(2,j) = 3 * cos(xsigma(1,j)) + sin(xsigma(2,j));
    end
    % ����������
    xhatminus = zeros(n,1);
    Pminus = zeros(n,n);
    for j = 1:2*n+1
        xhatminus = xhatminus + w(j)*xsigmaminus(:,j);
    end
    for j = 1:2*n+1
        Pminus = Pminus + w(j) * (xsigmaminus(:,j) - xhatminus) * (xsigmaminus(:,j)-xhatminus)';
    end
    Pminus = Pminus + Q;
    % Ԥ�ⲽ����
    % ���²���ʼ
    % �ٲ�sigma��
    xsigma = zeros(n,2*n + 1);
    xsigma(:,1) = xhatminus;
    xsigma(:,1) = xhatminus;
    L1 = chol(Pminus);
    for j = 1:n
        xsigma(:,j+1) = xsigma(:,1) + sqrt(n + lamda)*L1(:,j);
        xsigma(:,j+1+n) = xsigma(:,1) - sqrt(n + lamda)*L1(:,j);
    end
    % ����y,yhat
    yhat = zeros(n,1);
    y = zeros(n,1);
    for j = 1:2*n+1
%         y(1,j) = xsigma(1,j) + xsigma(2,j)^3;
%         y(2,j) = xsigma(1,j)^3 + xsigma(2,j);
        y(1,j) = xsigma(1,j) - xsigma(2,j);
        y(2,j) = xsigma(1,j) + xsigma(2,j);
        yhat = yhat + w(j)*y(:,j);
    end
    % ��Py Pxy
    Py = zeros(n,n);
    Pxy = zeros(n,n);
    for j = 1:2*n + 1
        Pxy = Pxy + w(j) * (xsigma(:,j) - xhatminus) * (y(:,j) - yhat)';
        Py = Py + w(j) * (y(:,j) - yhat) * (y(:,j) - yhat)';
    end
    Py = Py + R;
    % �󿨶�������
    K = Pxy / Py;
    % �۲�����
    Y = zeros(n,1);
    Y(1,1) = z(1,i);
    Y(2,1) = z(2,i);
    % ���²�
    X(:,i) = xhatminus + K*(Y - yhat);
    Pplus = Pminus + K * Py * K';    
end
figure;
plot(t,x(1,:),t,X(1,:));
legend('x','X');
figure;
plot(t,x(2,:),t,X(2,:));
legend('x','X');