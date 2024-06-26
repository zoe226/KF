% X(K) = F*X(K-1) + Q
% Y(K) = H*X(K) + R
% generate random signal and use KF filter

clc
clear
close all

%% data generate
t = 0.1:0.01:1;
L = length(t);

% generate real signal and measure x = t^2
x = zeros(1,L);
y = x;
for i = 1:L
    x(i) = t(i)^2;
    y(i) = x(i) + normrnd(0,0.1); % Normal distribution noise
end
figure;plot(t,x,t,y,'LineWidth',2);

 %% filter
 % Observation equation:Y(K) = X(K) + R,R~N(0,1)
 % Prediction equation: the datas are disarray for most situation  
 % for example: X(K) = X(k-1) + Q,Q~N(0,1);
 F1 = 1;
 H1 = 1;
 Q1 = 1;
 R1 = 1; % these initial value can be change
 % initial x(k)+
 Xplus1 = zeros(1,L);
 Xplus(1) = 0.01; % Xplus1(1)~N(0.01,0.01);
 Pplus1 = 0.01^2;
 % KF 
 % X(K)minus = F*X(K-1)plus;
 % P(K)minus = F*P(K-1)plus*F + Q;
 % K = P(K)minus * H' * inv(H * P(K)minus * H' + R); 
 % X(K)plus = X(K)minus + K*(y(K) - H*X(K)minus)
 % P(K)plus = (I - K*H)*P(K)minus;
 for i = 2:L
     % prediction
     Xminus1 = F1 * Xplus1(i-1);
     Pminus1 = F1*Pplus1*F1' + Q1;
     % renewal
     K1 = Pminus1 * H1 / (H1 * Pminus1 * H1' + R1);
     Xplus1(i) = Xminus1 + K1 * (y(i) - H1 * Xminus1);
     Pplus1 = (1 - K1 * H1) * Pminus1;
 end
 figure;plot(t,x,'g',t,y,'b',t,Xplus1,'r');
 
 % another example:X(K) = X(k-1) + X'(k-1)*dt + X''(k-1)*dt^2/(2!) + Q2
 % state is not a single value:X = [X(K) X'(k) X''(k)]T
 % and Y(K) = H2*X + R2, H2 = [1 0 0]
 % prediction equation:
 % X(k) = X(k-1) + X'(k-1)*dt + X''(k-1)*dt^2/(2!) + Q2
 % X'(k) = 0 * X(k-1) +  X'(k-1) + X''(k-1) * dt + Q3;
 % X''(k) = 0 * X(k-1) +  0 * X'(k-1) + X''(k-1) + Q4;
 % F = 1 dt 0.5dt^2
 %     0 1   dt
 %     0 0    1
 % H = [1 0 0]
 % Q = Q2 0 0
 %     0 Q3 0
 %     0 0 Q4 %cov 
 dt = t(2)-t(1);
 F2 = [1,dt,0.5*dt^2;0,1,dt;0,0,1];
 H2 = [1,0,0];
 Q2 = [1,0,0;0,0.01,0;0,0,0.0001];
 R2 = 20;
 Xplus2 = zeros(3,L);
 Xplus2(1,1) = 0.1^2;
 Xplus2(2,1) = 0;
 Xplus2(3,1) = 0;
 Pplus2 = [0.01,0,0;0,0.01,0;0,0,0.0001];
 for i = 2:L
    Xminus2 = F2 * Xplus2(:,i-1);
    Pminus2 = F2 * Pplus2 * F2' + Q2;
    K2 = (Pminus2*H2')/(H2 * Pminus2 * H2' + R2);
    Xplus2(:,i) = Xminus2 + K2*(y(i) - H2 * Xminus2);
    Pplus2 = (eye(3) + K2 * H2) + Pminus2;
 end
  figure;plot(t,x,'g',t,y,'b',t,Xplus1,'r',t,Xplus2(1,:),'k');
 
 
 
 
 
 