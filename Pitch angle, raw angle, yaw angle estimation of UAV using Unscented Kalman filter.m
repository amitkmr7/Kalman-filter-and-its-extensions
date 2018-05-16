clear;clc
N=80;
theeta=randi([0 180],1,N);
phi=randi([0 360],1,N);
shi=randi([0 90],1,N);
P=60;  %roll rate
Q=68;   %pitch rate
R=90;  %yaw rate

theeta_=(theeta(1)/2);
phi_=(phi(1)/2);
shi_=(shi(1)/2);
q0=(cos(phi_)*cos(theeta_)*cos(shi_))+(sin(phi_)*sin(theeta_)*sin(shi_));
q1=(sin(phi_)*cos(theeta_)*cos(shi_))-(cos(phi_)*sin(theeta_)*sin(shi_));
q2=(cos(phi_)*sin(theeta_)*cos(shi_))+(sin(phi_)*cos(theeta_)*sin(shi_));
q3=(cos(phi_)*cos(theeta_)*sin(shi_))-(sin(phi_)*sin(theeta_)*cos(shi_));
q=[q0; q1; q2; q3];   
q_=[(([0 -P -Q -R;P 0 R -Q;Q -R 0 P;R Q -P 0])/2)*q]; %linearised quarternion state

T = 0.1; % Sampling time
% Step 1: Define UT Scaling parameters and weight vectors
L = 4; % Size of state vector
alpha = 1; % Primary scaling parameter
beta = 2; % Secondary scaling parameter (Gaussian assumption)
kappa = 0; % Tertiary scaling parameter
lambda = alpha^2*(L+kappa) - L;
wm = ones(2*L + 1,1)*1/(2*(L+lambda));
wc = wm;
wm(1) = lambda/(lambda+L);
wc(1) = lambda/(lambda+L) + 1 - alpha^2 + beta;
% Step 2: Define noise assumptions
Q = diag([0 0 4 4]);
R = diag([1 1]);
% Step 3: Initialize state and covariance
x = zeros(4, N); % Initialize size of state estimate for all k
x(:,1) = [0; 0; 50; 50]; % Set initial state estimate
P0 = eye(4,4); % Set initial error covariance
% Simulation Only: Calculate true state trajectory for comparison
% Also calculate measurement vector
w = sqrt(Q)*randn(4, N); % Generate random process noise (from assumed Q)
v = sqrt(R)*randn(2, N); % Generate random measurement noise (from assumed R)
xt = zeros(4, N); % Initialize size of true state for all k
xt(:,1) = [0; 0; 50; 50] + sqrt(P0)*randn(4,1); % Set true initial state
yt = zeros(2, N); % Initialize size of output vector for all k
for k = 2:N
xt(:,k) = [1 0 T 0; 0 1 0 T; 0 0 1 0; 0 0 0 1]*xt(:,k-1) + w(:,k-1);
yt(:,k) = [sqrt((xt(1,k)-q_(1))^2 + (xt(2,k)-q_(2))^2); ...
sqrt((xt(1,k)-q_(3))^2 + (xt(2,k)-q_(4))^2)] + v(:,k);
end
P = P0; % Set first value of P to the initial P0
for k = 2:N
% Step 1: Generating the sigma-points
sP = chol(P,'lower'); % Calculate square root of error covariance
% chi_p = "chi previous" = chi(k-1)
chi_p = [x(:,k-1), x(:,k-1)*ones(1,L)+sqrt(L+lambda)*sP,x(:,k-1)*ones(1,L)-sqrt(L+lambda)*sP];
% Step 2: Prediction Transformation
% Propagate each sigma-point through prediction
% chi_m = "chi minus" = chi(k|k-1)
chi_m = [1 0 T 0; 0 1 0 T; 0 0 1 0; 0 0 0 1]*chi_p;
x_m = chi_m*wm; % Calculate mean of predicted state
% Calculate covariance of predicted state
P_m = Q;
for i = 1:2*L+1
P_m = P_m + wc(i)*(chi_m(:,i) - x_m)*(chi_m(:,i) - x_m)';
end
% Step 3: Observation Transformation
% Propagate each sigma-point through observation
psi_m = [sqrt((chi_m(1,:)-q_(1)).^2 + (chi_m(2,:)-q_(2)).^2); ...
sqrt((chi_m(1,:)-q_(3)).^2 + (chi_m(2,:)-q_(4)).^2)];
y_m = psi_m*wm; % Calculate mean of predicted output
% Calculate covariance of predicted output
% and cross-covariance between state and output
Pyy = R;
Pxy = zeros(L,2);
for i = 1:2*L+1
Pyy = Pyy + wc(i)*(psi_m(:,i) - y_m)*(psi_m(:,i) - y_m)';
Pxy = Pxy + wc(i)*(chi_m(:,i) - x_m)*(psi_m(:,i) - y_m)';
end
% Step 4: Measurement Update
K = Pxy/Pyy; % Calculate Kalman gain
x(:,k) = x_m + K*(yt(:,k) - y_m); % Update state estimate
P = P_m - K*Pyy*K'; % Update covariance estimate
end
% Display results
figure(1);
t = T*(1:N);
subplot(4,1,1)
plot(t,x(1,:),'b--.'); hold on
plot(t,xt(1,:),'r--.') ;hold off
legend('estimate','measurement','location','best');
ylabel('1st quaternion');


subplot(4,1,2);
plot(t,x(2,:),'b--.'); hold on
plot(t,xt(2,:),'r--.'); hold off
legend('estimate','measurement','location','best');
ylabel('2nd quaternion');
subplot(4,1,3)
plot(t,x(3,:),'b--.'); hold on
plot(t,xt(3,:),'r--.'); hold off
legend('estimate','measurement','location','best');
ylabel('3rd quaternion');

subplot(4,1,4)
plot(t,x(4,:),'b--.'); hold on
plot(t,xt(4,:),'r--.'); hold off
legend('estimate','measurement','location','best');
ylabel('4th quaternion');

