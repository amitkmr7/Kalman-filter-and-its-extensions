clear;clc
%%UAV using extended kalman 
%initial conditions
u=15.264;
to=0;
xu=-0.09426; %stability derevative
xw=5.026;    %Angle of Attack Derivative
g=9.81;      %gravity
zu=0.07848;  %Stability Derivative
zw=-5.31;    %Angle of Attack Derivative
uo=0.8305;
mupmwzu=0.2776;
mwpmwzw=-8.855;
mqpmwuo=-8.201;
del_u=[20 22 23 24 26 27]; %forward velocity
del_w=[60 60 64 67 68 69]; %angle of attack
del_q=[2 3 4 2 4 4];       %pitch rate
del_t=[21 23 21 22 32 32]; %pitch angle
del_de=[3 4 3 5 6 7];      %elevation deflection derivations
k=1;
p=[del_u(k);del_w(k);del_q(k);del_t(k)];
r=[xu xw 0 -g;zu zw uo 0;mupmwzu mwpmwzw mqpmwuo 0;0 0 1 0];
l=r*p;
s=([0 ;-0.2631; -18.59; 0])*(del_de(k));
r=l+s;
ls=[r(1);r(2);r(3);r(4)];   %linearised longitudinal uav system
y=[0 0 0 1]*p;
%% short period approximation can be choosen for the simulation...
%of the aircraft elevator transfer function.
del_a=[20 22 23 24 26 27]; %angle of attack
del_q=[2 3 4 2 4 4];       %pitch rate
del_de=[3 4 3 5 6 7];      %elevation deflection derivations
r=[-5.3293 1;-22.2728 -4.5916];
acc_q=2;
%Process Errors in Process Covaiance Matrix
del_px=20; %initial covariance matrix is choosen intuitively
del_pv=5;
uk=acc_q;
%obsevational error 
del_X=25;
del_VX=6;
Xk=[];
p=[];
for k=1:1:length(del_a)-1;
w=r*[del_a(k);del_q(k)];
s=([-0.5269;-32.9831]*del_de(k))+w;
%The Predicted State
A=[-19.92 -145.94 ;1 0  ;0 1 ];
B=[1;0;0];
C=[0 329.8 1640.4];
Xk_=[del_a(k);del_q(k)];
Xkp1=((A*Xk_));
Xkp2=((B*uk));   
Xkp=(Xkp1+Xkp2);
p=[p;Xkp];
%Initialising Process Covariance Matrix
Pk_=[((del_px).^2) 0;0 ((del_pv).^2)];
Pkp1=((A)*(Pk_));
Pkp=((Pkp1)*(A'));
%Calculating the Kalman gain
R=[((del_X)^2) 0;0 ((del_VX)^2)];
H=[1 0 0; 0 1 0];
K=((Pkp)*H')/((H*Pkp*H')+R);
k=k+1;
Ykm=[del_a(k);del_q(k)];
C=[1 0;0 1];
Yk=C*Ykm;
%Calculating the Current State
Xk=[Xk; Xkp + K*(Yk-(H*(Xkp)))];
%Updating the process covariance matrix
Pk=((eye)-(K*H))*Pkp;
k=k+1;
end
Xkf=[del_a(1)];
Vkf=[del_q(1)];
for k=1:3:3*(length(del_a)-1)
    Xkf=[Xkf;Xk(k)];
end
for k=2:3:3*(length(del_q)-1)
    Vkf=[Vkf;Xk(k)];
end
prx=[del_a(1)];
prv=[del_q(1)];
for i=1:3:(length(p)-1)
    prx=[prx;p(i)];
end
for i=2:2:(length(p))
    prv=[prv;p(i)];
end
t=1:1:length(Xkf);
subplot(2,1,1);
plot(t,Xkf,'r--o');hold on;
plot(t,del_a,'b--o'); hold on;
%plot(t,prx,'g--o');hold off;
xlabel('time(sec)');
ylabel('angle(degrees)');
legend('Estimation','measurement','location','best');
title('angle of attack estimation of UAV at each sec');
subplot (2,1,2);
plot(t,Vkf,'y--o');hold on;
plot(t,del_q,'g--o'); hold on;
%plot(t,prv,'m--o');hold off;
xlabel('time(sec)');
ylabel('angle(degrees)');
title('pitch angle estimation of uav at each second');
legend('Estimation','measurement','location','best');
suptitle('Unmanned Air Vehicle using Extended kalman filter')
