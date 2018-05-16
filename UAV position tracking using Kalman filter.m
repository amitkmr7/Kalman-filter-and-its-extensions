%%Target Tracking
clc;
clear all;
%plane is flying but we will be considering its position only in the x
%direction to make it a 2-d case.
%initial state
xo=4000;
vox=280;
%Observations 
X=[4000 4260 4550 4860 5410 5600 5990 6400 6790 7000];
V=[280 282 285 286 290 292 294 296 299 302];
%Process Errors in Process Covaiance Matrix
del_px=20; %initial covariance matrix is choosen intuitively
del_pv=5;
%initial conditions
acc_x=2;
del_t=1;
vx=280;
del_x=25; %uncertainity in the measurement
%Observation Error
del_X=25;
del_VX=6;
Xk=[];
%The Predicted State
A=[1 del_t;0 1];
B=[(0.5*((del_t).^2));del_t];
p=[];
for k=1:1:length(X)-1;
Xk_= [X(k);V(k)];
uk=acc_x;
Xkp1=((A*Xk_));
Xkp2=((B*uk));   
Xkp=(Xkp1+Xkp2);%this is our first estimation
p=[p;Xkp];
%Initialising Process Covariance Matrix
Pk_=[((del_px).^2) 0;0 ((del_pv).^2)];
%Predicted Process Covariance Matrix
Pkp1=((A)*(Pk_));
Pkp2=((Pkp1)*(A'));
pkp=(Pkp2-[0 Pkp2(2);Pkp2(3) 0]); %since the 2nd and 3rd term are not imp.
%Calculating the Kalman gain
R=[((del_X)^2) 0;0 ((del_VX)^2)];
H=[1 0 ; 0 1];
K=((pkp)*H')/((H*pkp*H')+R);
%The New Observation
k=k+1;
Ykm=[X(k);V(k)];
C=[1 0;0 1];
Yk=C*Ykm;
%Calculating the Current State
Xk=[Xk; Xkp + K*(Yk-(H*(Xkp)))];
%Updating the process covariance matrix
Pk1=((eye)-(K*H))*pkp;
pk=(Pk1-[0 Pk1(3);Pk1(2) 0]);
k=k+1;
end
Xkf=[X(1)];
Vkf=[V(1)];
for k=1:2:(length(Xk)-1)
    Xkf=[Xkf;Xk(k)];
end
for k=2:2:(length(Xk))
    Vkf=[Vkf;Xk(k)];
end
prx=[X(1)];
prv=[V(1)];
for i=1:2:(length(p)-1)
    prx=[prx;p(i)];
end
for i=2:2:(length(p))
    prv=[prv;p(i)];
end

Y1=[1200 1300 1480 1590 1700 1800 1990 2090 2200 2300];
V1=[20 22 23 25 27 29 32 34 37 40];
%Process Errors in Process Covaiance Matrix
del_py=20; %initial covariance matrix is choosen intuitively
del_pv1=5;
%initial conditions
acc_y=2;
del_t=1;
vy=280;
del_y=25; %uncertainity in the measurement
%Observation Error
del_Y=25;
del_VY=6;
Yk=[];
%The Predicted State
A=[1 del_t;0 1];
B=[(0.5*((del_t).^2));del_t];
q=[];
for k=1:1:length(Y1)-1;
Yk_= [Y1(k);V1(k)];
uk2=acc_x;
Ykp1=((A*Yk_));
Ykp2=((B*uk2));   
Ykp=(Ykp1+Ykp2);%this is our first estimation
q=[q;Ykp];
%Initialising Process Covariance Matrix
Pk2_=[((del_py).^2) 0;0 ((del_pv1).^2)];
%Predicted Process Covariance Matrix
Pkp12=((A)*(Pk_));
Pkp22=((Pkp1)*(A'));
pkp2=(Pkp22-[0 Pkp22(2);Pkp22(3) 0]); %since the 2nd and 3rd term are not imp.
%Calculating the Kalman gain
R1=[((del_Y)^2) 0;0 ((del_VY)^2)];
H=[1 0 ; 0 1];
K=((pkp2)*H')/((H*pkp2*H')+R1);
%The New Observation
k=k+1;
Ykm1=[Y1(k);V1(k)];
C=[1 0;0 1];
Yk1=C*Ykm1;
%Calculating the Current State
Yk=[Yk; Ykp + K*(Yk1-(H*(Ykp)))];
%Updating the process covariance matrix
Pk12=((eye)-(K*H))*pkp2;
pk2=(Pk12-[0 Pk12(3);Pk12(2) 0]);
k=k+1;
end
Ykf=[Y1(1)];
Vkf1=[V1(1)];
for k=1:2:(length(Yk)-1)
    Ykf=[Ykf;Yk(k)];
end
for k=2:2:(length(Xk))
    Vkf1=[Vkf1;Yk(k)];
end
pry=[Y1(1)];
prv1=[V1(1)];
for i=1:2:(length(q)-1)
    pry=[pry;q(i)];
end
for i=2:2:(length(q))
    prv1=[prv1;q(i)];
end
E=[];
for k=1:1:(length(X))
    E1=abs(((Xkf(k)-X(k))/Xkf(k))*100);
    E=[E;E1];
end

t=1:1:(length(X));
subplot(3,1,1);
plot(X,Y1,'r--o');hold on;
plot(Xkf,Ykf,'g--o'); hold on;
plot(prx,pry,'b--o');hold off;
xlabel('position in x direction (m)');
ylabel('pos in y dir (m)');
legend('measurement','estimation','predicted','location','best');
title('Position estimation of plane in x-y plane');
subplot(3,1,2);
plot(t,V,'m--o'); hold on;
plot(t,Vkf,'b--o'); hold on;
plot(t,prv,'c--o'); hold off;
xlabel('time(per sec)');
ylabel('vel(m/s) in x dir');
legend('measurement','estimation','predicted','location','best');
%title('velocity estimation of plane in x direction');
%subplot(4,1,3);
%plot(t,V1,'m--o'); hold on;
%plot(t,Vkf1,'b--o'); hold on;
%plot(t,prv1,'c--o'); hold off;
%xlabel('time(per sec)');
%ylabel('velocity in y direction');
%legend('measurement','estimation','predicted','location','best');
%title('velocity estimation of plane in y direction');
subplot(3,1,3);
plot(t,E,'m--o'); hold on;
ylim([-2 5]);
xlabel('time(per sec)');
ylabel('Error(in %)');
legend('Error','location','best');
%title('Error in accurate tracking');

suptitle('Tracking of a plane')


