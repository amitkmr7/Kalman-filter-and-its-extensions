%% M estimation Robust cubature kalman filter using Bhattacharya distance
clear;clc
n=40;
del=1;
R=0.005*eye(4);
Wk=gaussmf(0,R);
q=9*10^-12;
gama=255;
RMSE=[];
xk=randi([20 360],1 ,n);
yk=randi([5 10],1 ,n);
xvk=randi([120 140],1 ,n);
yvk=randi([20 40],1 ,n);
xk_=[xk; yk; xvk; yvk];

xko=randi([15 45],1 ,n+1);
yko=randi([2 5],1 ,n+1);
xkvo=randi([80 100],1 ,n+1);
ykvo=randi([20 40],1 ,n+1);

uk1=[];
uk2=[];
uk3=[];
uk4=[]; 
for i=1:n
    uk_1=(xko(i+1)-xko(i)-del*(xkvo(i)));
    uk1=[uk1;uk_1];
    uk_2=(yko(i+1)-yko(i)-del*(ykvo(i)));
    uk2=[uk2;uk_2];
    uk_3=xkvo(i+1)-xkvo(i);
    uk3=[uk3;uk_3];
    uk_4=ykvo(i+1)-ykvo(i);
    uk4=[uk4;uk_4];
end
ukk=[uk1 uk2 uk3 uk4]';

F=[1 0 del 0; 0 1 0 del; 0 0 1 0; 0 0 0 1];
Q=([((del^3)/3) 0 ((del^2)/2) 0; 0 ((del^3)/3) 0 ((del^2)/2); ((del^2)/2) 0 0 0; 0 ((del^2)/2) 0 0])*q;
Vk=gaussmf(0,Q);

xk1=F*xk_-ukk+Vk;     %dynamic model

%measurement model
zk=[];

for k=1:4
    zk1=atan((xk1(k,:)./yk));
    zk=[zk;zk1];
end

npts=2*n;
zeeta=sqrt(npts/2)*[eye(4,n/2) -eye(4,n/2)];  %made it n/2 to adjust matrix

xkc=xk1;   
Skk = diag([0.9  pi/6 1 1]); % assumed
pkc = Skk*Skk';

xkck=F*xkc+ukk-270;  %estimated value
pkck=F*pkc*F'+Q; %q value can be reduced
[U, S, V] =  svd(pkck);
skk1 = (0.5*[U;V]*sqrt(S));

for i=1:n
xks=sqrt(pkck)*zeeta+xkck;  

zks=[];
for k=1:4
    zks1=atan(xks(k,:)./yk);  
    zks=[zks;zks1];
end
zkck=(1/(2*n))*zks;

pzz=(1/2*n)*(zks-zkck)*(zks-zkck)'+R;


Db=bhattacharyya(zk,zkck);

pkb=inv(((Db*pzz)/gama)-pzz+R);

pxz=(1/2*n)*(xkck-xks)*(zks-zkck)'+R;
pzz=(1/2*n)*(zks-zkck)*(zks-zkck)'+(pkb); %updated pzz
K=(pxz*(inv(pzz)));
s=K';

xkk=xkck+K*(zk-zkck);
pkk=pkck-K*pzz*s;
end

for nn=1:n
RMSE1=(sqrt(mean((xkk(1,nn)-xk(nn)).^2)));
RMSE=[RMSE;RMSE1];
end

figure;
plot(xk,'r--.'); hold on
plot(xkck(1,:),'b--o');hold on
plot(xkk(1,:),'g--o'); hold off
ylabel('direction of angle(degrees)');
xlabel('sampling interval(sec)');
legend('measurement','estimation','predicted','location','best');
suptitle('direction of angle tracking by cubature kalman filter using bhattacharyya distance ');

figure;
plot(RMSE);
ylabel('RMSE ');
xlabel('sampling interval');
title('RMSE value for each sample');





