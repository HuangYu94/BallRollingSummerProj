function mathsimFXYD(ang_x,ang_y)
% this is the math simulation of the function FindXYDirection
% Michael Williams 2015, Email:michael.williams.hy@gmail.com
format compact
clc
BallSet=[10,15,20];
numBall=numel(BallSet);
alfa1=zeros(1,numBall);
for n=1:numBall                %set initial rotation
    alfa1(n)=ang_x/BallSet(n);
end
alfa2=zeros(1,numBall);
for n=1:numBall
    alfa2(n)=ang_y/BallSet(n);
end
d_step=zeros(1,4001);
psi0=zeros(1,numBall);
psi1=zeros(1,numBall);
for n=1:numBall
    psi0(n)=acos(cos(alfa1(n))*cos(alfa2(n)));
end
z_ori_init=sum(abs(psi0));
errorX=zeros(1,4001);
errorY=zeros(1,4001);
for i=1:4001
    d_step(i)=-pi*20/2+pi*(i-1)/100;
    for n=1:numBall   %error record for the X axis rotation
        psi1(n)=acos(-sin(alfa1(n))*sin(d_step(i)/BallSet(n))+...
            cos(alfa1(n))*cos(alfa2(n))*cos(d_step(i)/BallSet(n)));
    end
    errorX(i)=sum(abs(psi1))-z_ori_init;
    errorX(i)=errorX(i)*180/pi;
    for n=1:numBall   %error record for the Y axis rotaion
        psi1(n)=acos(-cos(alfa1(n))*sin(alfa2(n))*sin(d_step(i)/BallSet(n))+...
            cos(alfa1(n))*cos(alfa2(n))*cos(d_step(i)/BallSet(n)));
    end
    errorY(i)=sum(abs(psi1))-z_ori_init;
    errorY(i)=errorY(i)*180/pi;
end
d_degree=d_step*180/pi/20;
plot(d_degree,errorX,d_degree,errorY);
title('rotate around the X and Y axis to find the optimal d(mathsim)');
ylabel('error(/deg)');
xlabel('value of d(/deg) which means the rotation angle of the biggest ball');
legend('rotate around X axis','rotate around Y axis',...
    'location','Northeastoutside');
grid on;
end
