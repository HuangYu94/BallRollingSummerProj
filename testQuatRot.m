function testQuatRot(ang_X,ang_Y)
%This one is written to test the feasibility of using quaternions to make
%rotation simulation
%Micahel Williams, 2015, email:michael.williams.hy@gmail.com
clc
format compact
tic
rads=[10,13,17,19]; %define the rads of balls
numBall=numel(rads);
MaxRad=max(rads);
X=repmat(eye(4),[1,1,numBall]); %Those matrixes stand for balls
psi=zeros(1,numBall);
for n=1:numBall  %make the initial rotation
    q1=angle2quat(0,ang_X/rads(n),ang_Y/rads(n)); 
    q1_inv=quatinv(q1);
    X(:,:,n)=quatmultiply(q1,X(:,:,n));
    X(:,:,n)=quatmultiply(X(:,:,n),q1_inv);
    zaxis=X(4,:,n);
    psi(n)=acos(zaxis(4));
end
err_init=sum(psi);
steps=2000;
errorX=zeros(1,steps);
errorY=zeros(1,steps);
trial_ang=zeros(1,steps);
tempX=repmat(eye(4),[1,1,numBall]);
for i=1:steps  %here we make the data for the plot
    trial_ang(i)=-pi*MaxRad+i*4*pi*MaxRad/steps;
    for n=1:numBall  %rotate the balls around X
        q2X=angle2quat(0,0,trial_ang(i)/rads(n)); %rotate around X axis
        q2X_inv=quatinv(q2X);
        tempX(:,:,n)=quatmultiply(q2X,X(:,:,n));
        tempX(:,:,n)=quatmultiply(tempX(:,:,n),q2X_inv);
        zaxis=tempX(4,:,n);
        psi(n)=acos(zaxis(4)); %calculate the orientation of the Z axis
    end
    errorX(i)=(sum(psi)-err_init)*180/pi;
    for n=1:numBall  %rotate the balls around Y
        q2Y=angle2quat(0,trial_ang(i)/rads(n),0); %rotate around X axis
        q2Y_inv=quatinv(q2Y);
        tempX(:,:,n)=quatmultiply(q2Y,X(:,:,n));
        tempX(:,:,n)=quatmultiply(tempX(:,:,n),q2Y_inv);
        zaxis=tempX(4,:,n);
        psi(n)=acos(zaxis(4)); %calculate the orientation of the Z axis
    end
    errorY(i)=(sum(psi)-err_init)*180/pi;
    trial_ang(i)=trial_ang(i)*180/pi/MaxRad; %convert to degree
end
figure(1);
plot(trial_ang,errorX,trial_ang,errorY);
title('rotation around the X and Y axis to find the optimal distance')
ylabel('error(/deg)');
xlabel('value of d(/deg) which means the rotation angle of the biggest ball');
legend('rotate around X axis','rotate around Y axis',...
    'location','Northeastoutside');
grid on;
toc
end
    
