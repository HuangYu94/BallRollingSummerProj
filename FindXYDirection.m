function FindXYDirection(ang_x,ang_y)
%Make the balls move along X or Y axis, while judge whether the error is
%decreasing.
%This one is mere a simulation.
%Michael Williams 2015, Email: michael.williams.hy@gmail.com
format compact
clc
BallSet=[10,15,20];
numBall=numel(BallSet);
temp=repmat(eye(3),[1,1,numBall]);
z_ori_init=ones(1,numBall);
TurnZ = repmat(eye(3),[1,1,numBall]);
for n = 1:numBall
    theta = pi/12/BallSet(n);
    TurnZ(:,:,n) =RotateZ(theta);
end
X = repmat(eye(3),[1,1,numBall]);
psi = zeros(numBall,1);
%initial rotation of each sphere
for n=1:numBall
    X(:,:,n)=RotateY(ang_y/BallSet(n))*...
        RotateX(ang_x/BallSet(n))*X(:,:,n); %make the initial disturbance
end
for n=1:numBall
    zTurn=pi/6/BallSet(n);
    X(:,:,n)=RotateZ(zTurn)*X(:,:,n);
    zaxis=X(:,:,n)*[0;0;1];
    psi(n) = acos(zaxis(3));
    z_ori_init(n)=abs(psi(n))*360/2/pi; %get the initial orientation
    if(abs(zaxis(3)>1))
        display('pji(n) too big!');
    end
end
d_step=zeros(1,4001);
errorX=zeros(1,4001);
errorY=zeros(1,4001);
for i=1:1:4001
    d_step(i)=-pi*20/2+pi*(i-1)/100;
    for n=1:numBall   %try rotation around the X-axis
        temp(:,:,n)=RotateX(d_step(i)/BallSet(n))*X(:,:,n);
        zaxis=temp(:,:,n)*[0;0;1];
        psi(n)=acos(zaxis(3));
        errorX(i)=sum(psi)*180/pi-sum(z_ori_init); %calculate the error sum of every step
    end
    for n=1:numBall   %try rotation around the Y-axis
        temp(:,:,n)=RotateY(d_step(i)/BallSet(n))*X(:,:,n);
        zaxis=temp(:,:,n)*[0;0;1];
        psi(n)=acos(zaxis(3));
        errorY(i)=sum(psi)*180/pi-sum(z_ori_init);
    end
end
d_degree=d_step*180/pi/20;
plot(d_degree,errorX,'.-',d_degree,errorY);
title('rotate around the X and Y axis to find the optimal d');
ylabel('error(/deg)');
xlabel('value of d(/deg) which means the rotation angle of the biggest ball');
legend('rotate around X axis','rotate around Y axis',...
    'location','Northeastoutside');
grid on;
end




function RxTh = RotateX(theta)
RxTh = [1,  0,  0;
    0, cos(theta), -sin(theta);
    0, sin(theta),  cos(theta)];
end
function RyTh = RotateY(theta)
RyTh = [ cos(theta), 0, sin(theta);
    0,  1,  0;
    -sin(theta), 0, cos(theta)];
end
function RzTh = RotateZ(theta)
RzTh = [cos(theta),  -sin(theta),0;
    sin(theta),   cos(theta),0;
    0,  0,  1];
end
