function QuatMakeErrDecrease(ang_X,ang_Y,nr)
% based on the Balls2Sim, we use this one to insure the error function to
% decrease every step
format compact
clc
tic
numStep =100000;
rads = [10,13,29,17];
numBall=numel(rads);
stp=zeros(numStep,1);
path=zeros(2,numStep);
CircleRoll=zeros(2,nr);
error_rec=zeros(1,numStep);
X = repmat(eye(3),[1,1,numBall]);
psi = zeros(numBall,1);
k=0;
%this one signify the direction to take for every step
%1st one for direction, second one for step length
%initial rotation of each sphere, take the parameter
Zturn=600*pi/nr; %nr denote how many steps you need to rotate 
% to turn the bigest ball around
for n=1:numBall
    qrotX=angle2quat(0,0,ang_X/rads(n));%generate the quaternions for rotation
    qrotY=angle2quat(0,ang_Y/rads(n),0);
    X(:,:,n)=quatrotate(qrotX,X(:,:,n));
    X(:,:,n)=quatrotate(qrotY,X(:,:,n));
    zaxis=X(:,:,n)*[0;0;1];
    psi(n)=acos(zaxis(3));
end
nn=10*nr;
for i=1:nn
    err_delta=sum(psi);
    error_new=err_delta;
    while err_delta>0.0000001
        error_old=error_new;
        k=k+1;
        %we here to rotate the ball around either X or Y axis whichever to
        %decrease the overall error
        DirAndLen=FindXYLen(X,rads);
        if DirAndLen(1)==0
            for n=1:numBall  %rotate around X axis
                qrotX=angle2quat(0,0,DirAndLen(2)/rads(n));
                X(:,:,n)=quatrotate(qrotX,X(:,:,n));
            end
            if k==1
                path(2,k)=DirAndLen(2);
                path(1,k)=0;
            else
                path(2,k)=DirAndLen(2)+path(2,k-1);
                path(1,k)=path(1,k-1);
            end
        else   %rotate around Y axis
            for n=1:numBall
                qrotY=angle2quat(0,DirAndLen(2)/rads(n),0);
                X(:,:,n)=quatrotate(qrotY,X(:,:,n));
            end
            if k==1
                path(1,k)=DirAndLen(2);
                path(2,k)=0;
            else
                path(1,k)=DirAndLen(2)+path(1,k-1);
                path(2,k)=path(2,k-1);
            end
        end
        stp(k)=k;
        for n=1:numBall
            zaxis=X(:,:,n)*[0;0;1];
            psi(n) = acos(zaxis(3));
        end
        error_rec(k)=sum(abs(psi));
        error_new=error_rec(k);
        err_delta=error_old-error_new;
    end
    for n=1:numBall %if the error cannot decrease any more, we rotate all around Z axis
        qrotZ=angle2quat(Zturn/rads(n),0,0);
        X(:,:,n)=quatrotate(qrotZ,X(:,:,n));
        zaxis=X(:,:,n)*[0;0;1];
        psi(n)=acos(zaxis(3));
    end
    CircleRoll(:,i)=path(:,k);
    
    
    
path1=path(:,1:k-1);
save('mydate.mat','path','path1','error_rec','CircleRoll');
stp1=stp(1:k-1);
error_rec1=error_rec(1:k-1)*180/pi;
figure(1);
plot(stp1,error_rec1);
title('Noiseless Ensemble Control of 4 Spheres Orientation');
xlabel('steps');
ylabel('overall error(degs)');
figure(2)
plot(path1(1,:),path1(2,:));
hold on;
plot(CircleRoll(1,:),CircleRoll(2,:),'ro');
title('trajectory of the panel for the control process');
xlabel('motion projected on the X axis')
ylabel('motion projected on the Y axis')
hold off 
end
toc
end


function DirAndLen=FindXYLen(X,rads)
%this function blindly try every step length and compare the minimums to
%determine whether to rotate around X axis or Y axis.
numBall=numel(rads);
Mrad=max(rads);
temp=repmat(eye(3),[1,1,numBall]);
d_step=zeros(1,2001);
errorX=zeros(1,2001);
errorY=zeros(1,2001);
psi = zeros(numBall,1);
for i=1:1:2001
    d_step(i)=-pi*Mrad+2*pi*(i-1)*Mrad/2000;
    for n=1:numBall   %try rotation around the X-axis
        qrotX=angle2quat(0,0,d_step(i)/rads(n));
        temp(:,:,n)=quatrotate(qrotX,X(:,:,n));
        zaxis=temp(:,:,n)*[0;0;1];
        psi(n)=acos(zaxis(3));
        errorX(i)=sum(abs(psi)); %calculate the error sum of every step
    end
    for n=1:numBall   %try rotation around the Y-axis
        qrotY=angle2quat(0,d_step(i)/rads(n),0);
        temp(:,:,n)=quatrotate(qrotY,X(:,:,n));
        zaxis=temp(:,:,n)*[0;0;1];
        psi(n)=acos(zaxis(3));
        errorY(i)=sum(abs(psi));
    end
end
[a,b]=min(errorX);
[c,d]=min(errorY);
if a<=c
    DirAndLen=[0,d_step(b)];
else
    DirAndLen=[1,d_step(d)];
end
end

% function RxTh = RotateX(theta)
%      RxTh = [1,  0,  0;
%              0, cos(theta), -sin(theta);
%              0, sin(theta),  cos(theta)];
%   end
%  function RyTh = RotateY(theta)
%      RyTh = [ cos(theta), 0, sin(theta);
%               0,  1,  0;
%              -sin(theta), 0, cos(theta)];
%  end
% function RzTh = RotateZ(theta)
% RzTh = [cos(theta),  -sin(theta),0;
%     sin(theta),   cos(theta),0;
%     0,  0,  1];
% end
    