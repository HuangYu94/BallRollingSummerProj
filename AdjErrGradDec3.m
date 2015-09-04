function AdjErrGradDec3(ang_X,ang_Y)
%Michael Williams 2015, Email: michael.williams.hy@gmail.com
%
%The gradient descent method uses a flexible step length
% for X-axis rotation and the Y axis rotation.
% Available at the github at:  git@github.com:HuangYu94/BallRollingSummerProj.git

% profile the code-- what is slow?
% error should monotonically decrease why doesn't it?  -- we will switch to
% plot the RMSE
% why do the balls never roll in the positive x or y direction? 
%    you don't have a unimodal function, so golden section may give the
%    wrong results
% 
%
%  How to compare methods: 
%      which has the fewest number of z-turns?  (or the shortest 
%         sum(abs(z-angle_rotations))
%      which has the shortest overall path length?  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format compact
syms theta  % this is needed for optimization
theta = 5;
if nargin<2
    ang_X = 7*pi;ang_Y = 9*pi;
end

clc
tic
numStep=100000;  %prelocation of numStep to increase the speed
rads = [10,13,17,19];
numBall=numel(rads);
stp=zeros(1,numStep);
path=zeros(2,numStep);
error_rec=zeros(1,numStep);

 RMSE=zeros(1,numStep);
 RMSEinDEG=zeros(1,numStep);

X = repmat(eye(3),[1,1,numBall]);
psi = zeros(numBall,1);
for n=1:numBall %set the initial angle of balls
    X(:,:,n)=RotateX(ang_X/rads(n))*RotateY(ang_Y/rads(n));
    zaxis=X(:,:,n)*[0;0;1];
    psi(n) = acos(zaxis(3));
end
k=0; % number of steps
% Zturn=pi*600/nr; %nr denote how many steps you need to rotate
% to turn the bigest ball around
precision_control=0.01*numBall;  %stop the process when this error is reached
err_new=10;%make the loop start
while err_new>precision_control
    for n=1:numBall  %calculate the initial error
        Z_orient=X(:,:,n)*[0;0;1];
        psi=acos(Z_orient(3));
    end
    err_delta=sum(psi.^2);
    while err_delta>0.001
        k=k+1;
        for n=1:numBall %calculate the current error
            zaxis=X(:,:,n)*[0;0;1];
            psi(n) = acos(zaxis(3));
        end
        err_old=sum(psi.^2);
        Xturn=GDFindLengthX(X,rads,numBall);%find best length for X rotation
        if k==1
            path(2,k)=Xturn;  %record the trajectories
            path(1,k)=0;
        else
            path(2,k)=Xturn+path(2,k-1);
            path(1,k)=path(1,k-1);
        end
        for n=1:numBall  %apply the rotation.   %GOAL: use quaternions
            X(:,:,n)=RotateX(Xturn/rads(n))*X(:,:,n);
            zaxis=X(:,:,n)*[0;0;1];
            psi(n)=acos(zaxis(3));
        end
        err_new=sum(psi.^2);
        error_rec(k)=err_new;
        RMSE(k) = sqrt(1/numBall*sum(psi.^2));
        RMSEinDEG(k) = RMSE(k)*180/pi;
        stp(k)=k;     
        k=k+1;

        Yturn=GDFindLengthY(X,rads,numBall);%find a proper length for Y rotation
        for n=1:numBall
            X(:,:,n)=RotateY(Yturn/rads(n))*X(:,:,n);
            zaxis=X(:,:,n)*[0;0;1];
            psi(n)=acos(zaxis(3));
        end
        if k==1
            path(1,k)=Yturn;
            path(2,k)=0;
        else
            path(1,k)=Yturn+path(1,k-1);
            path(2,k)=path(2,k-1);
        end
        err_new=sum(psi.^2);
        err_delta=abs(err_old-err_new);
        error_rec(k)=err_new;
        
        RMSE(k) = sqrt(1/numBall*sum(psi.^2));
        RMSEinDEG(k) = RMSE(k)*180/pi;
        
        stp(k)=k;
        %the following just for drawing
        %error_recDEGREE=error_rec(1:k-1)*180/pi;
        figure(1);
        plot(stp(1:k-1),RMSEinDEG(1:k-1));
        title('Noiseless Ensemble Control of 4 Spheres Orientation');
        xlabel('steps (x,y, or z rotations)');
        ylabel('sum of psi^2 error (degs)');  %  sqrt(1/numBall*sum( psi.^2))  %TODO: switch to this error metric (RMSE)
        path1=path(:,1:k-1);
        save('GDmyData.mat','error_rec','path','path1');
        figure(2)
        plot(path1(1,:),path1(2,:), path1(1,1),path1(2,1),'go',path1(1,end),path1(2,end),'rx' );
        title('movement of the panel for the control process');
        legend ('path','start','end')
        xlabel('motion projected on the X axis')
        ylabel('motion projected on the Y axis')
    end
    Zturn=Z_align(X,rads);  %find out the optimal angle to align the balls`
    % angle theta together.
    for n=1:numBall
        X(:,:,n)=RotateZ(Zturn/rads(n))*X(:,:,n);
        Z_orient=X(:,:,n)*[0;0;1];
        psi=acos(Z_orient(3));
    end
end
%we need to get the error_rec rid of zero, in order to make a plot
% we do the following steps:
stp1=stp(1:k-1);
%error_recDEGREE=error_rec(1:k-1)*180/pi; %convert from radians to degrees
figure(1);
plot(stp1,RMSEinDEG);
title('Noiseless Ensemble Control of 4 Spheres Orientation');
xlabel('steps');
ylabel('overall error(degs)');
path1=path(:,1:k-1);
save('GDmyData.mat','error_rec','path','path1');
figure(2)
plot(path1(1,:),path1(2,:), path1(1,1),path1(2,1),'go',path1(1,end),path1(2,end),'rx' );
title('movement of the panel for the control process');
legend ('path','start','end')
xlabel('motion projected on the X axis')
ylabel('motion projected on the Y axis')
toc

function Xturn=GDFindLengthX(X,rads,numBall)
%find a proper rotation angle for X axis use gradient descent method with a
%flexible step length
% INPUTS:
% X: Rotation matrices for each ball
% ... 

oldx=-1;
newx=0;
precision=1e-6;
ff=cell(1,numBall);

for nx=1:numBall
    temp=RotateX(theta/rads(nx))*X(:,:,nx); %effect of symbolic rotation -- theta is a free parameter
    ff{nx}=acos(temp*[0;0;1]).^2;  %psi as a function of theta.
end
while abs(newx-oldx)>precision
    oldx=newx;
    f_deriv=0;
    for nx=1:numBall
        temp=vpa(subs(ff{nx},theta,newx));
        f_deriv=f_deriv+temp(3);
    end
    gamma=FindGammaX(f_deriv,X,numBall,newx,rads);
    newx=newx-f_deriv*gamma;  %follow negative gradient
end
Xturn=newx;
end

function gamma=FindGammaX(f_deriv,X,numball,newVar,rads)
%I want to implement the 0.618 method, because I cannot add the function
%handle in matlab
%  MICHAEL -- the function is not strictly unimodal.
a=0;
b=2;
GR=(sqrt(5)-1)/2; %global variable for Golden Section Search (https://en.wikipedia.org/wiki/Golden_section_search)
d=GR*(b-a)+a;
c=b-GR*(b-a);
while abs(c-d)>0.0000000001
    ffc=0;
    ffd=0;
    for ng=1:numball
        temp1=RotateX((newVar-f_deriv*c)/rads(ng))*X(:,:,ng);
        temp=temp1*[0;0;1];
        obj_temp=acos(temp(3)).^2;
        ffc=ffc+obj_temp; %get the objective function value for later comparison
        temp1=RotateX((newVar-f_deriv*d)/rads(ng))*X(:,:,ng);
        temp=temp1*[0;0;1];
        obj_temp=acos(temp(3)).^2;
        ffd=ffd+obj_temp;
    end
    if ffc<=ffd
        b=d;
        d=c;
        c=b-GR*(b-a);
    else
        a=c;
        c=d;
        d=a+GR*(b-a);
    end
end
gamma=(a+b)/2;
end

function Yturn=GDFindLengthY(X,rads,numBall)
%find a proper rotation angle for Y axis
oldy=-1;
newy=0;
precision=0.000001; %TODO: make this a global variable -- this is too high
ff=cell(1,numBall);
for ny=1:numBall
    temp=RotateY(theta/rads(ny))*X(:,:,ny);
    ff{ny}=acos(temp*[0;0;1]).^2;
end
while abs(newy-oldy)>precision
    oldy=newy;
    f_deriv=0;
    for ny=1:numBall
        temp=vpa(subs(ff{ny},theta,newy));
        f_deriv=f_deriv+temp(3);
    end
    gamma=FindGammaY(f_deriv,X,numBall,newy,rads);
    newy=newy-f_deriv*gamma;
end
Yturn=newy;
end

function gamma=FindGammaY(f_deriv,X,numball,newVar,rads)
%I want to implement the 0.618 method, because I cannot add the function
%handle in matlab
GR=(sqrt(5)-1)/2; %global varible for Golden Section Search
a=0;
b=2;
d=GR*(b-a)+a;
c=b-GR*(b-a);
while abs(c-d)>0.0000000000001
    ffc=0;
    ffd=0;
    for ng=1:numball
        temp1=RotateY((newVar-f_deriv*c)/rads(ng))*X(:,:,ng);
        temp=temp1*[0;0;1];
        obj_temp=acos(temp(3)).^2;
        ffc=ffc+obj_temp; %get the objective function value for later comparison
        temp1=RotateY((newVar-f_deriv*d)/rads(ng))*X(:,:,ng);
        temp=temp1*[0;0;1];
        obj_temp=acos(temp(3)).^2;
        ffd=ffd+obj_temp;
    end
    if ffc<=ffd
        b=d;
        d=c;
        c=b-GR*(b-a);
    else
        a=c;
        c=d;
        d=a+GR*(b-a);
    end
end
gamma=(a+b)/2;
end
%==================NOTICE====================================
% I chose to not use the function in the follwing.
% function Zturn=Z_align(X,rads)
% % this function was written for finding the optimal angle for rotation
% % around Z axsis to make the balls align together to the positive direction
% % of X axsis.
% %Michael Williams 2015, e-mail:michael.williams.hy@gamil.com
% numBall=numel(rads);
% sin_sum=0;
% cos_sum=0;
% obj_f=cell(1,numBall);
% for n=1:numBall  %calculate the average angle
%     temp=X(:,:,n)*[0;0;1];
%     sin_sum=sin_sum+temp(2);
%     cos_sum=cos_sum+temp(1);
% end
% sin_av=sin_sum./numBall;
% cos_av=cos_sum./numBall;
% syms theta1
% for n=1:numBall   %calculate the objctive funtion for every ball and store them in obj_f
%     temp=RotateZ(theta1/rads(n))*X(:,:,n);
%     temp1=temp*[0;0;1];
%     obj_f{n}=atan2(real(temp1(2))-sin_av,real(temp(1))-cos_av).^2;
% end
% precision=0.00001; %implement the gradient descent method
% newz=0;
% oldz=-1;
% while abs(newz-oldz)>precision
%     oldz=newz;
%     f_deriv=0;
%     for n=1:numBall
%         temp=vpa(subs(obj_f{n},theta1,newz));
%         f_deriv=f_deriv+temp;
%     end
%     gamma=FindGammaZ(f_deriv,X,numBall,newz,rads);
%     newz=newz-f_deriv*gamma;
% end
% Zturn=newz;
% end
%
% function gamma=FindGammaZ(f_deriv,X,numBall,newz,rads)
% %still I use the golden section search to find the proper gamma
% sin_sum=0;
% cos_sum=0;
% for n=1:numBall  %calculate the average angle
%     temp=X(:,:,n)*[0;0;1];
%     sin_sum=sin_sum+temp(2);
%     cos_sum=cos_sum+temp(1);
% end
% sin_av=sin_sum./numBall;
% cos_av=cos_sum./numBall;
% GR=(sqrt(5)-1)/2; %global varible for Golden Section Search
% a=0;
% b=2;
% d=GR*(b-a)+a;
% c=b-GR*(b-a);
% while abs(c-d)>0.0000000001
%     ffc=0;
%     ffd=0;
%     for n=1:numBall
%         temp1=RotateZ((newz-f_deriv*c)/rads(n))*X(:,:,n);
%         temp=temp1*[0;0;1];
%         obj_temp=atan2(real(temp(2))-sin_av,real(temp(1))-cos_av).^2;
%         ffc=ffc+obj_temp; %get the objective function value for later comparison
%         temp1=RotateZ((newz-f_deriv*d)/rads(n))*X(:,:,n);
%         temp=temp1*[0;0;1];
%         obj_temp=atan2(real(temp(2))-sin_av,real(temp(1))-cos_av).^2;
%         ffd=ffd+obj_temp;
%     end
%     if ffc<=ffd
%         b=d;
%         d=c;
%         c=b-GR*(b-a);
%     else
%         a=c;
%         c=d;
%         d=a+GR*(b-a);
%     end
% end
% gamma=(a+b)/2;
% end
% the lesson is that gradient descent method is not suitatble for this task.
% ============================NOTICE==============================
%so I will switch to the following method:
function Zturn=Z_align(X,rads)
%for this function I still use the blindly trial method to find the
%alignment angle.
%Michael Williams 2015, e-mail:michael.williams.hy@gamil.com
StepLength=0.01*pi; %length for trial
trial=0;  %initialize
numBall=numel(rads);
Xtemp=repmat(eye(3),[1,1,numBall]);
theta2=zeros(1,numBall);
sin_sum=0;
cos_sum=0;
trial_span=10000;
trial_result=zeros(1,trial_span);
for nz=1:numBall  %calculate the average angle
    temp=X(:,:,nz)*[0;0;1];
    sin_sum=sin_sum+temp(2);
    cos_sum=cos_sum+temp(1);
end
sin_av=sin_sum./numBall;
cos_av=cos_sum./numBall;
for i=1:trial_span
    trial=trial+StepLength*i;
    for nz=1:numBall
        Xtemp(:,:,nz)=RotateZ(trial/rads(nz))*X(:,:,nz);
        temp=Xtemp(:,:,nz)*[0;0;1];
        theta2(nz)=atan2(real(temp(2))-sin_av,real(temp(1))-cos_av).^2;
    end
    trial_result(i)=sum(abs(theta2));
end
[~,I]=min(trial_result);
Zturn=StepLength*I;
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
%TODO:  write RotateQuaternionX(theta) 
end
