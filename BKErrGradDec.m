function BKErrGradDec(ang_X,ang_Y,nr)
%Michael Williams 2015, Email:michael.williams.hy@gmail.com
%I use the gradient decrease method to find the proper step lenth for the X
%axis rotation and the Y axis rotation
%=================================NOTE====================================
%this one is cheating because I manipulated the balls seperately to make
%them align together, with this test we can see if aligning them together
%is a good idea.
%IT¡¡DOES NOT WORK!!!
%=================================NOTE====================================
format compact
clc
tic
numStep=100000;  %prelocating of numStep to increase the speed
rads = [10,13,17,29];
numBall=numel(rads);
stp=zeros(1,numStep);
path=zeros(2,numStep);
error_rec=zeros(1,numStep);
X = repmat(eye(3),[1,1,numBall]);
psi = zeros(numBall,1);
for n=1:numBall %set the initial angle of balls
    qrot=angle2quat(0,ang_Y/rads(n),ang_X/rads(n));
    X(:,:,n)=quatrotate(qrot,X(:,:,n));
    zaxis=X(:,:,n)*[0;0;1];
    psi(n) = acos(zaxis(3));
end
k=0;
% Zturn=pi*600/nr; %nr denote how many steps you need to rotate 
% to turn the bigest ball around
nn=3*nr;
cheat_theta=zeros(1,numBall);
for i=1:nn
    for n=1:numBall
        Z_orient=X(:,:,n)*[0;0;1];
        cheat_theta(n)=atan2(real(Z_orient(2)),real(Z_orient(1)));
        qrot=angle2quat(cheat_theta(n),0,0);
        X(:,:,n)=quatrotate(qrot,X(:,:,n));
        psi=acos(Z_orient(3));
    end
    err_delta=sum(abs(psi));
    while err_delta>0.000001
        k=k+1;
        for n=1:numBall
            zaxis=X(:,:,n)*[0;0;1];
            psi(n) = acos(zaxis(3));
        end
        err_old=sum(abs(psi));
        Xturn=GDFindLenthX(X,rads,numBall);%find a proper lenth for X rotation
        if k==1
            path(2,k)=Xturn;
            path(1,k)=0;
        else
            path(2,k)=Xturn+path(2,k-1);
            path(1,k)=path(1,k-1);
        end
        for n=1:numBall
            Xturn=real(Xturn)/rads(n);
            qrot=angle2quat(0,0,Xturn);
            X(:,:,n)=quatrotate(qrot,X(:,:,n));
            zaxis=X(:,:,n)*[0;0;1];
            psi(n)=acos(zaxis(3));
        end
        err_new=sum(abs(psi));
        error_rec(k)=err_new;
        stp(k)=k;
        k=k+1;
        Yturn=GDFindLenthY(X,rads,numBall);%find a proper lenth for Y rotation
        for n=1:numBall
            Yturn=real(Yturn)/rads(n);
            qrot=angle2quat(0,Yturn,0);
            X(:,:,n)=quatrotate(qrot,X(:,:,n));
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
        err_new=sum(abs(psi));
        err_delta=abs(err_old-err_new);
        error_rec(k)=err_new;
        stp(k)=k;
        %the following just for drawing
        error_rec1=error_rec(1:k-1)*180/pi;
        figure(1);
        plot(stp(1:k-1),error_rec1);
        title('Noiseless Ensemble Control of 4 Spheres Orientation');
        xlabel('steps');
        ylabel('overall error(degs)');
        path1=path(:,1:k-1);
        save('GDmyData.mat','error_rec','path','path1');
        figure(2)
        plot(path1(1,:),path1(2,:));
        title('movement of the panel for the control process');
        xlabel('motion projected on the X axis')
        ylabel('motion projected on the Y axis')
    end
end
%we need to get the error_rec rid of zero, in order to make a plot
% we do the following steps:
stp1=stp(1:k-1);
error_rec1=error_rec(1:k-1)*180/pi;
figure(1);
plot(stp1,error_rec1);
title('Noiseless Ensemble Control of 4 Spheres Orientation');
xlabel('steps');
ylabel('overall error(degs)');
path1=path(:,1:k-1);
save('GDmyData.mat','error_rec','path','path1');
figure(2)
plot(path1(1,:),path1(2,:));
title('movement of the panel for the control process');
xlabel('motion projected on the X axis')
ylabel('motion projected on the Y axis')
toc
end
%Because the quaternions manipulation function in MATLAB cannot take in the
%symbolic parameters, we still need to use the traditional rotation matrix
%in the following function.
function Xturn=GDFindLenthX(X,rads,numBall)
%find a proper rotation angle for X axis use gradient descent method with a
%flexibale step lenth
oldx=-1;
newx=0;
precision=0.000001;
ff=cell(1,numBall);
syms theta
for n=1:numBall
    temp=RotateX(theta/rads(n))*X(:,:,n);
    ff{n}=acos(temp*[0;0;1]);
end
while abs(newx-oldx)>precision
    oldx=newx;
    f_deriv=0;
    for n=1:numBall
        temp=vpa(subs(ff{n},theta,newx));
        f_deriv=f_deriv+temp(3);
    end
    gamma=FindGammaX(f_deriv,X,numBall,newx,rads);
    newx=newx-f_deriv*gamma;
end
Xturn=newx;
end

function gamma=FindGammaX(f_deriv,X,numball,newVar,rads)
%I want to implement the 0.618 method, because I cannot add the function
%handle in matlab
a=0;
b=2;
GR=(sqrt(5)-1)/2; %global varible for Golden Section Search
d=GR*(b-a)+a;
c=b-GR*(b-a);
while abs(c-d)>0.0000000001
    ffc=0;
    ffd=0;
    for n=1:numball
        temp1=RotateX((newVar-f_deriv*c)/rads(n))*X(:,:,n);
        temp=temp1*[0;0;1];
        obj_temp=acos(temp(3));
        ffc=ffc+obj_temp; %get the objective function value for later comparison
        temp1=RotateX((newVar-f_deriv*d)/rads(n))*X(:,:,n);
        temp=temp1*[0;0;1];
        obj_temp=acos(temp(3));
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

function Yturn=GDFindLenthY(X,rads,numBall)
%find a proper rotation angle for Y axis
oldy=-1;
newy=0;
precision=0.000001;
ff=cell(1,numBall);
syms theta
for n=1:numBall
    temp=RotateY(theta/rads(n))*X(:,:,n);
    ff{n}=acos(temp*[0;0;1]);
end
while abs(newy-oldy)>precision
    oldy=newy;
    f_deriv=0;
    for n=1:numBall
        temp=vpa(subs(ff{n},theta,newy));
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
    for n=1:numball
        temp1=RotateY((newVar-f_deriv*c)/rads(n))*X(:,:,n);
        temp=temp1*[0;0;1];
        obj_temp=acos(temp(3));
        ffc=ffc+obj_temp; %get the objective function value for later comparison
        temp1=RotateY((newVar-f_deriv*d)/rads(n))*X(:,:,n);
        temp=temp1*[0;0;1];
        obj_temp=acos(temp(3));
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
% function RzTh = RotateZ(theta)
% RzTh = [cos(theta),  -sin(theta),0;
%     sin(theta),   cos(theta),0;
%     0,  0,  1];
% end
        
        