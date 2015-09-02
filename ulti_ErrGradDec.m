function ulti_ErrGradDec(ang_X,ang_Y,nr)
%Michael Williams 2015, Email:michael.williams.hy@gmail.com
%I use the gradient decrease method to find the proper step lenth for the X
%axis rotation and the Y axis rotation
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
    X(:,:,n)=RotateX(ang_X/rads(n))*RotateY(ang_Y/rads(n))*X(:,:,n);
    zaxis=X(:,:,n)*[0;0;1];
    psi(n) = acos(zaxis(3));
end
k=0;
Zturn=pi*600/nr; %nr denote how many steps you need to rotate 
% to turn the bigest ball around
nn=3*nr;
for i=1:nn
    for n=1:numBall %if the error cannot decrease any more, we rotate all around Z axis
        X(:,:,n)=RotateZ(Zturn/rads(n))*X(:,:,n);
        zaxis=X(:,:,n)*[0;0;1];
        psi(n) = acos(zaxis(3));
    end
    err_delta=sum(abs(psi));
    while err_delta>0.0001
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
            X(:,:,n)=RotateX(Xturn/rads(n))*X(:,:,n);
            zaxis=X(:,:,n)*[0;0;1];
            psi(n)=acos(zaxis(3));
        end
        err_new=sum(abs(psi));
        error_rec(k)=err_new;
        stp(k)=k;
        k=k+1;
        Yturn=GDFindLenthY(X,rads,numBall);%find a proper lenth for Y rotation
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
        err_new=sum(abs(psi));
        err_delta=err_old-err_new;
        error_rec(k)=err_new;
        stp(k)=k;
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
        hold off
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
function Xturn=GDFindLenthX(X,rads,numBall)
%find a proper rotation angle for X axis use gradient descent method with a
%flexibale step lenth
oldx=-1;
newx=0;
temp=newx;
precision=0.0001;
syms theta;
while abs(newx-oldx)>precision
    f_deriv=0;
    oldx=temp;
    for n=1:numBall %make the objective function adjustable to the number of balls.
        ff=acos(RotateX(theta/rads(n))*X(:,:,n)*[0;0;1]);
        %calculate the gradient`s sum of the objective function
        f_deriv=f_deriv+vpa(subs(ff,theta,newx)); 
    end
    f_deriv=f_deriv(3);
    gamma=0.01;
%     FindGammaX(f_deriv,X,numBall,newx,rads);
    temp=newx;
    newx=newx-f_deriv*gamma;
end
Xturn=newx;
end

% function gamma=FindGammaX(f_deriv,X,numball,newVar,rads)
% %I want to implement the 0.618 method, because I cannot add the function
% %handle in matlab
% a=0;
% b=2;
% GR=(sqrt(5)-1)/2; %global varible for Golden Section Search
% d=GR*(b-a)+a;
% c=b-GR*(b-a);
% while abs(c-d)>0.0000001
%     ffc=0;
%     ffd=0;
%     for n=1:numball
%         temp=RotateX(newVar+f_deriv*c/rads(n))*X(:,:,n)*[0;0;1];
%         ff=acos(temp);
%         ffc=ffc+ff(3); %get the objective function value for later comparison
%         temp=RotateX(newVar+f_deriv*d/rads(n))*X(:,:,n)*[0;0;1];
%         ff=acos(temp);
%         ffd=ffd+ff(3);
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
% 
% function Yturn=GDFindLenthY(X,rads,numBall)
% %find a proper rotation angle for Y axis
% oldy=-1;
% newy=0;
% temp=newy;
% precision=0.0001;
% syms theta;
% while abs(newy-oldy)>precision
%     f_deriv=0;
%     oldy=temp;
%     for n=1:numBall
%         ff=acos(RotateY(theta/rads(n))*X(:,:,n)*[0;0;1]);
%         f_deriv=f_deriv+vpa(subs(ff,theta,newy));
%     end
%     f_deriv=f_deriv(3);
%     gamma=0.01;
% %     FindGammaY(f_deriv,X,numBall,newy,rads);
%     temp=newy;
%     newy=newy-f_deriv*gamma;
% end
% Yturn=newy;
% end
% function gamma=FindGammaY(f_deriv,X,numball,newVar,rads)
% %I want to implement the 0.618 method, because I cannot add the function
% %handle in matlab
% GR=(sqrt(5)-1)/2; %global varible for Golden Section Search
% a=0;
% b=2;
% d=GR*(b-a)+a;
% c=b-GR*(b-a);
% while abs(c-d)>0.0000001
%     ffc=0;
%     ffd=0;
%     for n=1:numball
%         temp=RotateY(newVar+f_deriv*c/rads(n))*X(:,:,n)*[0;0;1];
%         ff=acos(temp);
%         ffc=ffc+ff(3); %get the objective function value for later comparison
%         temp=RotateY(newVar+f_deriv*d/rads(n))*X(:,:,n)*[0;0;1];
%         ff=acos(temp);
%         ffd=ffd+ff(3);
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
        
        