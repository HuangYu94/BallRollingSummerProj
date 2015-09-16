function AdjErrGradDec5(ang_X,ang_Y)
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
% This file align the balls` z axises into the positive direction of the X-
% axis, so after each rotation about Z axis, moving toward the negative
% direction of the X axis will always be the best choice.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format compact
syms theta 
theta=5;  % this is needed for optimization
if nargin<2
    ang_X = 7*pi;ang_Y = 9*pi;
end

clc
tic
numStep=10000;  %prelocation of numStep to increase the speed
rads = [10,13,17,19];
numBall=numel(rads);
stp=zeros(1,numStep);
path=zeros(2,numStep); %record the trajectory
%path(1,:)for recording X
%path(2,:)for recording Y
error_rec=zeros(1,numStep);
ZrotRec=zeros(3,numStep);%record the position for z rotation and the rotation angle
RMSE=zeros(1,numStep);
RMSEinDEG=zeros(1,numStep);
global precision;
precision=0.00001;
X = repmat(eye(3),[1,1,numBall]);
psi = zeros(numBall,1);
for n=1:numBall %set the initial angle of balls
    X(:,:,n)=RotateX(ang_X/rads(n))*RotateY(ang_Y/rads(n));
    zaxis=X(:,:,n)*[0;0;1];
    psi(n) = acos(zaxis(3));
end
k=0; % number of steps
Zround=0; %index to record the Z angle;
% Zturn=pi*600/nr; %nr denote how many steps you need to rotate
% to turn the bigest ball around
precision_control=0.01;  %stop the process when this error is reached
err_new=10;%make the loop start
while err_new>precision_control
    for n=1:numBall  %calculate the initial error
        Z_orient=X(:,:,n)*[0;0;1];
        psi=acos(Z_orient(3));
    end
    err_delta=sum(psi.^2);
    while err_delta>0.0005
        k=k+1;
        for n=1:numBall %calculate the current error
            zaxis=X(:,:,n)*[0;0;1];
            psi(n) = acos(zaxis(3));
        end
        err_old=sqrt(sum(psi.^2)/numBall);
        Xturn=GDFindLengthX(X,rads,numBall);%find best length for X rotation
        path(2,k+1)=Xturn+path(2,k); %record the trajectory
        path(1,k+1)=path(1,k);
        for n=1:numBall  %apply the rotation.   %GOAL: use quaternions
            X(:,:,n)=RotateX(Xturn/rads(n))*X(:,:,n);
            zaxis=X(:,:,n)*[0;0;1];
            psi(n)=acos(zaxis(3));
        end
        err_new=sqrt(sum(psi.^2)/numBall);
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
        path(1,k+1)=Yturn+path(1,k);
        path(2,k+1)=path(2,k);
        err_new=sqrt(sum(psi.^2)/numBall);
        err_delta=abs(err_old-err_new);
        error_rec(k)=err_new;
        
        RMSE(k) = sqrt(1/numBall*sum(psi.^2));
        RMSEinDEG(k) = RMSE(k)*180/pi;
        stp(k)=k;
        %the following just for drawing
        %error_recDEGREE=error_rec(1:k-1)*180/pi;
        figure(1);
        plot(stp(1:k),RMSEinDEG(1:k));
        title('Noiseless Ensemble Control of 4 Spheres Orientation');
        xlabel('steps (x,y, or z rotations)');
        ylabel('standard deviation from the Z-Axis of the WOLRD coordinates (degs)');  %  sqrt(1/numBall*sum( psi.^2))  
        path1=path(:,1:k);
        if Zround>=2
            figure(2)
            plot(path1(1,:),path1(2,:), path1(1,1),path1(2,1),'go',path1(1,end),path1(2,end),'ro',ZrotRec(1,1:Zround),ZrotRec(2,1:Zround),'bo');
            legend ('path','start','end','rotate around z','Location','Southeastoutside');
            title('movement of the panel for the control process');
            xlabel('motion projected on the X axis')
            ylabel('motion projected on the Y axis')
        else
            break
        end
    end
    Zround=Zround+1;
    Zturn=Z_align(X,rads);  %find out the optimal angle to align the balls`
    % angle theta together.
    ZrotRec(1,Zround)=path(1,k);
    ZrotRec(2,Zround)=path(2,k);
    ZrotRec(3,Zround)=Zturn;
    for n=1:numBall
        X(:,:,n)=RotateZ(Zturn/rads(n))*X(:,:,n);
        Z_orient=X(:,:,n)*[0;0;1];
        psi=acos(Z_orient(3));
    end
end
%we need to get the error_rec rid of zero, in order to make a plot
% we do the following steps:
%error_recDEGREE=error_rec(1:k-1)*180/pi; %convert from radians to degrees
figure(1);
plot(stp(1:k),RMSEinDEG(1:k));
title('Noiseless Ensemble Control of 4 Spheres Orientation');
xlabel('steps');
ylabel('overall error(degs)');
path1=path(:,1:k);
save('GDmyData1.mat','error_rec','path','path1','ZrotRec','X');
figure(2)
plot(path1(1,:),path1(2,:), path1(1,1),path1(2,1),'go',path1(1,end),path1(2,end),'rx' );
hold on
plot(ZrotRec(1,1:Zround),ZrotRec(2,1:Zround),'ro');
title('movement of the panel for the control process');
legend ('path','start','end','rotate around z','location','Northeastoutside')
xlabel('motion projected on the X axis')
ylabel('motion projected on the Y axis')
toc

    function Xturn=GDFindLengthX(X,rads,numBall)
        %find a proper rotation angle for X axis use gradient descent method with a
        %flexible step length
        % INPUTS:
        % X: Rotation matrices for each ball
        % rads: array recording the radius for each ball
        % numBall: number of balls
        ff=cell(1,numBall);
        newx=MakeUnimodalX(X,rads); %find a unimodal region for the gradient descent method
        oldx=newx-1;  % make the loop start
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

    function startpoint=MakeUnimodalX(X,rads)
        %Previously, I implemented the blindly trial method which seems effective
        %in solving the problem, but it has an obvious defect that is too mush
        %unnecessary computation. However, now, I use it to find a proper
        %startpoint for the gradient descent method, therefore, I may not search
        %so intensively as the former one.
        %I guess that there will be no more than 1 local extremum within the period
        %of Max(rads)*pi/2 which can be used as the maximum step length. However,
        %for better performance, I want to use smaller one.
        % INPUTS:
        % X: Rotation matrices for each ball
        % rads: array recording the radius for each ball
        %Yu Huang 2015, e-mail:Michael.Williams.hy@gmail.com
        MAXrad=max(rads);
        step_length=MAXrad*pi/10;
        trial0=-2*pi*MAXrad;
        Mstep=41; %make the maximum trial value at pi*Mrad
        tempPsi=zeros(1,numBall); %store the psi.^2 for each ball
        errorX=zeros(1,Mstep);
        for i=1:Mstep+1
            trial=trial0+(i-1)*step_length;
            for nx=1:numBall
                temp=RotateX(trial/rads(nx))*X(:,:,nx); %effect of symbolic rotation
                temp1=acos(temp*[0;0;1]);
                tempPsi(nx)=temp1(3);
            end
            errorX(i)=sum(tempPsi.^2);
        end
        [~,minI]=min(errorX);
        startpoint=trial0+minI*step_length;
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
        while abs(c-d)>0.000000000001
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
        ff=cell(1,numBall);
        newy=MakeUnimodalY(X,rads);
        oldy=newy-1;
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

    function startpoint=MakeUnimodalY(X,rads)
        %Previously, I implemented the blindly trial method which seems effective
        %in solving the problem, but it has an obvious defect that is too mush
        %unnecessary computation. However, now, I use it to find a proper
        %startpoint for the gradient descent method, therefore, I may not search
        %so intensively as the former one.
        %I guess that there will be no more than 1 local extremum within the period
        %of Max(rads)*pi/2 which can be used as the maximum step length. However,
        %for better performance, I want to use smaller one.
        % INPUTS:
        % X: Rotation matrices for each ball
        % rads: array recording the radius for each ball
        %Yu Huang 2015, e-mail:Michael.Williams.hy@gmail.com
        MAXrad=max(rads);
        step_length=MAXrad*pi/10;
        trial0=-2*pi*MAXrad;
        Mstep=41; %make the maximum trial value at pi*Mrad
        tempPsi=zeros(1,numBall); %store the psi.^2 for each ball
        errorY=zeros(1,Mstep);
        for i=1:Mstep+1
            trial=trial0+(i-1)*step_length;
            for nx=1:numBall
                temp=RotateY(trial/rads(nx))*X(:,:,nx); %effect of symbolic rotation
                temp1=acos(temp*[0;0;1]);
                tempPsi(nx)=temp1(3);
            end
            errorY(i)=sum(tempPsi.^2);
        end
        [~,minI]=min(errorY);
        startpoint=trial0+minI*step_length;
    end
        
        function gamma=FindGammaY(f_deriv,X,numball,newVar,rads)
            %I want to implement the 0.618 method, because I cannot add the function
            %handle in matlab
            GR=(sqrt(5)-1)/2; %global varible for Golden Section Search
            a=0;
            b=2;
            d=GR*(b-a)+a;
            c=b-GR*(b-a);
            while abs(c-d)>0.000000000001
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
    function Zturn=Z_align(X,rads)
        % this function was written for finding the optimal angle for rotation
        % around Z axsis to make the balls align together to the positive direction
        % of X axsis.
        %I still extend the previous
        %Michael Williams 2015, e-mail:michael.williams.hy@gamil.com
        %         sin_sum=0;
        %         cos_sum=0;
        %         for Balln=1:numBall  %calculate the average angle
        %             temp=X(:,:,Balln)*[0;0;1];
        %             sin_sum=sin_sum+temp(2);
        %             cos_sum=cos_sum+temp(1);
        %         end
        %         sin_av=sin_sum./numBall;
        %         cos_av=cos_sum./numBall;
        obj_f2=0;
        obj_f1=0;
        for Balln=1:numBall   %calculate the objctive funtion for every ball and store them in obj_f
            temp=RotateZ(theta/rads(Balln))*X(:,:,Balln);
            temp1=temp*[0;0;1];
            %I align them to the so-called "avrage direction" direction
            obj_f2=obj_f2+atan2(temp1(2),temp(1)).^2;
            obj_f1=obj_f1+atan2(temp1(2),temp(1));
        end
        obj_f1=obj_f1/numBall;
        obj_f2=obj_f2/numBall;
        obj_function=obj_f1-obj_f2.^2;
        %implement the gradient descent method
        newz=MakeUnimodalZ(X,rads);
        oldz=newz-1;   %start the loop
        while abs(newz-oldz)>precision
            oldz=newz;
            f_deriv1=vpa(subs(obj_function,theta,newz));
            gamma=FindGammaZ(f_deriv1,X,numBall,newz,rads);
            newz=newz-f_deriv1*gamma;
        end
        Zturn=newz;
    end

    function gamma=FindGammaZ(f_deriv,X,numBall,newz,rads)
        %still I use the golden section search to find the proper gamma
        %         sin_sum=0;
        %         cos_sum=0;
        %         for Balln=1:numBall  %calculate the average angle
        %             temp=X(:,:,Balln)*[0;0;1];
        %             sin_sum=sin_sum+temp(2);
        %             cos_sum=cos_sum+temp(1);
        %         end
        %         sin_av=sin_sum./numBall;
        %         cos_av=cos_sum./numBall;
        GR=(sqrt(5)-1)/2; %global varible for Golden Section Search
        a=0;
        b=2;
        d=GR*(b-a)+a;
        c=b-GR*(b-a);
        obj_f2=0;
        obj_f1=0;
        for Balln=1:numBall   %calculate the objctive funtion for every ball and store them in obj_f
            temp=RotateZ(theta/rads(Balln))*X(:,:,Balln);
            temp1=temp*[0;0;1];
            %I align them to the so-called "avrage direction" direction
            obj_f2=obj_f2+atan2(temp1(2),temp(1)).^2;
            obj_f1=obj_f1+atan2(temp1(2),temp(1));
        end
        obj_f1=obj_f1/numBall;
        obj_f2=obj_f2/numBall;
        obj_function=sqrt(obj_f1-obj_f2.^2);
        while abs(c-d)>0.000000000001
            ffc=subs(obj_function,theta,newz-f_deriv*c);
            ffd=subs(obj_function,theta,newz-f_deriv*d);
            %             for Balln=1:numBall
            %                 temp1=RotateZ((newz-f_deriv*c)/rads(Balln))*X(:,:,Balln);
            %                 temp=temp1*[0;0;1];
            %                 obj_temp=atan2(real(temp(2))-sin_av,real(temp(1))-cos_av).^2;
            %                 ffc=ffc+obj_temp; %get the objective function value for later comparison
            %                 temp1=RotateZ((newz-f_deriv*d)/rads(Balln))*X(:,:,Balln);
            %                 temp=temp1*[0;0;1];
            %                 obj_temp=atan2(real(temp(2))-sin_av,real(temp(1))-cos_av).^2;
            %                 ffd=ffd+obj_temp;
            %             end
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
    function Zturn=MakeUnimodalZ(X,rads)
        %for this function I still use the blindly trial method to find the
        %alignment angle.
        %Still we can use the gradient descent method to search for a
        %more precise solution to the rotation angle around Z axis
        %Michael Williams 2015, e-mail:michael.williams.hy@gamil.com
        MAXrad=max(rads);
        StepLength=MAXrad*pi/10; %length for trial
        trial0=-4*pi*MAXrad;  %initialize
        %             Xtemp=repmat(eye(3),[1,1,numBall]);
        %             theta2=zeros(1,numBall);
        trial_span=81;
        trial_result=zeros(1,trial_span);
        %             sin_sum=0;
        %             cos_sum=0;
        %             for nz=1:numBall  %calculate the average angle
        %                 temp=X(:,:,nz)*[0;0;1];
        %                 sin_sum=sin_sum+temp(2);
        %                 cos_sum=cos_sum+temp(1);
        %             end
        %             sin_av=sin_sum./numBall;
        %             cos_av=cos_sum./numBall;
        obj_f2=0;
        obj_f1=0;
        for Balln=1:numBall   %calculate the objctive funtion for every ball and store them in obj_f
            temp=RotateZ(theta/rads(Balln))*X(:,:,Balln);
            temp1=temp*[0;0;1];
            %I align them to the so-called "avrage direction" direction
            obj_f2=obj_f2+atan2(temp1(2),temp(1)).^2;
            obj_f1=obj_f1+atan2(temp1(2),temp(1));
        end
        obj_f1=obj_f1/numBall;
        obj_f2=obj_f2/numBall;
        obj_function=sqrt(obj_f1-obj_f2.^2);
        for i=1:trial_span
            trial=trial0+StepLength*i;
            trial_result(i)=subs(obj_function,theta,trial);
            %                 for nz=1:numBall
            %                     Xtemp(:,:,nz)=RotateZ(trial/rads(nz))*X(:,:,nz);
            %                     temp=Xtemp(:,:,nz)*[0;0;1];
            %                     theta2(nz)=atan2(real(temp(2))-sin_av,real(temp(1))-cos_av).^2;
            %                 end
            %                 trial_result(i)=sum(abs(theta2));
        end
        [~,I]=min(trial_result);
        Zturn=StepLength*I+trial0;
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
