function AdjErrGradDec6(ang_X,ang_Y)
%Michael Williams 2015, Email: michael.williams.hy@gmail.com
%
%The gradient descent method uses a flexible step length
% for X-axis rotation and the Y axis rotation.
% Available at the github at:  git@github.com:HuangYu94/BallRollingSummerProj.git

% profile the code-- what is slow?
% error should monotonically decrease why doesn't it?  -- we will switch to
% plot the RMSE
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
rads = [10,11,13,15];
numBall=numel(rads);
RadAv=sum(rads)/numBall;
stp=zeros(1,numStep);
path=zeros(2,numStep); %record the trajectory of X-Y rolling
%path(1,:)for recording X
%path(2,:)for recording Y
error_rec=zeros(1,numStep);
ZrotRec=zeros(3,numStep);%record the position for z rotation and the rotation angle
error_sep=zeros(numBall,numStep);
RMSE=zeros(1,numStep);
RMSEinDEG=zeros(1,numStep);
global precision;
precision=0.00001;
global color_arr;
color_arr=['y','m','c','r','g','b','w','k'];
X = repmat([0;0;1],[1,1,numBall]);
psi = zeros(numBall,1);
for n=1:numBall %set the initial angle of balls
    X(:,:,n)=RotateX(ang_X/rads(n))*RotateY(ang_Y/rads(n))*X(:,:,n);
%     zaxis=X(:,:,n)*[0;0;1];
    psi(n) = acos(X(3,:,n));
end
k=0; % number of steps
Zround=0; %index to record the Z angle;
% Zturn=pi*600/nr; %nr denote how many steps you need to rotate
% to turn the bigest ball around
precision_control=0.01;  %stop the process when this error is reached
err_new=10;%make the loop start
while err_new>precision_control && k<=1000
    psi=zeros(1,numBall);
    for n=1:numBall  %calculate the initial error
%         Z_orient=X(:,:,n)*[0;0;1];
        psi(n)=acos(X(3,:,n));
    end
    err_delta=sum(psi.^2);
    while err_delta>0.00001
        k=k+1;
        for n=1:numBall %calculate the current error
%             zaxis=X(:,:,n)*[0;0;1];
            psi(n) = acos(X(3,:,n));
            error_sep(n,k)=psi(n);
        end
        err_old=sqrt(sum(psi.^2)/numBall);
        Xturn=GDFindLengthX(X,rads,numBall);%find best length for X rotation
        path(2,k+1)=Xturn+path(2,k); %record the trajectory
        path(1,k+1)=path(1,k);
        for n=1:numBall  %apply the rotation.   %GOAL: use quaternions
            X(:,:,n)=RotateX(Xturn/rads(n))*X(:,:,n);
%             zaxis=X(:,:,n)*[0;0;1];
            psi(n)=acos(X(3,:,n));
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
%             zaxis=X(:,:,n)*[0;0;1];
            psi(n)=acos(X(3,:,n));
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
        hold on
        legendinfo=cell(1,numBall);
        legendinfo{1}='error sum';
        for ii=1:numBall
            plot(stp,180/pi*error_sep(ii,:),color_arr(ii));
            legendinfo{ii+1}=['rad:',num2str(rads(ii))];
        end
        legend(legendinfo);
        legend('boxoff');
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
    SumTheta1=0;
    SumTheta2=0;
    for n=1:numBall %calculate the current error
        %             zaxis=X(:,:,n)*[0;0;1];
        psi(n) = acos(X(3,:,n));
        SumTheta1=SumTheta1+X(1,:,n);
        SumTheta2=SumTheta2+X(2,:,n);
    end
    OrientTheta=atan2(SumTheta2,SumTheta1);
    PsiAv=sum(psi)/numBall;  %calculate the average psi error
    rctemp=RadAv*tan(PsiAv/2);
    if rctemp>20   %practical restriction
        rc=20;
    else
        rc=rctemp;
    end
    RotPoint=[-rc*cos(OrientTheta),-rc*sin(OrientTheta),0];
    alpha=FindAlpha(RotPoint,X,rads,rc,PsiAv);
    BallConfig=repmat([0,0;0,0;0,0],[1,1,numBall]);
    for n=1:numBall
        rcRoll=RadAv*rads(n)/sqrt(RadAv.^2+rads(n).^2);
        BallConfig(:,2,n)=X(:,:,n)*rads(n);
        BallConfig(:,:,n)=ArbAxisRotate(RotPoint,[0,0,rads(n)],-alpha*rc/rcRoll,BallConfig(:,:,n));
        BallConfig(:,:,n)=RotateZ(alpha)*BallConfig(:,:,n);
        X(:,:,n)=(BallConfig(:,2,n)-BallConfig(:,1,n))/rads(n);
    end
    % angle theta together.
    ZrotRec(1,Zround)=path(1,k);
    ZrotRec(2,Zround)=path(2,k);
    ZrotRec(3,Zround)=alpha;
end
%we need to get the error_rec rid of zero, in order to make a plot
% we do the following steps:
%error_recDEGREE=error_rec(1:k-1)*180/pi; %convert from radians to degrees
figure(1);
plot(stp(1:k),RMSEinDEG(1:k));
title('Noiseless Ensemble Control of 4 Spheres Orientation');
xlabel('steps');
ylabel('overall error(degs)');
hold on
legendinfo=cell(1,numBall);
legendinfo{1}='error sum';
for ii=1:numBall
    plot(stp,180/pi*error_sep(ii,:),color_arr(ii));
    legendinfo{ii+1}=['rad:',num2str(rads(ii))];
end
legend(legendinfo);
legend('boxoff');
path1=path(:,1:k);
save('GDmyData1.mat','error_rec','path','path1','ZrotRec','X');
figure(2)
plot(path1(1,:),path1(2,:), path1(1,1),path1(2,1),'go',path1(1,end),path1(2,end),'rx' );
hold on
plot(ZrotRec(1,1:Zround),ZrotRec(2,1:Zround),'ro');
title('movement of the panel for the control process');
legend ('path','start','end','rotate around z','location','Northeastoutside')
legend('boxoff')
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
            ff{nx}=acos(temp(3)).^2;  %psi as a function of theta.
        end
        while abs(newx-oldx)>precision
            oldx=newx;
            f_deriv=0;
            for nx=1:numBall
                temp=vpa(subs(ff{nx},theta,newx));
                f_deriv=f_deriv+temp;
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
%                 temp1=acos(temp(3));
                tempPsi(nx)=acos(temp(3));
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
                temp=RotateX((newVar-f_deriv*c)/rads(ng))*X(:,:,ng);
%                 temp=temp1*[0;0;1];
                obj_temp=acos(temp(3)).^2;
                ffc=ffc+obj_temp; %get the objective function value for later comparison
                temp=RotateX((newVar-f_deriv*d)/rads(ng))*X(:,:,ng);
%                 temp=temp*[0;0;1];
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
            ff{ny}=acos(temp(3)).^2;
        end
        while abs(newy-oldy)>precision
            oldy=newy;
            f_deriv=0;
            for ny=1:numBall
                temp=vpa(subs(ff{ny},theta,newy));
                f_deriv=f_deriv+temp;
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
%                 temp1=acos(temp*[0;0;1]);
                tempPsi(nx)=acos(temp(3));
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
                    temp=RotateY((newVar-f_deriv*c)/rads(ng))*X(:,:,ng);
                    obj_temp=acos(temp(3)).^2;
                    ffc=ffc+obj_temp; %get the objective function value for later comparison
                    temp=RotateY((newVar-f_deriv*d)/rads(ng))*X(:,:,ng);
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
    
    function alpha=FindAlpha(RotPoint,X,rads,rc,psi_ang)
        %This function will determine for what angle should the spheres
        %rotate, the metric for determining whether a turn is good or not
        %is that good turning angle always make the Z axis of the balls
        %align together. Additionally, I observed that when the average psi
        %error decreased to a small value, the metric we use here can`t
        %make a wise choice about alpha, so I want this function to
        %generate a really big value for the circular motion, because when
        %the error is small, the lollipop rotation and the rotation about
        %the Z axis caused by the translation will somewhat neutralize to
        %each other.
        %Yu Huang 2015, E-mail: Michael.Williams.hy@gmail.com
        error_in_deg=180*psi_ang/pi;
        if error_in_deg<=5
            alpha=60*pi;
        else
            configX=repmat([0,0;0,0,;0,0],[1,1,numBall]);
            TryTime=10000;
            trial_result=zeros(1,TryTime);
            %those two values are for storing the information for calculating
            %the Theta:
            temptheta1=zeros(1,numBall);
            temptheta2=zeros(1,numBall);
            for numStp=1:TryTime
                tryAlpha=10*pi+numStp*40*pi/TryTime;
                for numLoop=1:numBall
                    configX(:,2,numLoop)=X(:,:,numLoop)*rads(numLoop);
                    rcBall=rc*rads(numLoop)/sqrt(rc.^2+rads(numLoop).^2);
                    configX(:,:,numLoop)=ArbAxisRotate(RotPoint,[0,0,rads(numLoop)],-tryAlpha*rc/rcBall,configX(:,:,numLoop));
                    configX(:,:,numLoop)=RotateZ(tryAlpha)*configX(:,:,numLoop);
                    Zaxis=(configX(:,2,numLoop)-configX(:,1,numLoop))/rads(numLoop);
                    temptheta1(numLoop)=Zaxis(1);
                    temptheta2(numLoop)=Zaxis(2);
                end
                theta_av=atan2(sum(temptheta2),sum(temptheta1));
                alpha_metric=0;
                for numB=1:numBall
                    alpha_metric=alpha_metric+(atan2(temptheta2(numB),temptheta1(numB))-theta_av).^2;
                end
                trial_result(numStp)=alpha_metric;
            end
            [~,IN]=min(trial_result);
            alpha=10*pi*IN*40*pi/TryTime;
        end
    end
    

%     function Zturn=FindAlpha(RotPoint,X,rads,rc)
%         % this function was written for finding the optimal angle for rotation
%         % around Z axsis to make the balls align together to the positive direction
%         % of X axsis.
%         %I still extend the previous
%         %Michael Williams 2015, e-mail:michael.williams.hy@gamil.com
%         %         for Balln=1:numBall  %calculate the average angle
%         %             temp=X(:,:,Balln)*[0;0;1];
%         %             sin_sum=sin_sum+temp(2);
%         %             cos_sum=cos_sum+temp(1);
%         %         end
%         sin_sum=0;
%         cos_sum=0;        
%         tempX=repmat([0,0;0,0;0,1],[1,1,numBall]);
%         obj_function=0;
%         for Balln=1:numBall   %calculate the objctive funtion for every ball and store them in obj_f
%             tempX(:,2,Balln)=X(:,:,Balln)*rads(Balln);
%             r_locus=rads(Balln)*rc/sqrt(rads(Balln).^2+rc.^2);
%             tempX(:,:,Balln)=ArbAxisRotate(RotPoint,[0,0,rads(Balln)],-theta*rc/r_locus,tempX(:,:,Balln));
%             tempX(:,:,Balln)=RotateZ(theta)*tempX(:,:,Balln);
%             zaxis=(tempX(:,:,Balln)-tempX(:,:,Balln))/rads(Balln);
%             sin_sum=sin_sum+zaxis(2);
%             cos_sum=cos_sum+zaxis(1);
%         end
%         Theta_Av=atan2(sin_sum,cos_sum);
%         for Balln=1:numBall
%             zaxis=(tempX(:,:,Balln)-tempX(:,:,Balln))/rads(Balln);
%             theta_ori=(atan2(zaxis(2),zaxis(1))-Theta_Av).^2;
%             obj_function=obj_function+theta_ori;
%         end
%         %implement the gradient descent method
%         newz=MakeUnimodalZ(RotPoint,X,rads,rc);
%         oldz=newz-1;   %start the loop
%         c=0;
%         while abs(newz-oldz)>precision && c<20
%             oldz=newz;
%             f_deriv1=vpa(subs(obj_function,theta,newz));
%             gamma=FindGammaZ(f_deriv1,X,numBall,newz,rads,rc);
%             newz=newz-f_deriv1*gamma;
%             c=c+1;
%         end
%         Zturn=newz;
%     end

%     function gamma=FindGammaZ(f_deriv,X,numBall,newz,rads,rc)
%         %still I use the golden section search to find the proper gamma
%         %         sin_sum=0;
%         %         cos_sum=0;
%         %         for Balln=1:numBall  %calculate the average angle
%         %             temp=X(:,:,Balln)*[0;0;1];
%         %             sin_sum=sin_sum+temp(2);
%         %             cos_sum=cos_sum+temp(1);
%         %         end
%         %         sin_av=sin_sum./numBall;
%         %         cos_av=cos_sum./numBall;
%         GR=(sqrt(5)-1)/2; %global varible for Golden Section Search
%         a=0;
%         b=2;
%         d=GR*(b-a)+a;
%         c=b-GR*(b-a);
%         sin_sum=0;
%         cos_sum=0;        
%         tempX=repmat([0,0;0,0;0,1],[1,1,numBall]);
%         obj_function=0;
%         for Balln=1:numBall   %calculate the objctive funtion for every ball and store them in obj_f
%             tempX(:,2,Balln)=X(:,:,Balln)*rads(Balln);
%             r_locus=rads(Balln)*rc/sqrt(rads(Balln).^2+rc.^2);
%             tempX(:,:,Balln)=ArbAxisRotate(RotPoint,[0,0,rads(Balln)],-theta*rc/r_locus,tempX(:,:,Balln));
%             tempX(:,:,Balln)=RotateZ(theta)*tempX(:,:,Balln);
%             zaxis=(tempX(:,:,Balln)-tempX(:,:,Balln))/rads(Balln);
%             sin_sum=sin_sum+zaxis(2);
%             cos_sum=cos_sum+zaxis(1);
%         end
%         Theta_Av=atan2(sin_sum,cos_sum);
%         for Balln=1:numBall
%             zaxis=(tempX(:,:,Balln)-tempX(:,:,Balln))/rads(Balln);
%             theta_ori=(atan2(zaxis(2),zaxis(1))-Theta_Av).^2;
%             obj_function=obj_function+theta_ori;
%         end
%         while abs(c-d)>0.000000000001
%             ffc=subs(obj_function,theta,newz-f_deriv*c);
%             ffd=subs(obj_function,theta,newz-f_deriv*d);
%             %             for Balln=1:numBall
%             %                 temp1=RotateZ((newz-f_deriv*c)/rads(Balln))*X(:,:,Balln);
%             %                 temp=temp1*[0;0;1];
%             %                 obj_temp=atan2(real(temp(2))-sin_av,real(temp(1))-cos_av).^2;
%             %                 ffc=ffc+obj_temp; %get the objective function value for later comparison
%             %                 temp1=RotateZ((newz-f_deriv*d)/rads(Balln))*X(:,:,Balln);
%             %                 temp=temp1*[0;0;1];
%             %                 obj_temp=atan2(real(temp(2))-sin_av,real(temp(1))-cos_av).^2;
%             %                 ffd=ffd+obj_temp;
%             %             end
%             if ffc<=ffd
%                 b=d;
%                 d=c;
%                 c=b-GR*(b-a);
%             else
%                 a=c;
%                 c=d;
%                 d=a+GR*(b-a);
%             end
%         end
%         gamma=(a+b)/2;
%     end

%     function alpha=MakeUnimodalZ(RotPoint,X,rads,rc)
%         %this function will determine for what angle should the spheres
%         %rotate
%         %Yu Huang 2015, E-mail: Michael.Williams.hy@gmail.com
%         TryTime=10000;
%         trial_result=zeros(1,TryTime);
%         sin_sum=0;
%         cos_sum=0;        
%         tempX=repmat([0,0;0,0;0,1],[1,1,numBall]);
%         obj_function=0;
%         for Balln=1:numBall   %calculate the objctive funtion for every ball and store them in obj_f
%             tempX(:,2,Balln)=X(:,:,Balln)*rads(Balln);
%             r_locus=rads(Balln)*rc/sqrt(rads(Balln).^2+rc.^2);
%             tempX(:,:,Balln)=ArbAxisRotate(RotPoint,[0,0,rads(Balln)],-theta*rc/r_locus,tempX(:,:,Balln));
%             tempX(:,:,Balln)=RotateZ(theta)*tempX(:,:,Balln);
%             zaxis=(tempX(:,:,Balln)-tempX(:,:,Balln))/rads(Balln);
%             sin_sum=sin_sum+zaxis(2);
%             cos_sum=cos_sum+zaxis(1);
%         end
%         Theta_Av=atan2(sin_sum,cos_sum);
%         for Balln=1:numBall
%             zaxis=(tempX(:,:,Balln)-tempX(:,:,Balln))/rads(Balln);
%             theta_ori=(atan2(zaxis(2),zaxis(1))-Theta_Av).^2;
%             obj_function=obj_function+theta_ori;
%         end
%         for numStp=1:TryTime
%             tryAlpha=numStp*20*pi/TryTime;
%             trial_result(numStp)=subs(obj_function,theta,tryAlpha);
%         end
%         [~,IN]=min(trial_result);
%         alpha=IN*2*pi/TryTime;
%     end

    function Y=ArbAxisRotate(initP,endP,RotAng,X)
        %this function using quaternions to make rotation about arbitrary axis in
        %the 3D space.
        %INPUT:
        %initP: initial point of the rotation axis
        %endP: end point of the rotation axis
        %RotAng:rotation angle
        %X:matirx to be rotated
        %Yu Huang 2015, Email: michael.williams.hy@gmail.com
        axis_module=sqrt((endP(1)-initP(1)).^2+(endP(2)-initP(2)).^2+(endP(3)-initP(3)).^2); %normalize
        % q=[cos(RotAng./2) RotAxis(1).*sin(RotAng./2)/axis_module RotAxis(2).*sin(RotAng./2)/axis_module RotAxis(3).*sin(RotAng./2)/axis_module];
        % %making the quaternions for rotation
        % Y=quatrotate(q,X);
        a=initP(1);
        b=initP(2);
        c=initP(3);
        u=(endP(1)-initP(1))/axis_module;
        v=(endP(2)-initP(2))/axis_module;
        w=(endP(3)-initP(3))/axis_module;
        % x=X(1);
        % y=X(2);
        % z=X(3);
        % Y=[(a*(v.^2+w.^2)-u*(b*v+c*w-u*x-v*y-w*z))*(1-cos(RotAng))+x*cos(RotAng)+(-c*v+b*w-w*y+v*z)*sin(RotAng);
        %    (b*(u.^2+w.^2)-v*(a*u+c*w-u*x-v*y-w*z))*(1-cos(RotAng))+y*cos(RotAng)+(c*u-a*w+w*x-u*z)*sin(RotAng);
        %    (c*(u.^2+v.^2)-w*(a*u+b*v-u*x-v*y-w*z))*(1-cos(RotAng))+z*cos(RotAng)+(-b*u+a*v-v*x+u*y)*sin(RotAng)];
        RotMat=[u.^2+(v.^2+w.^2)*cos(RotAng),u*v*(1-cos(RotAng))-w*sin(RotAng),u*w*(1-cos(RotAng))+v*sin(RotAng),(a*(v.^2+w.^2)-u*(b*v+c*w))*(1-cos(RotAng))+(b*w-c*v)*sin(RotAng);
            u*v*(1-cos(RotAng))+w*sin(RotAng),v.^2+(u.^2+w.^2)*cos(RotAng),v*w*(1-cos(RotAng))-u*sin(RotAng),(b*(u.^2+w.^2)-v*(a*u+c*w))*(1-cos(RotAng))+(c*u-a*w)*sin(RotAng);
            u*w*(1-cos(RotAng))-v*sin(RotAng),v*w*(1-cos(RotAng))+u*sin(RotAng),w.^2+(u.^2+v.^2)*cos(RotAng),(c*(u.^2+v.^2)-w*(a*u+b*v))*(1-cos(RotAng))+(a*v-b*u)*sin(RotAng);
            0,0,0,1];
        temp1=[X(1,1),X(1,2);X(2,1),X(2,2);X(3,1),X(3,2);1,1];
        temp=RotMat*temp1;
        Y=[temp(1,1),temp(1,2);temp(2,1),temp(2,2);temp(3,1),temp(3,2)];
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
