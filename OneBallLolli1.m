function OneBallLolli1(ang_x,ang_y)
%this function simulate the lollipop rolling of one ball whose Z-axis is
%randomly oriented. Our purpose of rolling the ball like this is to make
%the z-axis of the ball aligned to the world coordinates` z-axis.
%INPUT:
%ang_x: initial rotation around x-axis
%ang_y: initial rotation around y-axis
%================Optimization Varibles==================%
%!!!arc:the angle related to the arc along which the ball will roll
%================Optimization Varibles==================%
%alfac: we let the ball centered at the origin and pick up a point as
%the center of the circle the ball rolled along and this is the angle between 
%line segment of circle center and ball`s start point and the X axis of the
%world.
%R: radus of the circle.
%Yu Huang, 2015, E-mail: Michael.Williams.hy@gmail.com
format compact
if nargin<2
    ang_x = 0.7*pi;ang_y = 0.9*pi;
end
rad=1;
precision=0.01;
R=sqrt(5)/2;
rc=R*rad/sqrt(R.^2+rad.^2);
configBall=[0,0;0,0;0,rad]; %configuration of the ball`s z axis(2 point).
error_rec=zeros(1,100);
% arc_rec=zeros(1,1000);
%make the initial rotation:
configBall=RotateX(ang_x/rad)*RotateY(ang_y/rad)*configBall;
Zerror=1;
stp=0;
step1=zeros(1,100);
tra_arc=zeros(1,100);
while Zerror>precision
    zaxis=(configBall(:,2)-configBall(:,1))/rad;
    psi=acos(zaxis(3));
    Rc=rad*tan(psi/2);
    stp=stp+1;
    step1(stp)=stp;
    StepArc=FindArc(configBall,Rc,rc,rad,zaxis);
%     StepArc=pi*rc/R;
    tra_arc(stp)=StepArc;
    %make the lollipop rotation:
    configX=ArbAxisRotate([-Rc*zaxis(1),-Rc*zaxis(2),0],[0,0,rad],-StepArc*R/rc,configBall);
    configX=ArbAxisRotate([0,0,0],[0,0,rad],StepArc,configX);
    %calculate the zaxis of the ball after rotation
    zaxis=configX(:,2)-configX(:,1);
    %make the objective function:
    Zerror=acos(zaxis(3));
    error_rec(stp)=Zerror;
    figure(1);
plot(step1(1:stp),180*error_rec(1:stp)/pi);
title('Lollipop Rolling of one Ball to Align Z-axis of the Ball');
xlabel('steps');
ylabel('error in deg');
end
%we first try several values of arc and choose one that could make the GDM
%work better!
%plot the error curve:
figure(1);
plot(step1(1:stp),180*error_rec(1:stp)/pi);
title('Lollipop Rolling of one Ball to Align Z-axis of the Ball');
xlabel('steps');
ylabel('error in deg');
%plot the trajectory:
interP=200; %we use the interpolation method to plot the arcs.
initP=[0,0];
for j=1:stp
    Xvar=zeros(1,interP+1);
    Yvar=zeros(1,interP+1);
    distanceP=tra_arc(j)/interP;
    for k=0:interP
        CtrAng=k*distanceP;
        Xvar(k+1)=R*cos(CtrAng)-R+initP(1);
        Yvar(k+1)=R*sin(CtrAng)+initP(2);
    end
    figure(2);
    title('Trajectory of the Ball');
    axis('equal');
    xlabel('X (cm)');
    ylabel('Y (cm)');
    plot(Xvar,Yvar,'r',Xvar(end),Yvar(end),'b.',0,0,'go');
    legend('circular trajectory','turning point','start point');
    hold on
    initP(1)=Xvar(end);
    initP(2)=Yvar(end);
end
end

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
%     function RzTh = RotateZ(theta)
%         RzTh = [cos(theta),  -sin(theta),0;
%             sin(theta),   cos(theta),0;
%             0,  0,  1];
%     end

function stepArc=FindArc(configBall,Rc,rc,rad,zaxis)
trial_num=1000; %number of values for us to try.
trial_result=zeros(1,trial_num);
% syms arc
% arc=5;
% configXtemp=ArbAxisRotate([0,-R,0],[0,0,rad],-arc*R/rc,configBall);
% configXtemp=ArbAxisRotate([0,-R,0],[0,-R,1],arc,configXtemp);
% zaxistemp=configXtemp(:,2)-configXtemp(:,1);
% obj_fun=acos(zaxistemp(3));
for i=1:trial_num
    try_value=2*pi*i/trial_num;
    configXtemp=ArbAxisRotate([-Rc*zaxis(1),-Rc*zaxis(2),0],[0,0,rad],-try_value*Rc/rc,configBall);
    configXtemp=ArbAxisRotate([0,0,0],[0,0,rad],try_value,configXtemp);
    %calculate the zaxis of the ball after rotation
    zaxistemp=configXtemp(:,2)-configXtemp(:,1);
    trial_result(i)=acos(zaxistemp(3));
end
[~,in]=min(trial_result);
startP=2*pi*in/trial_num;
%do the gradient descent method:
% X_old=100;
stepArc=startP;

% while abs(X_old-X_new)>0.0000001
%     X_old=X_new;
%     gradF=vpa(subs(obj_fun,arc,X_new));
%     gamma=FindGamma(gradF,X_new);
%     X_new=X_old-gamma*gradF;
% end
% stepArc=X_new;
% % display(X_new);
% % zaxis0=configBall(:,2)-configBall(:,1);
% % error0=acos(zaxis0(3));
% % error1=subs(obj_fun,arc,X_new);
% % display('error0:');
% % display(error0*180/pi)
% % display('error1:');
% % display(error1*180/pi);
%     function gamma=FindGamma(gradF,X_new)
%         %I want to implement the 0.618 method, because I cannot add the function
%         %handle in matlab
%         a1=0;
%         b1=2;
%         GR=(sqrt(5)-1)/2; %global variable for Golden Section Search (https://en.wikipedia.org/wiki/Golden_section_search)
%         d1=GR*(b1-a1)+a1;
%         c1=b1-GR*(b1-a1);
%         while abs(c1-d1)>0.000000000001
%             ffc=subs(obj_fun,arc,X_new-gradF*c1);
%             ffd=subs(obj_fun,arc,X_new-gradF*d1);
%             if ffc<=ffd
%                 b1=d1;
%                 d1=c1;
%                 c1=b1-GR*(b1-a1);
%             else
%                 a1=c1;
%                 c1=d1;
%                 d1=a1+GR*(b1-a1);
%             end
%         end
%         gamma=(a1+b1)/2;
%     end
end