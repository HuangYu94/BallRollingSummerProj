function LollipopRoll(thetaL,R)
%this one demonstrate the so-called Lollipop Roll, which means the ball are
%rolling along a circular trajectory. This can be implemented on the CNC
%work panel easily. And by doing this, we can alter the orientation of the
%Z-axis of a ball.
%Yu Huang, 2015, Email:Michael.Willians.hy@gmail.com
%INPUT:
%thetaL: the corresponding angle of the arc along which the ball rotated.
%R:the circle trajectory for rotation.
format compact
clc
rads =1; %balls` radii for experiment
numBall=numel(rads);
configX=repmat([0,0;0,0;0,1;1,1],[1,1,numBall]); %configuration for each ball.
% radmin=min(rads);
if nargin<2
    thetaL = 2*pi; R =sqrt(5)/2 ;
end
rc=zeros(1,numBall); %make the rc for lollipop
for n=1:numBall  
    rc(n)=rads(n)*R/sqrt(R.^2+rads(n).^2);
    configX(:,1,n)=[R;0;rads(n);1]; %denote the center of sphere
    configX(:,2,n)=[R;0;2*rads(n);1]; %denote the orientation of the Z axis
end
trial_times=4000;
stp=zeros(1,trial_times+1);
error_rec=zeros(numBall,trial_times+1);%record the error repectively
error_total=zeros(1,trial_times+1);
% increamental=thetaL/trial_times;
for i=0:trial_times
    trial_ang=i*thetaL./trial_times;
    stp(i+1)=trial_ang*180/pi;
    sumz1=0;
    sumz2=0;
    for n=1:numBall
%         qq=angle2quat(trial_ang,0,0);
%         RotVec=quatrotate(qq,[R 0 rads(n)]);
%make the lollipop rotation
%first rotate the end point of the Z axis:
        temp=ArbQuatRotate([0;0;0],[R;0;rads(n)],-trial_ang*R/rc(n),configX(:,:,n));
        configX(:,:,n)=ArbQuatRotate([0;0;0],[0;0;1],trial_ang,temp);%first rotate the end point of the Z axis
        %then rotate the center of the sphere:
%         configX(:,1,n)=ArbQuatRotate([0;0;0],[0;0;1],trial_ang,configX(:,1,n));
        zaxis=configX(:,2,n)-configX(:,1,n);
        zaxisModule=sqrt(sum(zaxis.^2));
        error_rec(n,i+1)=acos(zaxis(3)/zaxisModule);
%         error_rec(n,i+1)=atan2(zaxis(2),zaxis(1));
%         sumz1=sumz1+zaxis(1);
%         sumz2=sumz2+zaxis(2);
    end
    error_total(i+1)=atan2(sumz2,sumz1);
end
save('LollipopData','error_rec');
%then we make the plot
color_arr=['y','m','c','r','g','b','w','k'];
figure(1);
plot(stp, 180/pi*error_total,'k' );
hold on
legendinfo=cell(1,numBall);
legendinfo{1}='error sum';
for n=1:numBall
    plot(stp,180/pi*error_rec(n,:),color_arr(n));
    legendinfo{n+1}=['rad:',num2str(rads(n))];
end
legend(legendinfo);
legend('boxoff');
title('Noiseless Ensemble Lollipop Rolling Control of 4 Spheres Orientation');
xlabel('rotation angle with respect to the circle center');
ylabel('abs(psi) Degs');
end

function Y=ArbQuatRotate(initP,endP,RotAng,X)
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
u=(endP(1)-initP(1))./axis_module;
v=(endP(2)-initP(2))./axis_module;
w=(endP(3)-initP(3))./axis_module;
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
Y=RotMat*X;
end