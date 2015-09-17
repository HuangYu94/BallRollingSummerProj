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
rads = [10,13,17,19]; %balls` radii for experiment
numBall=numel(rads);
configX=repmat(eye(4),[1,1,numBall]); %configuration for each ball.
radmin=min(rads);
if nargin<2
    thetaL = 2*pi; R = radmin;
end
rc=zeros(1,numBall); %make the rc for lollipop
for n=1:numBall  
    rc(n)=rads(n).^2/sqrt(R.^2+rads(n).^2);
end
trial_times=1000;
stp=zeros(1,trial_times+1);
error_rec=zeros(numBall,trial_times+1);%record the error repectively
error_total=zeros(1,trial_times+1);
% increamental=thetaL/trial_times;
for i=0:trial_times
    trial_ang=i*thetaL./trial_times;
    stp(i+1)=trial_ang*180/pi;
    for n=1:numBall
%         qq=angle2quat(trial_ang,0,0);
%         RotVec=quatrotate(qq,[R 0 rads(n)]);
        temp=ArbQuatRotate([R,0,rads(n)],-trial_ang*R/rads(n),configX(:,:,n));
        configX(:,:,n)=ArbQuatRotate([0,0,1],trial_ang,temp);%make the lollipop rotation
        zaxis=[0,0,0,1]*configX(:,:,n);
        error_rec(n,i+1)=acos(zaxis(4));
    end
    error_total(i+1)=sum(abs(error_rec(:,i+1)));
end
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

function Y=ArbQuatRotate(RotAxis,RotAng,X)
%this function using quaternions to make rotation about arbitrary axis in
%the 3D space.
%INPUT:
%axis: central rotation axis which are not required to be a unite vector
%RotAng:rotation angle
%X:matirx to be rotated
%Yu Huang 2015, Email: michael.williams.hy@gmail.com
axis_module=sqrt(RotAxis(1)+RotAxis(2)+RotAxis(3)); %normalize
q=[cos(RotAng./2) RotAxis(1).*sin(RotAng./2)/axis_module RotAxis(2).*sin(RotAng./2)/axis_module RotAxis(3).*sin(RotAng./2)/axis_module];
q_inv=quatinv(q);
%making the quaternions for rotation
temp=quatmultiply(q,X);
Y=quatmultiply(temp,q_inv);
end