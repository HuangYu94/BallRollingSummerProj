function FindTheOptimalDirction
% This simulation is to find the optimal dirction 
% in order to decrease the error
% Michael Williams 2015, email:michael.williams.hy@gmail.com
format compact
clc
rads =[10,15,20];
numBall=numel(rads);
vector=repmat([0;0;0],[1,1,360]);
error_rec=ones(1,360);
temp=repmat(eye(3),[1,1,numBall]);
z_ori_init=ones(1,numBall);
stp=1:1:360;
iC=1:1:numBall;
TurnZ = repmat(eye(3),[1,1,numBall]);
for n = 1:numBall
    theta = pi/12/rads(n);
    TurnZ(:,:,n) =RotateZ(theta);
end
X = repmat(eye(3),[1,1,numBall]);
psi = zeros(numBall,1);
theta = zeros(numBall,1);
%initial rotation of each sphere
for n=1:numBall
    X(:,:,n)=RotateY(pi*2/rads(n))*RotateX(pi/rads(n))*X(:,:,n);
end

    %display(k)
    %rotate all balls around world z-axis, then measure error psi, and 
    for n=1:numBall
        zTurn=pi/6/rads(n);
        temp(:,:,n)=RotateZ(zTurn)*X(:,:,n);
        zaxis=temp(:,:,n)*[0;0;1];
        psi(n) = acos(zaxis(3));
        z_ori_init(n)=abs(psi(n))*360/2/pi; %get the initial orientation
        if(abs(zaxis(3)>1))
            display('pji(n) too big!');
        end
        theta(n) = atan2(zaxis(2),zaxis(1));
    end
    d = rads(iC)*psi(iC);
    for i=1:360    %try every dirction by increase 1 degree
        for n=1:numBall
            alfa=2*pi*i/360;
            vector(:,:,i)=[cos(alfa);sin(alfa);0];
            temp(:,:,n)=rodRotM(vector(:,:,i), -d/rads(n)/8)*X(:,:,n);
            zaxis = temp(:,:,n)*[0;0;1];
            psi(n) = acos(zaxis(3));
            theta(n) = atan2(zaxis(2),zaxis(1));
        end
        error_rec(i)=sum(psi)*360/2/pi-sum(z_ori_init); %see which direction will make the error decrease
    end                                                 %I mean the overall error
    figure(1)
    plot(stp,error_rec);
    title('Which direction will decrease the overall error?')
    xlabel('direction(deg)')
    ylabel('overall error(deg)')
    %display(psi)
    % calculate ....
%     pow =2;
%     ux = -rads.^(-1)*(psi.^pow.*cos(theta))/numBall;
%     uy = -rads.^(-1)*(psi.^pow.*sin(theta))/numBall;
%     d=sqrt(ux^2+uy^2);
%     if(d>0)
%         vector=[uy;-ux;0]/d;
%     else
%         vector=[1,0,0]; %why?
%         display('we are going to make a big mve.  WHy?')
%     end
%     
    %apply the control law
%     for n = 1:numBall
        % the rotation is scaled by the radius

%     end

end
function vrot = rodRotM(k, theta)
%If v is a vector in \mathbb{R}^3 and k is a unit vector describing an axis of rotation
% about which we want to rotate v by an angle θ (in a right-handed sense), the Rodrigues
%formula is:
vrot = [cos(theta) + (1-cos(theta))*k(1)^2,     (1-cos(theta))*k(1)*k(2),               sin(theta)*k(2) ;
    (1-cos(theta))*k(1)*k(2),               cos(theta) + (1-cos(theta))*k(2)^2,   -sin(theta)*k(1) ;
    -sin(theta)*k(2),                       sin(theta)*k(1),                     cos(theta) ];


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