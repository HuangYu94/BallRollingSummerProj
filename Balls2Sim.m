function Balls2Sim
% Michael Williams 2015 1st simulation, email:michael.williams.hy@gmail.com
% about 2 balls convergence without pertubution
format compact
clc
tic
numStep =8000;
rads = [10,15,25,50];
numBall=numel(rads);
stp=zeros(numStep,1);
error_rec=zeros(numStep,numBall);


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
    X(:,:,n)=RotateX(pi)*X(:,:,n);
end

for k=1:numStep
    %display(k)
    %rotate all balls around world z-axis, then measure error psi, and 
    for n=1:numBall
        zTurn=pi/6/rads(n);
        X(:,:,n)=RotateZ(zTurn)*X(:,:,n);
        zaxis=X(:,:,n)*[0;0;1];
        psi(n) = acos(zaxis(3));
        if(abs(zaxis(3)>1))
            display('pji(n) too big!');
        end
        theta(n) = atan2(zaxis(2),zaxis(1));
    end
    %display(psi)
    % calculate ....
    pow =2;
    ux = -rads.^(-1)*(psi.^pow.*cos(theta))/numBall;
    uy = -rads.^(-1)*(psi.^pow.*sin(theta))/numBall;
    d=sqrt(ux^2+uy^2);
    if(d>0)
        vector=[uy;-ux;0]/d;
    else
        vector=[1,0,0]; %why?
        display('we are going to make a big mve.  WHy?')
    end
    
    %apply the control law
    for n = 1:numBall
        % the rotation is scaled by the radius
        X(:,:,n)=rodRotM(vector, -d/rads(n)/8)*X(:,:,n);
        zaxis = X(:,:,n)*[0;0;1];
        psi(n) = acos(zaxis(3));
        theta(n) = atan2(zaxis(2),zaxis(1)); 
        error_rec(k,n)=abs(psi(n));
        stp(k)=k;
    end
end
color_arr=['y','m','c','r','g','b','w','k'];
figure(1);
plot(stp, 180/pi*sum(error_rec,2),'k' );
hold on
legendinfo=cell(1,numBall);
legendinfo{1}='error sum';
for i=1:numBall
    plot(stp,180/pi*error_rec(:,i),color_arr(i));
    legendinfo{i+1}=['rad:',num2str(rads(i))];
end
legend(legendinfo);
legend('boxoff');
title('Noiseless Ensemble Control of two Spheres Orientation');
xlabel('steps');
ylabel('abs(psi) Degs');
toc
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
%  function RyTh = RotateY(theta)
%      RyTh = [ cos(theta), 0, sin(theta);
%               0,  1,  0;
%              -sin(theta), 0, cos(theta)];
%  end
function RzTh = RotateZ(theta)
RzTh = [cos(theta),  -sin(theta),0;
    sin(theta),   cos(theta),0;
    0,  0,  1];
end
end