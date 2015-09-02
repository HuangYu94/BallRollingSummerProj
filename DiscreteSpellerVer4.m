% Drive 5 spheres in  a discrete time fashion from starting to final
% orientation.

% Version 3 rolls the ball back and forth along the y axis, and controls
% the x axis movement to make the balls align


% DiscreteSpeller.m works, in that it steers the balls to positions to form
% U, I, U, C.  However, it relies on perfect turning in place around the
% z-axis and sometimes increases the error (should asymptotically decrease)

% DiscreteSpellerVer4  follows the gradient.

function DiscreteSpellerVer4
format compact
clc
close all

XUP = [1;0;0];  % white
XDN = [-1;0;0]; % black



%Initialize 5 spheres and give them goals
BallSet = [25.6,24.3, 22.43, 19.21,18.23];
numBall = numel(BallSet);
BallAve = mean(BallSet);
rs = BallAve./BallSet;

numStep =  20000;

GoalU = [XDN, XUP, XDN, XUP,XDN];
GoalI = [XDN, XDN, XUP, XDN,XUP];
GoalC = [XDN, XDN, XDN, XUP,XUP];

% starting orientation is all up
X = repmat(eye(3),[1,1,numBall]);


phi = zeros(numBall,1);
theta = zeros(numBall,1);

PHI = zeros(numBall,2*numStep);
THETA = zeros(numBall,2*numStep);
PATH = zeros(2,2*numStep+1);
numneg = 0;
numpos = 0;
numzero = 0;

scaler = 32;
for k = 1:numStep
    

    
    if(     k < 5000)
        Goal = GoalU;
    elseif (k < 10000)
        Goal = GoalI;
    elseif (k < 15000)
        Goal = GoalU;
    elseif (k < 20000)
        Goal = GoalC;
    end
    grad = zeros(numBall,1);
    
    
    for n = 1:numBall
    	% 'turn' about the x-axis (doesn't increase the error)   
        X(:,:,n) =  RotateX(pi/(16*rs(n)))*X(:,:,n);
        
        % Measure the resulting error
        xaxis = Goal(1,n)*X(:,:,n)*[1;0;0];
        phi(n) = acos(xaxis(1));  %distance error (latitude)
        
        theta(n) = atan2(xaxis(2),xaxis(3)); %longitude
        
        if(phi(n) > 0 && cos(phi(n)/rs(n))^2 < 1)
            grad(n) = (cos(theta(n)/rs(n))*sin(phi(n)/rs(n)))/(rs(n)*(1 - cos(phi(n)/rs(n))^2)^.5);
        end
       
    end
    %Record the phi and theta
    PATH(:,2*k) = PATH(:,2*k-1)+[0;pi/16];
    PHI(:,2*k-1) = phi;
    THETA(:,2*k-1) = theta;
    
  
    
    %find the error (large balls have larger errors, then we average the
    %results*)
    %ux = sum(-grad)/10;
    ux = -(rs.^(0))*(phi.*cos(theta))/numBall;
    uy = 0;%-(rs.^(-1))*(phi.*sin(theta))/numBall; %.^(-1)
    
    %If v is a vector in \mathbb{R}^3 and k is a unit vector describing an axis of rotation about which we want to rotate v by an angle θ (in a right-handed sense), the Rodrigues formula is:
    %   v_{rot} = v \cos\theta + (\mathbf{k} \times \mathbf{v})\sin\theta + \mathbf{k} (\mathbf{k} \cdot \mathbf{v}) (1 - \cos\theta). 
    
    %alpha = atan2(uy,ux);
    d = sqrt(ux^2+uy^2);
    if(d < 1.1525e-15)
        ux = 1/scaler;
        d = sqrt(ux^2+uy^2);
    end
    vector = [uy;-ux;0]/d;
  
     for n = 1:numel(BallSet)
        % the rotation is scaled by the radius
        %X(:,:,n) =  rodRotM(vector, -d/rs(n)/scaler)*X(:,:,n); 
        %X(:,:,n) =  RotateY(-ux/rs(n)/scaler)*X(:,:,n);
        X(:,:,n) =  RotateY(-ux/rs(n)/scaler)*X(:,:,n); 
        % Measure the resulting error
        xaxis = Goal(1,n)*X(:,:,n)*[1;0;0];
        phi(n) = acos(xaxis(1));
        theta(n) = atan2(xaxis(2),xaxis(3));
    end

    %Record the phi and theta
     PATH(:,2*k+1) = PATH(:,2*k)+[ux;uy]/scaler;
    PHI(:,2*k) = phi;
    THETA(:,2*k) = theta;
    
end

display(['numpos=',num2str(numpos),', numneg=',num2str(numneg),', numzero=',num2str(numzero)])
% Plot the spheres position in phi and theta

colors = ['.-b';'.-g';'.-m';'.-c';'.-k';'.-y'];
figure(1)
clf
    for n = 1:numel(BallSet)
        plot( rs(n)*PHI(n,:),colors(n,:))
        hold on
    end
plot( sum(repmat(rs,size(PHI,2),1)'.*PHI,1),'.-r')
title('Phi error*radius')

figure(2)
clf
    for n = 1:numel(BallSet)
        plot( PHI(n,:),colors(n,:))
        hold on
    end
plot( sum(PHI,1),'.-r')
title('Phi error')
    
%     figure(2)
%     clf
%     for n = 1:numel(BallSet)
%         plot( THETA(n,:),colors(n,:))
%         hold on
%     end
%     title('theta')
    
    figure(3)
    clf
    plot(PATH(1,:),PATH(2,:),'.-r',PATH(1,1),PATH(2,1),'gx',PATH(1,end),PATH(2,end),'bo');
    axis equal
    legend('path','start','end')
    xlabel('units (radians)')
    ylabel('units (radians)')
    title('Path of ball center')
    
    
    figure(4)
    clf
    for n = 1:numel(BallSet)
        polar( THETA(n,:), rs(n)*PHI(n,:),colors(n,:))
        hold on
    end
    title('\theta vs \psi Polar Plot of ball positions')
    
%     figure(5)
%     clf
%     plot(THETA(1,:)-THETA(2,:),'.-')
%     title('Delta between thetas')
%     hold on
%    plot( signnum,'o-g')
% plot total error as a function of time

figure(1)
end

function vrot = rodRot(v, k, theta)
 %If v is a vector in \mathbb{R}^3 and k is a unit vector describing an axis of rotation 
 % about which we want to rotate v by an angle θ (in a right-handed sense), the Rodrigues
 %formula is:
    vrot = v *cos(theta) + cross(k,v)*sin(theta) + k* (dot(k,v)) *(1 - cos(theta)); 
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
        

