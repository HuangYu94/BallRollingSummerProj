function SimMovieBWspheresTest
% makes a movie simulating prototype 3 rolling 8 spheres from "U" to "I"
% to "U" to "C"
%RandStream.setDefaultStream(RandStream('mt19937ar','seed',1));
%TODO: fix perspective.
global FrameCount MOVIE_NAME MAKE_MOVIE
clc
format compact

MAKE_MOVIE = 0; 

randPert  = 0.05;
iC =[1,2,4,5,6];
if randPert == 0
    iC = [1,2,4,5,6];
end

scaler = 8;
    
%mogrify - crop <width>x<height><offset x><offset y>
%mogrify -crop 901x726+490+93 *.png
%cd ~/Documents/motionlab/trunk/documents/Aaron/MagneticSphere/Matlab
%ffmpeg -r 50 -f image2 -i BlackWhite12spheres%06d.png -b 8000k -r 30 BlackWhite12spheres.mp4
%ffmpeg -r 30 -f image2 -i BlackWhite12spheres%06d.png -b 8000k -r 60 BlackWhite12spheres.mp4

MOVIE_NAME = 'BlackWhite12spheres';
numStep =800;
FrameCount = 0;


BallSet = [25.6,24.3, 22.43, 19.21,18.23];
BallAve = mean(BallSet);
rs = BallAve./BallSet;


ballTypes = [     1, 2, 1;   
                 3, 4, 5;
                 3, 4, 5;
                 1, 1, 1];
ZUP = [0;0;1];  % white
ZDN = [0;0;-1]; % black
%         1    2    3    4   5
%GoalS = [ZUP, ZUP, ZUP, ZUP,ZUP];
GoalU = [ZDN, ZUP, ZDN, ZUP,ZDN];
GoalI = [ZDN, ZDN, ZUP, ZDN,ZUP];
GoalC = [ZDN, ZDN, ZDN, ZUP,ZUP];


             
numBall =numel(ballTypes); % number of spheres

%%%% set up plots
figure(1);
close(1);
figure(1);
set(gcf,'Color',[202 	225 	255]/255);   %cornflower blue background
clf;
set(gcf, 'position',[6          57        1670         915]);  % full screen
ax = axes;%('XLim',[-3,3],'YLim',[0,18.5],'ZLim',[-2,3]);

axis equal
hold on
set(gcf,'Renderer','opengl')

hold on
axis off
%make spheres
[X,Y,Z]=sphere(40);
s = zeros(numBall);

% Make 12 spheres
Rsphere = cell(numBall,1);
Tsphere = cell(numBall,1);
rads = ones(numBall,1);

sc = 1;
spacing =2.4/sc;


c = 1;
for i = 1:size(ballTypes,1)
    for j = 1:size(ballTypes,2)
        %                    [x,y,w,h]
        rad = 1/rs(ballTypes(i,j))/sc;
        rads(c) = rad;
        Xs = rad*X;
        Ys = rad*Y;
        Zs = rad*Z;
            img=imread('BWsphere2.png'); %img=imread('BWsphere.png');
        spherepic = warp(Xs,Ys,Zs,img);
        s(c) = hgtransform('Parent',ax);
        set(spherepic,'Parent',s(c))
        Rsphere{c} = makehgtform('zrotate',pi/2,'yrotate',0);
 

        Tsphere{c} = makehgtform('translate',[spacing*(j-1),spacing*(i-1),0]);
        set(s(c),'Matrix',Tsphere{c}*Rsphere{c}) 
        c=c+1;
    end
end




%lightangle(-45,30)
 camlight(-59,154)
 camlight(30,180,'infinite')
 camlight(90,90)
%set(gcf,'Renderer','zbuffer')
set(findobj(gca,'type','surface'),...
    'FaceLighting','phong',...%'AmbientLightColor',[.5,0,1],...
    'AmbientStrength',1.3,'DiffuseStrength',1,...
    'SpecularStrength',1.0,'SpecularExponent',10,'SpecularColorReflectance',.8,...
    'BackFaceLighting','lit')
material shiny



view(0,90)
axis equal

axisDef = axis;
%view(-21,38);
%axis equal
%axisDef = [-3.5   7   -.5    6.76   -1.    1.];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SIMULATION %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% plot the balls moving.


TurnZ = repmat(eye(3),[1,1,numBall]);
for n = 1:numBall
        theta = pi/12/rads(n);
        TurnZ(:,:,n) =RotateZ(theta);  
end
X = repmat(eye(3),[1,1,numBall]);
psi = zeros(numBall,1);
theta = zeros(numBall,1);

GoalU = [ZDN, ZUP, ZDN,    ZDN, ZUP, ZDN,    ZDN, ZUP, ZDN,      ZDN, ZDN,ZDN];
GoalI = [ZDN, ZDN, ZDN,    ZUP, ZDN, ZUP,    ZUP, ZDN, ZUP,      ZDN, ZDN,ZDN];
GoalC = [ZDN, ZDN, ZDN,    ZDN, ZUP,ZUP,     ZDN, ZUP,ZUP,       ZDN, ZDN,ZDN];

goalTimes = [2,numStep/4,2*numStep/4,3*numStep/4,numStep];
SPEED = 0.05;%0.05; %0.01
for k = 1:numStep
    
    if sum( k ==  goalTimes)
        saveLotsFrames()
        %pause(2)
    end
    
    if(     k < goalTimes(2))
        Goal = GoalU;
    elseif (k < goalTimes(3))
        Goal = GoalI;
    elseif (k < goalTimes(4))
        Goal = GoalU;
    elseif (k < goalTimes(5))
        Goal = GoalC;
    end
    
    % 'turn' in place by moving in a square designed so that it moves the
    % nominal sphere approx pi/4 about the z-axis
    for cnt = 1:1
        for n = 1:numBall
            rM = RandRot(randPert);
            zTurn = pi/6/rads(n);
            
            hgtRot = diag([1,1,1,1]);
            hgtRot(1:3,1:3) = rM*RotateZ(zTurn);
            
            X(:,:,n) =   rM*RotateZ(zTurn)*X(:,:,n);
            %Rsphere{n} = makehgtform('zrotate',zTurn)*Rsphere{n};
            Rsphere{n} = hgtRot*Rsphere{n};
            
            Tx = makehgtform('translate',[0,0,0])*Tsphere{n};  %??? why not Tx=Tsphere{n}?
            set(s(n),'Matrix',Tx*Rsphere{n})

            % Measure the resulting error
            zaxis = Goal(3,n)*X(:,:,n)*[0;0;1];     %??? what do you mean by "measure the resulting error"?
            psi(n) = acos(zaxis(3));

            if(abs(zaxis(3))>1)
                display('pji(n) too big!');
            end
            theta(n) = atan2(zaxis(2),zaxis(1));
        end
        axis(axisDef);
        if cnt  == 1
            updateDrawing
        end
    end
    
    %find the error (large balls have larger errors, then we average the
    %results*)
    ux = -(rads(iC).^(-1))'*(psi(iC).^2.*cos(theta(iC)))/numel(iC);
    uy = -(rads(iC).^(-1))'*(psi(iC).^2.*sin(theta(iC)))/numel(iC); %.^(-1)
    
    % move in a straight line to minimize this error  
    d = sqrt(ux^2+uy^2);
    if(d > 0)
        vector = [uy;-ux;0]/d;
    else
        vector = [1,0,0];
    end
   
    for n = 1:numBall
        % the rotation is scaled by the radius
        X(:,:,n) =  rodRotM(vector, -d/rads(n)/scaler)*X(:,:,n); %???what`s the function of the scaler?
        %makehgtform('axisrotate',[ax,ay,az],t)
        Rsphere{n} = makehgtform('axisrotate',vector,-d/rads(n)/scaler)*Rsphere{n};
        Tx = makehgtform('translate',[0,0,0])*Tsphere{n};
        set(s(n),'Matrix',Tx*Rsphere{n})
        
        zaxis = Goal(3,n)*X(:,:,n)*[0;0;1];
        psi(n) = acos(zaxis(3));
        theta(n) = atan2(zaxis(2),zaxis(1));
    end
    axis(axisDef);
    updateDrawing
end


saveLotsFrames()
set(gcf,'name','done')


function saveLotsFrames()
    for i = 1:60  % save final position
        updateDrawing
    end

function updateDrawing
global FrameCount MOVIE_NAME MAKE_MOVIE
drawnow
if(MAKE_MOVIE)
     
            FrameCount=FrameCount+1;
            if FrameCount >= 1237
                F = getframe(gcf);
                fname = sprintf('%s%06d.png',MOVIE_NAME,FrameCount);
                imwrite(F.cdata, fname,'png'); 
            end
            while(FrameCount < 10)
                updateDrawing
            end
            

end


function vrot = rodRotM(k, theta)
 %If v is a vector in \mathbb{R}^3 and k is a unit vector describing an axis of rotation 
 % about which we want to rotate v by an angle θ (in a right-handed sense), the Rodrigues
 %formula is:
    vrot = [cos(theta) + (1-cos(theta))*k(1)^2,     (1-cos(theta))*k(1)*k(2),               sin(theta)*k(2) ;
            (1-cos(theta))*k(1)*k(2),               cos(theta) + (1-cos(theta))*k(2)^2,   -sin(theta)*k(1) ;
            -sin(theta)*k(2),                       sin(theta)*k(1),                     cos(theta) ];


% function RxTh = RotateX(theta) 
%     RxTh = [1,  0,  0;
%             0, cos(theta), -sin(theta);
%             0, sin(theta),  cos(theta)];
%  end   
% function RyTh = RotateY(theta)
%     RyTh = [ cos(theta), 0, sin(theta);
%              0,  1,  0;
%             -sin(theta), 0, cos(theta)];
% end      
function RzTh = RotateZ(theta) 
    RzTh = [cos(theta),  -sin(theta),0;
            sin(theta),   cos(theta),0;
            0,  0,  1];
function M = RandRot(d)
% Uniformly distributed random 3D rotation matrix using Arvo's method.
% "Fast Random Rotation Matrices", by James Arvo. In Graphics Gems III, 1992. http://www.ics.uci.edu/~arvo/papers/RotationMat.ps http://www.ics.uci.edu/~arvo/code/rotate.c (2011-06-22)

% To generate uniformly distributed random rotations of a unit sphere, first perform
% a random rotation about the vertical axis, then rotate the north pole to a random
% position.

% /*======================================================================*
%  *  R A N D _ R O T A T I O N      Author: Jim Arvo, 1991               *
%  *                                                                      *
%  *  This routine maps three values (x[0], x[1], x[2]) in the range      *
%  *  [0,1] into a 3x3 rotation matrix M.  Uniformly distributed random   *
%  *  variables x0, x1, and x2 create uniformly distributed random        *
%  *  rotation matrices.  To create small uniformly distributed           *
%  *  "perturbations", supply samples in the following ranges             *
%  *                                                                      *
%  *      x[0] in [ 0, d ]                                                *
%  *      x[1] in [ 0, 1 ]                                                *
%  *      x[2] in [ 0, d ]                                                *
%  *                                                                      *
%  * where 0 < d < 1 controls the size of the perturbation.  Any of the   *
%  * random variables may be stratified (or "jittered") for a slightly    *
%  * more even distribution.                                              *
%  *                                                                      *
%  *======================================================================*/
%void rand_rotation( float x[], Matrix3 *M )
    x = rand(3,1);

    theta = x(1)*d * 2 *pi; %/* Rotation about the pole (Z).      */
    phi   = x(2)   * 2 *pi; %/* For direction of pole deflection. */
    z     = x(3)*d * 2.0;      %/* For magnitude of pole deflection. */

    %/* Compute a vector V used for distributing points over the sphere  */
    %/* via the reflection I - V Transpose(V).  This formulation of V    */
    %/* will guarantee that if x[1] and x[2] are uniformly distributed,  */
    %/* the reflected points will be uniform on the sphere.  Note that V */
    %/* has length sqrt(2) to eliminate the 2 in the Householder matrix. */

    r  = sqrt( z );
    Vx = sin( phi ) * r;
    Vy = cos( phi ) * r;
    Vz = sqrt( 2.0 - z );    

    %/* Compute the row vector S = Transpose(V) * R, where R is a simple */
    %/* rotation by theta about the z-axis.  No need to compute Sz since */
    %/* it's just Vz.                                                    */

    st = sin( theta );
    ct = cos( theta );
    Sx = Vx * ct - Vy * st;
    Sy = Vx * st + Vy * ct;

    %/* Construct the rotation matrix  ( V Transpose(V) - I ) R, which   */
   % /* is equivalent to V S - R.                                        */

   M = [-Vx * Sx + ct,  -Vx * Sy + st, -Vx * Vz;
         -Vy * Sx - st,  -Vy * Sy + ct, -Vy * Vz;
        Vz * Sx,        Vz * Sy,      1.0 - z];   %ATB: fixed this
    
    M = RotateZ(rand(1)*d) ;






