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