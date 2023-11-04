function [gm]=forward_gravity(rho,model)
G=6.6732e-11;   % Kontanta Gravirtasi damlam Mks

l0 =0.01; 
nx=model(1);    %Banyak cell lateral
nz=model(2);    %Banyak cell vetikal
dx =model(3); %Dimensi cell lateral (m)
dh =model(4); %Dimensi cell vetikal (m) 

V =reshape(rho,nx,nz);  V=V';

VV = V(:);

d=dx; %station spacing
h=dh; %depth spacing (thickness)
 
x1 = [0:dx:(nx-1)*dx];
x = x1+d/2;
 
z1 = [0:h:(nz-1)*h]';
z = z1+h/2;
zz = repmat(z,1,nx);
xx = repmat(x,nz,1);
 
nb = nx*nz;
 
for i=1:nx
    for j = 1:nb
        r1 = ((zz(j)-h/2).^2 + (x(i)-xx(j)+d/2).^2).^0.5;
        r2 = ((zz(j)+h/2).^2 + (x(i)-xx(j)+d/2).^2).^0.5;
        r3 = ((zz(j)-h/2).^2 + (x(i)-xx(j)-d/2).^2).^0.5;
        r4 = ((zz(j)+h/2).^2 + (x(i)-xx(j)-d/2).^2).^0.5;
        theta1 = atan((x(i)-xx(j)+d/2)/(zz(j)-h/2));
        theta2 = atan((x(i)-xx(j)+d/2)/(zz(j)+h/2));
        theta3 = atan((x(i)-xx(j)-d/2)/(zz(j)-h/2));
        theta4 = atan((x(i)-xx(j)-d/2)/(zz(j)+h/2));
 
        %theta1 = theta1*180/pi; theta2 = theta2*180/pi; theta3 = theta3*180/pi; theta4 = theta4*180/pi;
        
        AA(i,j) = 2*G*(...
            (x(i)-xx(j)+d/2)*log((r2*r3)/(r1*r4))...
            + d*log(r4/r3)...
            - (zz(j)+h/2)*(theta4-theta2)...
            + (zz(j)-h/2)*(theta3-theta1));
        
    end
     
end

GG = AA*VV;
gm = GG*1e5;
