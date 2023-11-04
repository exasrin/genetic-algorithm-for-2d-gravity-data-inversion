%-----------------------------------------------------------------------
% Program Forward Model Gravitasi 2D
% Refererensi Last Kubik
% Dibuat Oleh Dr. Eng. Jamhir Safani
% Dikembangkan oleh Al Rubaiyn
% Kendari 28 Mei 2015
%-----------------------------------------------------------------------
clc; clear all;
G=6.6732e-11;   % Kontanta Gravirtasi damalm Mks

l0 =0.01; 
dx = 10; %(m)
dh =10; %(m)
rho = 1000; %(kg/m^3)
nx = 9; nz =4 ;

V = zeros(nz,nx);
% V(2,3)=rho;
% V(3,3)=rho;
% V(3,4)=rho;
% V(4,4)=rho;
% V(4,5)=rho;
% V(5,5)=rho;
V(2,4:6)=rho;
% V(2,5)=rho;
% V(2,6)=rho;
V(3,4:6)=rho;
% V(3,5)=rho;
% V(3,6)=rho;


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
figure 
subplot(2,1,1); plot(x,GG*1e5,'ro-','LineWidth',2.0,'MarkerSize',3);
title('PENAMPANG SINTETIK 2D GRAVITASI', 'fontweight','bold','fontsize',11)
ylabel('Anomali Medan Gravitasi [mGal]','fontsize',7);
xlabel('Spasi [m]','fontsize',7);
set(gca,'fontsize',7,'XMinorTick','off');
hleg= legend('g-obs','Location','NorthEast');
set(hleg,'FontAngle','italic','FontSize',7);

subplot(2,1,2); imagesc(x,z,V);
% subplot(2,1,2); contourf(x,-z,V);
ylabel('Kedalaman [m]','fontsize',7);
% ax=gca;
set(gca,'XAxisLocation','top','fontsize',7,'XMinorTick','on','YMinorTick','on');

colorbar('horiz');
colormap('default');
% colorbar('vert'); colormap(cool);
% view(0,270)
% grid on
% set(gca,'YLim',[-nz*dh 0])
% set(gca,'YTick',[-nz*dh:dh:0])
set(gca, 'GridLineStyle', '-');
% set(gca,'XGrid','on','YGrid','on','ZGrid','on')
% grid(gca,'minor');


hasil= [x' GG*1e5];%satuan dalam mGal
% hasil= [x' GG*1e8]; %satuan dalam mikroGal
dlmwrite('grav_obs.dat',hasil,'newline','pc','delimiter','\t','precision','%10.5f');