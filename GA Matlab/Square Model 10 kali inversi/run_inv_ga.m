%------------------------------------------------------
% Inversion of Graity Data using GA
% Created by Asrin dan Al rubaiyn
% GA alfgorithm adopted from Yang, W.Y., dkk., (2020)
% Forward adopted from Last Kubik
%------------------------------------------------------

clc; clear; close all;
% Input Parameter Model
[infileextn,inpath] = uigetfile('*.xlsx','Select File InputInversion.xls...');
itin=xlsread(infileextn);
[ni mi]=size(itin);
% Anomali Gravitasi Data sintetik
[infileextn,inpath] = uigetfile('*.dat','Select grav_obs...');
grav=load(infileextn);
[md nd]=size(grav);
f=grav(:,2);        
f_space=grav(:,1);

% Parameter Model Blok 2D
nlay=itin(1,1);   % Banyak grid
nxg=itin(1,2);    % Banyak cell lateral
nzg=itin(1,3);    % Banyak cell vetikal
dx =itin(1,4);    % Dimensi cell lateral (m)
dh =itin(1,5);    % Dimensi cell vetikal (m)
model=[nxg nzg dx dh];
l(1,1:nlay)=itin(2,1:nlay); % Batas bawah pencarian
u(1,1:nlay)=itin(3,1:nlay); % Batas atas pencarian
x0=l+(u-l)./2;    % Tebakan model awal
Np=30;            % Jumlah populasi
Nb=itin(4,1:nlay);% Jumlah bit
Pc=1;             % Probabilitas crossover
Pm=0.01;          % Probabilitas mutasi
eta=0.1;          % Learning rate and the maximum
kmax=1000;        % Jumlah Iterasi
% [xo_gen,fo_gen,gap]=genetic(f,x0,l,u,Np,Nb,Pc,Pm,eta,kmax,model);

% =======================Inversi 10 Kali=========================%
color = ['r-','g-','r-','c-','m-','y-','k-','w-','c-','m-'];
xo = [];
foo = [];
gapa = [];
for i=1:10
    [xo_gen,fo_gen,gap]=genetic(f,x0,l,u,Np,Nb,Pc,Pm,eta,kmax,model);
    xo(i,:) =xo_gen; 
    foo(i,:) =fo_gen;
    gapa(i,:) =gap;
    [gm_cal]=forward_gravity(xo_gen,model); %Fungsi Forward
    VV2=reshape(xo_gen,nxg,nzg);
    VV=VV2';
    zSA1 = [0:dh:(nzg-1)*dh]';
    zSA = zSA1+10/2;
    
    figure(1)
    plot(gap,color(i));hold on
    title('Kurva Misfit Rata-rata','fontweight','bold','fontsize',8)
    ylabel('Misfit Rata-rata','fontsize',7);
    xlabel('Iterasi','fontsize',7);
    figure(2)
    subplot(2,1,1); plot(f_space,f,'s','LineWidth',5.0,'MarkerSize',3);hold on
    subplot(2,1,1); plot(f_space,gm_cal,color(i),'LineWidth',2.0,'MarkerSize',3);
    % xlim([0 max(x)]); 
    title('Hasil Inversi 10 Kali GA','fontweight','bold','fontsize',8)
    ylabel('Anomali Medan Gravitasi [mGal]','fontsize',7);
    xlabel('Spasi [m]','fontsize',7);
    % set(gca,'XLim',[0 21700])
    set(gca,'fontsize',7);
%     hleg= legend('G-obs','Inversi GA','Location','northeast');
%     set(hleg,'FontAngle','italic','FontSize',6);

    subplot(2,1,2); imagesc(f_space,zSA,VV);
    set(gca,'XAxisLocation','top','fontsize',7,'XMinorTick','on');
    ylabel('Kedalaman [m]','fontsize',7);
    % xlabel('Spasi [m]');

    colorbar('horiz');
    colormap('default');
    grid on
    % set(gca, 'GridLineStyle', '-');
    set(gca,'XGrid','on','YGrid','on','ZGrid','on')
    grid(gca,'minor');

end

rata2 = [];
for i=1:nxg*nzg
    rata2(i)=mean(xo(:,i));
end
mas2=reshape(rata2,nxg,nzg);
mas=mas2';


hasil = forward_gravity(rata2,model);
figure(3)
subplot(2,1,1); plot(f_space,f,'r-o','LineWidth',2.0,'MarkerSize',3);hold on
subplot(2,1,1); plot(f_space,hasil,'b-o','LineWidth',2.0,'MarkerSize',3);
% xlim([0 max(x)]); 
title('Hasil Inversi 10 Kali GA Rata-rata','fontweight','bold','fontsize',8)
ylabel('Anomali Medan Gravitasi [mGal]','fontsize',7);
xlabel('Spasi [m]','fontsize',7);
% set(gca,'XLim',[0 21700])
set(gca,'fontsize',7);
hleg= legend('G-obs','Inversi GA','Location','northeast');
set(hleg,'FontAngle','italic','FontSize',6);

subplot(2,1,2); imagesc(f_space,zSA,mas);
set(gca,'XAxisLocation','top','fontsize',7,'XMinorTick','on');
ylabel('Kedalaman [m]','fontsize',7);
% xlabel('Spasi [m]');

colorbar('horiz');
colormap('default');
grid on
% set(gca, 'GridLineStyle', '-');
set(gca,'XGrid','on','YGrid','on','ZGrid','on')
grid(gca,'minor');
