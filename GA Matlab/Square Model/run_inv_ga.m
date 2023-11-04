%------------------------------------------------------
% Inversion of Graity Data using GA
% Created by asrin dan Al rubaiyn
% GA alfgorithm adopted from Yang, W.Y., dkk., (2020)
% Forward adopted from Last Kubik
%------------------------------------------------------
clc; clear all; close all;

% Calling the input model parameters
[infileextn,inpath] = uigetfile('*.xlsx','Select File InputInversion.xls...');
itin=xlsread(infileextn);
[ni mi]=size(itin);
% Calling the Gravity Anomaly
[infileextn,inpath] = uigetfile('*.dat','Select grav_obs...');
grav=load(infileextn);
[md nd]=size(grav);
f=grav(:,2); % observation data
f_space=grav(:,1);

% Parameter Model Blok 2-D
nlay=itin(1,1);   %Banyak grid
nxg=itin(1,2);    %Banyak cell lateral
nzg=itin(1,3);    %Banyak cell vetikal
dx =itin(1,4); %Dimensi cell lateral (m)
dh =itin(1,5); %Dimensi cell vetikal (m)
model=[nxg nzg dx dh];

l(1,1:nlay)=itin(2,1:nlay);% lower earch of paramater
u(1,1:nlay)=itin(3,1:nlay); % upper search of paramter

x0=l+(u-l)./2; % initial guess
%do_genetic.m
Np=30; % Population size
Nb=itin(4,1:nlay);% Numbers of bits for representing each variable
% Nb=[10 10]; % Numbers of bits for representing each variable
Pc=1; Pm=0.01; % Probability of crossover/mutation
eta=0.1; kmax=100; % Learning rate and the maximum # of iterations
[xo_gen,fo_gen]=genetic(f,x0,l,u,Np,Nb,Pc,Pm,eta,kmax,model);
[gm_cal]=forward_gravity(xo_gen,model); %Fungsi Forward
VV2=reshape(xo_gen,nxg,nzg);
VV=VV2';
zSA1 = [0:dh:(nzg-1)*dh]';
zSA = zSA1+10/2;

figure(2)
subplot(2,1,1); plot(f_space,f,'b-o','LineWidth',2.0,'MarkerSize',3);hold on
subplot(2,1,1); plot(f_space,gm_cal,'r-o','LineWidth',2.0,'MarkerSize',3);
% xlim([0 max(x)]); 
title('Penampang Hasil Inversi GA 2D Gravitasi','fontweight','bold','fontsize',8)
ylabel('Anomali Medan Gravitasi [mGal]','fontsize',7);
xlabel('Spasi [m]','fontsize',7);
% set(gca,'XLim',[0 21700])
set(gca,'fontsize',7);
hleg= legend('g-obs','Inversi GA','Location','NorthEast');
set(hleg,'FontAngle','italic','FontSize',6);

subplot(2,1,2); imagesc(f_space,zSA,VV);
set(gca,'XAxisLocation','top','fontsize',7,'XMinorTick','on');
ylabel('Kedalaman [m]','fontsize',7);
% xlabel('Spasi [m]');

colorbar('horiz');
colormap('default');
grid on
set(gca, 'GridLineStyle', '-');
set(gca,'XGrid','on','YGrid','on','ZGrid','on')
grid(gca,'minor');
