%{
!*********************************************************************/
!** This code has been done in the Barcelona Center for Subsurface 
!** Imaging (BCSI).
!** Goal: Set of tools to analyse the FWI results.
!** Authors: Daniel Dagnino.
!*********************************************************************/
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;
set(0,'Units','pixels');
scnsize = get(0,'ScreenSize');

pos = [520, 100, 1800, 1000];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
[ strat, freq, iter ] = get_inv();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
[ illu, nx,ny,nz, dx, x,y,z, illu_min,illu_max ] = read_model( ['illu/illu_strat_',int2str(strat),'_freq_',num2str(freq, '%3.2f') ,'_iter_freq_',int2str(iter)] );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure(3);
clf(h); set(h,'OuterPosition',pos);

illu = illu/illu_max;
% illu(illu>0.1) = 1;
% illu(illu<1) = 0;

pcolor( x, y, illu );
colorbar('location','EastOutside'); set(gca,'linewidth',1,'fontsize',16); 
shading flat; axis('ij');
title('Illumination','fontsize',12); xlabel('x (km)','fontsize',12); ylabel('z (km)','fontsize',12);
caxis([0 0.5]); 

% pcolor( x, y, log(illu) );
% colorbar('location','EastOutside'); set(gca,'linewidth',1,'fontsize',16); 
% shading flat; axis('ij');
% title('Illumination','fontsize',12); xlabel('x (km)','fontsize',12); ylabel('z (km)','fontsize',12);





