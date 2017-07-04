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

pos1 = [1, 1, 1920, 1080];
pos2 = [1, 221, 1920, 860];
pos3 = [1, 521, 1920, 560];

pos1b = [1920+1, 1, 1920, 1080];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
[ strat, freq, iter ] = get_inv();

y0 = 0.;
y1 = 3.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
[ vp_f2, nx_f2,ny_f2,nz_f2, dx_f2, x_f2,y_f2,z_f2, vf_min2,vf_max2 ] = read_model( '../data_input/model/vP_FAST_3' );
% [ vp_f, nx_f,ny_f,nz_f, dx_f, x_f,y_f,z_f, vf_min,vf_max ] = read_model( '../data_input/model/vP_FAST_3' );
[ vp_f, nx_f,ny_f,nz_f, dx_f, x_f,y_f,z_f, vf_min,vf_max ] = read_model( ['model/vP_strat_',int2str(strat),'_freq_',num2str(freq,'%3.2f'),'_iter_freq_',int2str(iter)] );
vp_f2  = vp_f2/1000;
vp_f  = vp_f/1000;
disp(['[vf_min,vf_max] = ',num2str([vf_min,vf_max]/1000)]);
% vi_min = vf_min/1000;
% vi_max = vf_max/1000;
vi_min = 1.5;
vi_max = 4.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
[ mig, nx_m,ny_m,nz_m, dx_m,dy_m,dz_m, x_m,y_m,z_m, mig_min,mig_max ] = read_migration( '../../depth_migration_Ranero' );

x_m=x_m/1000;
y_m=y_m/1000;
[ vp_mig ] = vp_plus_migration( 1.0, vp_f, nx_f,ny_f, x_f,y_f, mig, nx_m,ny_m, x_m,y_m );

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% [ illu, nx,ny,nz, dx, x,y,z, illu_min,illu_max ] = read_model( ['illu/illu_strat_',int2str(strat),'_freq_',num2str(freq, '%3.2f') ,'_iter_freq_',int2str(iter)] );
% 
% % 
% illu = illu/illu_max;
% illu(illu>0.05) = 1;
% illu(illu<1) = 0;
% 
% vp_f = vp_f.*illu;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clc;

%--------------------------------------------------------------%
hall = figure('Name','10','renderer','opengl'); clf(hall); set(hall,'OuterPosition',pos1);

% 
pcolor( 135.74+x_f, y_f, vp_f ); hold on;
colorbar('location','EastOutside'); set(gca,'linewidth',1,'fontsize',16); shading flat; axis('ij');
caxis([vi_min vi_max]);
freezeColors;
ylim([0 max(y_m)]);
% xlim([150 170]);

% Overlay semitransparent ice speed:
h = pcolor( 135.74+x_m, y_m, mig/mig_max );
caxis(0.1*[-1 1]);
colormap(gray(256)); c=colorbar('location','WestOutside'); colorbar(c,'off');
shading interp; set(h,'facealpha',.4);


%--------------------------------------------------------------%
hall = figure('Name','11','renderer','opengl'); clf(hall); set(hall,'OuterPosition',pos1b);

% 
pcolor( 135.74+x_f2, y_f2, vp_f2 ); hold on;
colorbar('location','EastOutside'); set(gca,'linewidth',1,'fontsize',16); shading flat; axis('ij');
caxis([vi_min vi_max]);
freezeColors;
ylim([0 max(y_m)]);
% xlim([150 170]);

% Overlay semitransparent ice speed:
h = pcolor( 135.74+x_m, y_m, mig/mig_max );
caxis(0.1*[-1 1]);
colormap(gray(256)); c=colorbar('location','WestOutside'); colorbar(c,'off');
shading interp; set(h,'facealpha',.4);





