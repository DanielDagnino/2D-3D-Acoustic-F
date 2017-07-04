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

pos1b = [1, 1920+1, 1920, 1080];
pos2b = [1, 1920+221, 1920, 860];
pos3b = [1, 1920+521, 1920, 560];
pos1 = [1, 1, 1920, 1080];
pos2 = [1, 221, 1920, 860];
pos3 = [1, 521, 1920, 560];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
[ strat, freq, iter ] = get_inv();

y0 = 0.;
y1 = 3.0;

x0 = 135.74;
x0 = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% [ vp_i, nx_i, ny_i, nz_i, dx_i, x_i,y_i,z_i, vi_min,vi_max ] = read_model( '../data_input/model/vP_new' );
[ vp_i, nx_i, ny_i, nz_i, dx_i, x_i,y_i,z_i, vi_min,vi_max ] = read_model( '../data_input/model/vP_FAST_3' );
% [ vp_i, nx_i, ny_i, nz_i, dx_i, x_i,y_i,z_i, vi_min,vi_max ] = read_model( '../data_input/model/vP_FAST_2' );
% [ vp_i, nx_i, ny_i, nz_i, dx_i, x_i,y_i,z_i, vi_min,vi_max ] = read_model( '../data_input/model/vP' );
% [ vp_i, nx_i,ny_i,nz_i, dx_i, x_i,y_i,z_i, vi_min,vi_max ] = read_model( '../data_output/model/vp' );
vp_i  = vp_i/1000;
disp(['[vi_min,vi_max] = ',num2str([vi_min,vi_max]/1000)]);
vi_min = vi_min/1000;
vi_max = vi_max/1000;
vi_min = 1.5;
vi_max = 4.5;
% vi_max = 5.0;

% 
[ vp_f, nx_f,ny_f,nz_f, dx_f, x_f,y_f,z_f, vf_min,vf_max ] = read_model( ['model/vP_strat_',int2str(strat),'_freq_',num2str(freq,'%3.2f'),'_iter_freq_',int2str(iter)] );
vp_f  = vp_f/1000;
disp(['[vf_min,vf_max] = ',num2str([vf_min,vf_max]/1000)]);

% 
[ searc, nx_s,ny_s,nz_s, dx_s, x_s,y_s,z_s, s_min,s_max ] = read_model( ['grad/search_strat_',int2str(strat),'_freq_',num2str(freq,'%3.2f'),'_iter_freq_',int2str(iter)] );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Old grid.
x1_i = zeros(nx_i,ny_i);
y1_i = zeros(nx_i,ny_i);

for ix=1:nx_i
  x1_i(ix,:) = x_i(ix);
end

for iy=1:ny_i
  y1_i(:,iy) = y_i(iy);
end

% New grid.
xI_f = zeros(nx_f,ny_f);
yI_f = zeros(nx_f,ny_f);

for ix=1:nx_f
  xI_f(ix,:) = x_f(ix);
end

for iy=1:ny_f
  yI_f(:,iy) = y_f(iy);
end

vpI_i = interp2(x1_i',y1_i',vp_i,xI_f',yI_f','linear',NaN);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure(2); clf(h); set(h,'OuterPosition',pos1);

%--------------------------------------------------------------%
subplot(221);
% pcolor( x_i, y_i, fliplr(vp_i) );
% pcolor( x_i, y_i, vp_i );
pcolor( x0+xI_f, yI_f, vpI_i' );
colorbar('location','EastOutside'); set(gca,'linewidth',1,'fontsize',16); shading flat; axis('ij');
caxis([vi_min vi_max]);
% caxis([1470 1550]/1000);
title('Initial Vp (km/s)')
ylim([y0 y1]);

%--------------------------------------------------------------%
subplot(222);
% pcolor( x_f, y_f, fliplr(vp_f) );
pcolor( x0+x_f, y_f, vp_f );
colorbar('location','EastOutside'); set(gca,'linewidth',1,'fontsize',16); shading flat; axis('ij');
caxis([vi_min vi_max]);
title('Final Vp (km/s)')
ylim([y0 y1]);

%--------------------------------------------------------------%
subplot(223);
% pcolor( x_f, y_f, fliplr(searc) );
pcolor( x0+x_f, y_f, searc );
colorbar('location','EastOutside'); set(gca,'linewidth',1,'fontsize',16); shading flat; axis('ij');
% caxis([0.1 0.1]);
title('Search');
ylim([y0 y1]);

%--------------------------------------------------------------%
subplot(224);
pcolor( x0+x_f, y_f, vp_f-vpI_i );
% pcolor( x_f, y_f, 100*(vp_f-vpI_i)./vpI_i );
colorbar('location','EastOutside'); set(gca,'linewidth',1,'fontsize',16); shading flat; axis('ij');
% caxis(0.5*[-1 1]);
title('Final-Initial Vp (km/s)');
% ylim([y0 y1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure(21); clf(h); set(h,'OuterPosition',pos2);

%--------------------------------------------------------------%
xp = [13:10:53];
ixp = 1 + floor(xp/dx_f);

CM = jet(length(xp));
for ip=1:length(xp)
  plot( y_f, vp_f(:,ixp(ip)), 'color', CM(ip,:),'linewidth',2 ); hold on;
  legendInfo{ip} = num2str(xp(ip));
end
set(gca,'linewidth',2,'fontsize',16);
% caxis([v_min v_max]);
caxis([1 5.5]);
title('Vp','fontsize',16); ylabel('v (km/s)','fontsize',16); xlabel('x (km)','fontsize',16);
legend(legendInfo,'Location','NorthEast');


