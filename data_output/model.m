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
[ vel, nx,ny,nz, dx, x,y,z, v_min,v_max ] = read_model( '../data_input/model/vP' );
% [ vel, nx,ny,nz, dx, x,y,z, v_min,v_max ] = read_model( '../data_input/model/vP_FAST' );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure(1); clf(h); set(h,'OuterPosition',pos2);

%--------------------------------------------------------------%
pcolor( x, y, vel );
colorbar('location','EastOutside'); set(gca,'linewidth',1,'fontsize',16); shading flat; axis('ij');
% caxis([1500 4200]);
% caxis([v_min v_max]);
% ylim([0 3000]);
title('Vp','fontsize',16); ylabel('z (km)','fontsize',16); xlabel('x (km)','fontsize',16);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure(2); clf(h); set(h,'OuterPosition',pos2);

%--------------------------------------------------------------%
xp = [10:10:50];
ixp = 1 + floor(xp/dx);

CM = jet(length(xp));
for ip=1:length(xp)
  plot( y, vel(:,ixp(ip)), 'color', CM(ip,:),'linewidth',2 ); hold on;
  legendInfo{ip} = num2str(xp(ip));
end
set(gca,'linewidth',2,'fontsize',16);
% caxis([v_min v_max]);
title('Vp','fontsize',16); ylabel('v (km/s)','fontsize',16); xlabel('x (km)','fontsize',16);
legend(legendInfo,'Location','NorthEast');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sum(vel(:,1)/nz
% sum(vel(:,floor(nx/2))/nz
% sum(vel(:,nx)/nz

two_t_depth = zeros(length(xp),1);
for ip=1:length(xp)
  two_t_depth(ip) = 2*sum( 1000*dx*ones(nz,1)./vel(:,ixp(ip)) );
end

disp('two way time = ');
disp(num2str(two_t_depth));

