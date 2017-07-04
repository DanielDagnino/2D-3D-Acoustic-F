%{
!*********************************************************************/
!** This code has been done in the Barcelona Center for Subsurface 
!** Imaging (BCSI).
!** Goal: Set of tools to analyse the FWI results.
!** Authors: Daniel Dagnino.
!*********************************************************************/
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(0,'Units','pixels') 
scnsize = get(0,'ScreenSize');
pos = [520, 400, 1800, 950];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
strat = 1;
freq = 10.00;
iter = 1;

dx = 15/1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pre = load(['illu/prec_mod_strat_',int2str(strat),'_freq_',num2str(freq, '%3.2f') ,'_iter_freq_',int2str(iter),'.txt']);
% pre = load(['illu/illu_strat_',int2str(strat),'_freq_',num2str(freq, '%3.2f') ,'_iter_freq_',int2str(iter),'.txt']);

ny = length(pre(:,1));
nx = length(pre(1,:));

x = dx*(0:nx-1);
y = dx*(0:ny-1);

% premax = 0.1*max(max(pre));
% premin = min(min(pre));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure(4); clf(h); set(h,'OuterPosition',pos);

% for ix=1:nx
%   pre(:,ix) = (1+(0:ny-1)/(ny-1))*pre(:,ix);
% end

% for ix=1:nx
%   pre(:,ix) = (1+1000*(ix-1)/(nx-1))*pre(:,ix);
% end

pcolor( x, y, pre );
colorbar('location','EastOutside'); set(gca,'linewidth',1,'fontsize',16); 
% caxis([premin premax]); 
shading flat; axis('ij');
title('premination','fontsize',12); xlabel('x (km)','fontsize',12); ylabel('z (km)','fontsize',12);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h = figure(4); clf(h); set(h,'OuterPosition',pos);
% 
% plot( y, pre(:,20) );


