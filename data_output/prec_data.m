%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(0,'Units','pixels');
scnsize = get(0,'ScreenSize');
pos = [520, 100, 1800, 1000];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
strat = 0;
freq = 2.00;
iter = 1;

tfin = 21.;
dx_off = 150/1000;
x_rec0 = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
adj = load(['adj_sou/prec_data_strat_',int2str(strat),'_freq_',num2str(freq,'%3.2f'),'_iter_freq_',int2str(iter),'.txt']);

coef  =  0.4;
Amax = coef*max(max(adj)); Amin = coef*min(min(adj));

nt = length(adj(:,1));
n_rec = length(adj(1,:));
dt = tf/nt;

h = figure(41); clf(h); set(h,'OuterPosition',pos);

set(h,'OuterPosition',pos);
pcolor( x_rec0/1000 + dx_off*(0:n_rec-2)/1000, tfin*(0:nt-1)/(nt-1), adj(:,2:n_rec) );
colorbar('location','EastOutside'); set(gca,'linewidth',1,'fontsize',16); shading flat; axis('ij');
% caxis([Amin Amax]);
title('adj','fontsize',12); ylim([0 tfin]); xlabel('x_{Rec} (km)','fontsize',12); ylabel('t (s)','fontsize',12);


