clear functions;
close all; clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(0,'Units','pixels');
scnsize = get(0,'ScreenSize');
pos = [50, 200, 1870, 890];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tr_obs = load('../data_input/shot_gather/aux.txt');
Amax = 0.99*max(max(tr_obs));
Amin = 0.99*min(min(tr_obs));

NumRec = length(tr_obs(1,:));
nt = length(tr_obs(:,1));

dx_off = 150.;
x_rec0 = 0;
tf = 30.;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure(10); clf(h); set(h,'OuterPosition',pos);

pcolor( x_rec0/1000 + dx_off*(0:NumRec-1)/1000, tf*(0:nt-1)/(nt-1), tr_obs(1:nt,:) );
% pcolor( 1:NumRec, tf*(0:nt-1)/(nt-1), tr_obs(1:nt,:) );
colorbar('location','EastOutside'); set(gca,'linewidth',1,'fontsize',16); shading flat; axis('ij');
caxis([Amin Amax]);
ylim([3 17]);
title('p'); xlabel('x_{Rec} (km)'); ylabel('t (s)');


% for ir=1:NumRec
%   vect(ir) = sum( abs(tr_obs(:,ir)) )/NumRec;
% end
% 
% plot(1:NumRec,vect)

