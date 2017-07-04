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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
userpath('reset');
% addpath(genpath('/home/knife/dag_matlab'))
addpath(genpath('/home/ddagnino/dag_matlab'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
[ strat, freq, iter ] = get_inv();

pre = true;
% pre = false;

exa = false;
exa1 = 5;
exa2 = 5;
exa3 = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path = './';
% path = '../data_output_2a/';
iter_name = ['_strat_',int2str(strat),'_freq_',num2str(freq,'%3.2f'),'_iter_freq_',int2str(iter)];
[ sh_o, n_rec,nt, dt, check, time, sh_o_min,sh_o_max ] = read_data( [path,'trace/tr',iter_name] );
[ sh_s, ~,~, ~, ~, ~, sh_s_min,sh_s_max ] = read_data( [path,'trace/tr_s',iter_name] );
[ adj , ~,~, ~, ~, ~, adj_min,adj_max   ] = read_data( [path,'adj_sou/adj_sou',iter_name] );
[ prec, ~,~, ~, ~, ~, prec_min,prec_max ] = read_data( [path,'adj_sou/prec_data',iter_name] );

tf = max(time);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % AGC
% t_agc = 0.5/freq;
% it0 = 1*ones(n_rec,1);
% it1 = nt*ones(n_rec,1);
% % [ sh_o, ~ ] = agc( sh_o, t_agc, tf, it0, it1 );
% % [ sh_s, ~ ] = agc( sh_s, t_agc, tf, it0, it1 );
% [ sh_o, ~ ] = agc_reverse( sh_o, t_agc, tf, it0, it1 );
% [ sh_s, ~ ] = agc_reverse( sh_s, t_agc, tf, it0, it1 );
% [ sh_o, ~ ] = agc_reverse( sh_o, t_agc, tf, it0, it1 );
% [ sh_s, ~ ] = agc_reverse( sh_s, t_agc, tf, it0, it1 );

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% slopes = zeros(4,1);
% amps = zeros(4,1);
% slopes(1) = 1./5500.;   amps(1) = 0.;
% slopes(2) = 1./5000.;   amps(2) = 1.;
% slopes(3) = 1./1450.;   amps(3) = 1.;
% slopes(4) = 1./1350.;   amps(4) = 0.;
% % bias = 1/1500;
% bias = 0.;
% dx = 25;
% 
% sh_o = dipfilter( sh_o, dt, nt, slopes, amps, bias, dx );
% sh_s = dipfilter( sh_s, dt, nt, slopes, amps, bias, dx );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
sh_o = sh_o';
sh_s = sh_s';
[ sh_o ] = filt_butterworth_all( sh_o, dt, 3, 0.5, 2*freq, 6, true );
[ sh_s ] = filt_butterworth_all( sh_s, dt, 3, 0.5, 2*freq, 6, true );
sh_o = sh_o';
sh_s = sh_s';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% adj = prec.*(sh_s - sh_o);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
x_rec0 = 250/1000;
dx_off = 25/1000;
line = x_rec0 + dx_off*(0:n_rec-1);

%
t0 = 0;
t1 = max(time);

% 
if pre
  sh_o = prec.*sh_o;
  sh_s = prec.*sh_s;
end

% 
u1max = max(max(abs(sh_o)));
u2max = max(max(abs(sh_s)));
u3max = max(abs([adj_min,adj_max]));
if exa
  u1max = u1max/exa1;
  u2max = u2max/exa2;
  u3max = u3max/exa3;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure(5); clf(h); set(h,'OuterPosition',pos1);

% 
tl1 = 1.75;
tl2 = 1.95;
tl3 = 2.10;

% 
dist_sou_rec = 187.9; dist_GPS_sou = 0.; doff = 25;
c_water = 1500;
dist = dist_GPS_sou + dist_sou_rec + doff*((1:n_rec)-1);
line_dw = dist/c_water;

% 
dist_sou_rec = 187.9; dist_GPS_sou = 0.; doff = 25;
vel = 2190;
dist = dist_GPS_sou + dist_sou_rec + doff*((1:n_rec)-1);
line_v = 0.09 + dist/vel;

% 
dist_sou_rec = 187.9; dist_GPS_sou = 0.; doff = 25;
vel = 2190;
dist = dist_GPS_sou + dist_sou_rec + doff*((1:n_rec)-1);
line_v2 = 0.39 + dist/vel;

subplot(131);
% pcolor( line, time, sh_o ); hold on;
pcolor( (1:n_rec), time, sh_o' ); hold on;
colorbar('location','SouthOutside'); set(gca,'linewidth',1,'fontsize',16); shading flat; axis('ij');
title('sh_o','fontsize',12); xlabel('x_{Rec} (km)','fontsize',12); ylabel('t (s)','fontsize',12);
% xlim([1 240]); 
ylim([t0 t1]); caxis(u1max*[-1 1]);
% p = plot([1 n_rec],[tl1 tl1],'-k'); set(p,'LineWidth',2);
% p = plot([1 n_rec],[tl2 tl2],'-k'); set(p,'LineWidth',2);
% p = plot([1 n_rec],[tl3 tl3],'-k'); set(p,'LineWidth',2);
p = plot((1:n_rec),line_dw,'-k'); set(p,'LineWidth',3);
p = plot((1:n_rec),line_v,'-k'); set(p,'LineWidth',3);
p = plot((1:n_rec),line_v2,'-k'); set(p,'LineWidth',2);
p = plot((1:n_rec),line_v2+0.25,'-k'); set(p,'LineWidth',2);

subplot(132);
% pcolor( line, time, sh_s ); hold on;
pcolor( (1:n_rec), time, sh_s' ); hold on;
colorbar('location','SouthOutside'); set(gca,'linewidth',1,'fontsize',16); shading flat; axis('ij');
title('sh_s','fontsize',12); xlabel('x_{Rec} (km)','fontsize',12); ylabel('t (s)','fontsize',12);
% xlim([1 240]); 
ylim([t0 t1]); caxis(u2max*[-1 1]);
% p = plot([1 n_rec],[tl1 tl1],'-k'); set(p,'LineWidth',2);
% p = plot([1 n_rec],[tl2 tl2],'-k'); set(p,'LineWidth',2);
% p = plot([1 n_rec],[tl3 tl3],'-k'); set(p,'LineWidth',2);
p = plot((1:n_rec),line_dw,'-k'); set(p,'LineWidth',3);
p = plot((1:n_rec),line_v,'-k'); set(p,'LineWidth',3);
p = plot((1:n_rec),line_v2,'-k'); set(p,'LineWidth',2);
p = plot((1:n_rec),line_v2+0.25,'-k'); set(p,'LineWidth',2);

subplot(133);
% pcolor( line, time, adj ); hold on;
pcolor( (1:n_rec), time, adj' ); hold on;
% pcolor( (1:n_rec), time, -(sh_o-sh_s)' ); hold on;
colorbar('location','SouthOutside'); set(gca,'linewidth',1,'fontsize',16); shading flat; axis('ij');
title('adj','fontsize',12); xlabel('x_{Rec} (km)','fontsize',12); ylabel('t (s)','fontsize',12);
% xlim([1 240]); 
ylim([t0 t1]); caxis(u3max*[-1 1]);
% p = plot([1 n_rec],[tl1 tl1],'-k'); set(p,'LineWidth',2);
% p = plot([1 n_rec],[tl2 tl2],'-k'); set(p,'LineWidth',2);
% p = plot([1 n_rec],[tl3 tl3],'-k'); set(p,'LineWidth',2);
p = plot((1:n_rec),line_dw,'-k'); set(p,'LineWidth',3);
p = plot((1:n_rec),line_v,'-k'); set(p,'LineWidth',3);
p = plot((1:n_rec),line_v2,'-k'); set(p,'LineWidth',2);
p = plot((1:n_rec),line_v2+0.25,'-k'); set(p,'LineWidth',2);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h = figure(51); clf(h); set(h,'OuterPosition',pos);
% 
% % ir = 1;
% % ir = 50;
% % ir = 100;
% % ir = 150;
% ir = 200;
% 
% plot( time, sh_s(:,ir)/max(abs([sh_s_min,sh_s_max])), ...
%       time, sh_o(:,ir)/max(abs([sh_o_min,sh_o_max])), ...
%       time, -adj(:,ir)/max(abs([adj_min,adj_max])) );
% legend('synt','obs','adj',3,'Location','NorthEast');


