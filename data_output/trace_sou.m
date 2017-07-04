%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;

set(0,'Units','pixels');
scnsize = get(0,'ScreenSize');
pos = [50, 200, 1870, 890];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
strat = 0;
freq = 2.00;
iter = 1;

tf = 17.5;
tf_s = 17.5;

tw = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tr_obs = load(['trace/tr_strat_',int2str(strat),'_freq_',num2str(freq,'%3.2f'),'_iter_freq_',int2str(iter),'.txt']);
tr_syn = load(['trace/tr_s_strat_',int2str(strat),'_freq_',num2str(freq,'%3.2f'),'_iter_freq_',int2str(iter),'.txt']);
prec = load(['adj_sou/prec_data_strat_',int2str(strat),'_freq_',num2str(freq,'%3.2f'),'_iter_freq_',int2str(iter),'.txt']);

n_rec = length(tr_obs(1,:));
nt = length(tr_obs(:,1));
dt = tf/(nt-1);
time = dt*(0:nt-1);

% calib = ones(n_rec,2);
for ir=1:n_rec
  calib(ir,2) = max(abs(tr_syn(:,ir)))/max(abs(tr_obs(:,ir)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
coef  = 1;
shift = 0.;

% ntr = [ 220, 240, 260, 280, 400, 643 ];
ntr = [ 220, 240, 249, 250, 260, 280 ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure(2); clf(h); set(h,'OuterPosition',pos);

subplot(321);
Amax = coef*max(max(calib(ntr(1),2)*tr_obs(:,ntr(1)))); Amin = coef*min(min(calib(ntr(1),2)*tr_obs(:,ntr(1))));
plot( time, calib(ntr(1),2)*tr_obs(:,ntr(1)), time+shift, coef*tr_syn(:,ntr(1)),'linewidth',2 );
% ylim([Amin Amax]);
% xlim([0 tf]);
xlim(3*[1 1] + [0 tw]);
legend('1 tr_{obs}','1 tr_{syn}',2,'Location','NorthEast');

subplot(322);
Amax = coef*max(max(calib(ntr(2),2)*tr_obs(:,ntr(2)))); Amin = coef*min(min(calib(ntr(2),2)*tr_obs(:,ntr(2))));
plot( time, calib(ntr(2),2)*tr_obs(:,ntr(2)), time+shift, coef*tr_syn(:,ntr(2)),'linewidth',2 );
% ylim([Amin Amax]);
% xlim([0 tf]);
xlim(1.5*[1 1] + [0 tw]);
legend('2 tr_{obs}','2 tr_{syn}',2,'Location','NorthEast');

subplot(323);
Amax = coef*max(max(calib(ntr(3),2)*tr_obs(:,ntr(3)))); Amin = coef*min(min(calib(ntr(3),2)*tr_obs(:,ntr(3))));
plot( time, calib(ntr(3),2)*tr_obs(:,ntr(3)), time+shift, coef*tr_syn(:,ntr(3)),'linewidth',2 );
% ylim([Amin Amax]);
% xlim([0 tf]);
xlim(1*[1 1] + [0 tw]);
legend('3 tr_{obs}','3 tr_{syn}',2,'Location','NorthEast');

subplot(324);
Amax = coef*max(max(calib(ntr(4),2)*tr_obs(:,ntr(4)))); Amin = coef*min(min(calib(ntr(4),2)*tr_obs(:,ntr(4))));
plot( time, calib(ntr(4),2)*tr_obs(:,ntr(4)), time+shift, coef*tr_syn(:,ntr(4)),'linewidth',2 );
% ylim([Amin Amax]);
% xlim([0 tf]);
xlim(1*[1 1] + [0 tw]);
legend('4 tr_{obs}','4 tr_{syn}',2,'Location','NorthEast');

subplot(325);
Amax = coef*max(max(calib(ntr(5),2)*tr_obs(:,ntr(5)))); Amin = coef*min(min(calib(ntr(5),2)*tr_obs(:,ntr(5))));
plot( time, calib(ntr(5),2)*tr_obs(:,ntr(5)), time+shift, coef*tr_syn(:,ntr(5)),'linewidth',2 );
% ylim([Amin Amax]);
% xlim([0 tf]);
xlim(1*[1 1] + [0 tw]);
legend('5 tr_{obs}','5 tr_{syn}',2,'Location','NorthEast');

subplot(326);
Amax = coef*max(max(calib(ntr(6),2)*tr_obs(:,ntr(6)))); Amin = coef*min(min(calib(ntr(6),2)*tr_obs(:,ntr(6))));
plot( time, calib(ntr(6),2)*tr_obs(:,ntr(6)), time+shift, coef*tr_syn(:,ntr(6)),'linewidth',2 );
% ylim([Amin Amax]);
% xlim([0 tf]);
xlim(2*[1 1] + [0 tw]);
legend('6 tr_{obs}','6 tr_{syn}',2,'Location','NorthEast');




