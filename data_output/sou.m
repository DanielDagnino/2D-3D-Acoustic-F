%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(0,'Units','pixels');
scnsize = get(0,'ScreenSize');
pos = [50, 200, 1870, 890];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
strat = 1;
freq = 10.00;
iter = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source = load(['source/sou_strat_',int2str(strat),'_freq_',num2str(freq,'%3.2f'),'_iter_freq_',int2str(iter),'.txt']);

tf = 4.90;

nt = length(source(:,1));
dt = tf/nt;
time = tf*(0:nt-1)/(nt-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure(20); clf(h); set(h,'OuterPosition',pos);
plot( time, source(:,1));
% xlim([0 0.8]);
set(gca,'linewidth',1,'fontsize',12);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% source = load(['../data_input/sou_rec/swl_aux_1_1Hz.txt']);
% 
% tf = 3.;
% 
% nt = length(source(:,1));
% dt = tf/nt;
% time = tf*(0:nt-1)/(nt-1);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h = figure(20); clf(h); set(h,'OuterPosition',pos);
% plot( time, source(:,1));
% % xlim([0 0.8]);
% set(gca,'linewidth',1,'fontsize',12);
