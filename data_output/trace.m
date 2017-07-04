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

pos = [1, 201, 1920, 880];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
[ strat, freq, iter ] = get_inv();

% pre = false;
pre = true;

exa = true;
exa1 = 5;
exa2 = 5;
exa3 = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iter_name = ['_strat_',int2str(strat),'_freq_',num2str(freq,'%3.2f'),'_iter_freq_',int2str(iter)];
[ sh_o, n_rec,nt, dt, check, time, sh_o_min,sh_o_max ] = read_data( ['trace/tr',iter_name] );
[ sh_s, ~,~, ~, ~, ~, sh_s_min,sh_s_max ] = read_data( ['trace/tr_s',iter_name] );
[ adj , ~,~, ~, ~, ~, adj_min,adj_max   ] = read_data( ['adj_sou/adj_sou',iter_name] );
[ prec, ~,~, ~, ~, ~, prec_min,prec_max ] = read_data( ['adj_sou/prec_data',iter_name] );

% sh_o = sh_o';
% sh_s = sh_s';
% adj = adj';
% prec = prec';

tf = max(time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AGC
t_agc = 0.5/freq;
it0 = 1*ones(n_rec,1);
it1 = nt*ones(n_rec,1);
[ sh_o, ~ ] = agc_reverse( sh_o, t_agc, tf, it0, it1 );
[ sh_s, ~ ] = agc_reverse( sh_s, t_agc, tf, it0, it1 );
[ sh_o, ~ ] = agc_reverse( sh_o, t_agc, tf, it0, it1 );
[ sh_s, ~ ] = agc_reverse( sh_s, t_agc, tf, it0, it1 );

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
calib = ones(n_rec,2);

%   calib = load(['calib/a_calib_sou_strat_',int2str(strat),'_freq_',num2str(freq,'%3.2f'),'_iter_freq_',int2str(iter),'.txt']);
%   calib = load(['calib/a_calib_adj_strat_',int2str(strat),'_freq_',num2str(freq,'%3.2f'),'_iter_freq_',int2str(iter),'.txt']);
%   calib = load(['calib/a_calib_mf_strat_',int2str(strat),'_freq_',num2str(freq,'%3.2f'),'_iter_freq_',int2str(iter),'.txt']);

for ir=1:n_rec
  calib(ir,2) = max(abs(sh_o(ir,:)))/max(abs(sh_s(ir,:)));
end

% for ir=1:n_rec
%   calib(ir,2) = max(abs(sh_s(ir,:)))/max(abs(sh_o(ir,:)));
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
coef  = 1;
shift = 0.;

% 
tw = tf;
% ntr = [ 1, 5, 10, 15, 30, 240 ];
% tp = 0.0*[ 1, 1, 1, 1, 1, 1 ];

% ntr = [ 40, 60, 80, 100, 105, 113 ];
% tp = 0.5*[ 1, 1, 1, 1, 1, 1 ];

ntr = [ 40, 60, 80, 100, 120, 133 ];
tp = 0.5*[ 1, 1, 1, 1, 1, 1 ];

% ntr = [ 50, 100, 120, 130, 140, 193 ];
% tp = 0.5*[ 1, 1, 1, 1, 1, 1 ];

% ntr = [ 140, 160, 180, 200, 220, 240 ];
% tp = 1.5*[ 1, 1, 1, 1, 1, 1 ];

% % ntr = [ 40, 140, 160, 186, 200, 210 ];
% ntr = [ 60, 100, 125, 150, 180, 210 ];
% tp = 1.0*[ 1, 1, 1, 1, 1, 1 ];

tpp = tp'*[1 1] + ones(6,1)*[0 tw];
tpp(:,2) = min(tpp(:,2),tf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure(21); clf(h); set(h,'OuterPosition',pos);

for k=1:6
  
  subplot(3,2,k);
  
  itr = ntr(k);
  cal = calib(itr,2);
%   Amax = coef*max(max(cal*sh_o(itr,:))); Amin = coef*min(min(cal*sh_o(itr,:)));
  Amax = coef*max(max(cal*sh_s(itr,:))); Amin = coef*min(min(cal*sh_s(itr,:)));
  
%   plot( time, cal*sh_o(itr,:), time+shift, coef*sh_s(itr,:),'linewidth',2 );
%   plot( time, sh_o(itr,:), time+shift, coef*cal*sh_s(itr,:),'linewidth',2 );
%   legend([num2str(k),' tr_{obs}'],[num2str(k),' tr_{syn}'],2,'Location','NorthEast');
  
  a = cal*sh_s(itr,:);
  b = sh_o(itr,:);
  a = a./max(max(abs(a)));
  b = b./max(max(abs(b)));
%   a(floor(3.3/dt):nt) = 0;
%   b(floor(3.3/dt):nt) = 0;
  a(isnan(a))=0;
  b(isnan(b))=0;
  [c, i] = max( xcorr(a,b)/(sqrt(max(abs(xcorr(a,a))))*sqrt(max(abs(xcorr(b,b))))) );
  [c, (tf-dt*(i-1))/(0.5/freq)]
%   shift = (tf-dt*(i-1));
%   coef = max(xcorr(a,b))/max(abs(xcorr(b,b))) ;
  
  adj = adj/max(max(abs(adj)));
  plot( time+shift,coef*a, time,b, time,adj(itr,:), 'linewidth',2 );
  legend([num2str(k),' tr_{syn}'],[num2str(k),' tr_{obs}'],[num2str(k),' adj'],3,'Location','NorthEast');
  
  % ylim([Amin Amax]);
  % xlim([0 tf]);
  xlim(tpp(k,:));
  
end



