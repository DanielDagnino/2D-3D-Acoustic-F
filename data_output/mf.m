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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
path = './';
% path = '../data_output_2a/';

% 
iter_name = ['_strat_',int2str(strat),'_freq_',num2str(freq,'%3.2f'),'_iter_freq_',int2str(iter)];
file1 = [ path,'resume/mf1', iter_name, '.txt' ];
file2 = [ path,'resume/mf2', iter_name, '.txt' ];
data1 = load(file1);
data2 = load(file2);

iter_name = ['_strat_',int2str(strat),'_freq_',num2str(freq,'%3.2f'),'_iter_freq_',int2str(1)];
file1 = [ path,'resume/mf1', iter_name, '.txt' ];
data0 = load(file1);

ns = length(data1(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure(1); clf(h); set(h,'OuterPosition',pos1);

%--------------------------------------------------------------%
subplot(211);

lineap = 50*((1:ns)-1)/1000;
line2p = 50*([1 ns]-1)/1000;

% plot( lineap, data1(:,1)/max(data1(:,1)), lineap, data2/max(data2) ); hold on;
plot( lineap, data1(:,1), lineap, data2, lineap, data1(:,1)-data2 ); hold on;
set(gca,'linewidth',1,'fontsize',16);
xlim([0 max(lineap)]);
legend('data1','data2','data1-data2',3,'Location','NorthEast');
title('Last','fontsize',12);

p = plot(line2p,[0 0],'-k'); set(p,'LineWidth',2);

%--------------------------------------------------------------%
subplot(212);

lineap = 50*((1:ns)-1)/1000;
line2p = 50*([1 ns]-1)/1000;

% plot( lineap, data1(:,1)/max(data1(:,1)), lineap, data2/max(data2) ); hold on;
plot( lineap, data0(:,1), lineap, data2, lineap, data0(:,1)-data2 ); hold on;
set(gca,'linewidth',1,'fontsize',16);
xlim([0 max(lineap)]);
legend('data1','data2','data1-data2',3,'Location','NorthEast');
title('Total','fontsize',12);

p = plot(line2p,[0 0],'-k'); set(p,'LineWidth',2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dif = data1(:,1)-data2;
% prec_wrong = dif(dif<0);
% sum(abs(prec_wrong))/sum(abs(dif))

dif = data1(:,1)-data2;
prec_wrong = dif(dif<0);
disp(['last = ', num2str(100*sqrt(sum(prec_wrong.^2)/sum(dif.^2))) ,'%']);

dif = data0(:,1)-data2;
prec_wrong = dif(dif<0);
disp(['total = ', num2str(100*sqrt(sum(prec_wrong.^2)/sum(dif.^2))) ,'%']);








