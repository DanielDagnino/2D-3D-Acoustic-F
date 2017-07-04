%{
!*********************************************************************/
!** This code has been done in the Barcelona Center for Subsurface 
!** Imaging (BCSI).
!** Goal: Set of tools to analyse the FWI results.
!** Authors: Daniel Dagnino.
!*********************************************************************/
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% My Matlab script for the AGC of the SU.

function [ sh_o_agc, agc_factor ] = agc( sh_o, t_agc, tf )
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 
  agc_factor = zeros(size(sh_o));
  nt = length(sh_o(:,1));
  n_rec = length(sh_o(1,:));
  dt = tf/(nt-1);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for ir=1:n_rec

    %--------------------------------------------------------------%
    % 
    data = sh_o(:,ir)';

    %--------------------------------------------------------------%
    % 
    iwagc = 1+floor(t_agc/dt);

    % 
    agc_tr = zeros(1,nt);

    % compute initial window for first datum.
    sum = 0.0;
    for i=1:(iwagc+1)
      sum = sum + data(i)^2;
    end
    nwin = 2*iwagc+1;
    rms = sum/nwin;
    agc_tr(1) = 1/sqrt(rms);

    % ramping on.
    for i=2:(iwagc+1)
      sum = sum + data(i+iwagc)^2;
      nwin = nwin + 1;
      rms = sum/nwin;
      agc_tr(i) = 1/sqrt(rms);
    end

    % middle range -- full rms window.
    for i=(iwagc+2):(nt-iwagc)
      sum = sum + data(i+iwagc)^2 - data(i-iwagc)^2;
      rms = sum/nwin;
      agc_tr(i) = 1/sqrt(rms);
    end

    % ramping off.
    for i=(nt-iwagc+1):nt
      sum = sum - data(i-iwagc)^2;
      nwin = nwin - 1;
      rms = sum/nwin;
      agc_tr(i) = 1/sqrt(rms);
    end

    %--------------------------------------------------------------%
    % 
    agc_factor(:,ir) = agc_tr;

  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Apply AGC.
  sh_o_agc = agc_factor.*sh_o;
  
end




