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

function [ sh_o_agc, agc_factor ] = agc_reverse( sh_o, t_agc, tf, it0, it1 )
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 
  agc_factor = zeros(size(sh_o));
  nt = length(sh_o(1,:));
  n_rec = length(sh_o(:,1));
  dt = tf/(nt-1);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for ir=1:n_rec

    %--------------------------------------------------------------%
    % 
    data = sh_o(ir,:);

    %--------------------------------------------------------------%
    % 
    iwagc = 1+floor(t_agc/dt);

    % 
    agc_tr = zeros(1,nt);

    % compute initial window for first datum.
    sum = 0.0;
    for it=it1(ir):-1:(it1(ir)-iwagc)
      sum = sum + data(it)^2;
    end
    nwin = 2*iwagc+1;
    rms = sum/nwin;
    if rms<=0
      agc_tr(1) = 0;
    else
      agc_tr(1) = 1/sqrt(rms);
    end

    % ramping on.
    for it=(it1(ir)-1):-1:(it1(ir)-iwagc)
      sum = sum + data(it-iwagc)^2;
      nwin = nwin + 1;
      rms = sum/nwin;
      if rms<=0
        agc_tr(it) = 0;
      else
        agc_tr(it) = 1/sqrt(rms);
      end
    end
    
    % middle range -- full rms window.
    for it=(it1(ir)-iwagc-1):-1:(it0(ir)+iwagc)
      sum = sum + data(it-iwagc)^2 - data(it+iwagc)^2;
      rms = sum/nwin;
      if rms<=0
        agc_tr(it) = 0;
      else
        agc_tr(it) = 1/sqrt(rms);
      end
    end
    
    % ramping off.
   for it=(it0(ir)+iwagc-1):-1:it0(ir)
      sum = sum - data(it+iwagc)^2;
      nwin = nwin - 1;
      rms = sum/nwin;
      if rms<=0
        agc_tr(it) = 0;
      else
        agc_tr(it) = 1/sqrt(rms);
      end
    end

    %--------------------------------------------------------------%
    % 
    agc_factor(ir,:) = agc_tr;

  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Apply AGC.
  sh_o_agc = agc_factor.*sh_o;
  
end




