%{
!*************************************************************************/
!** This code has been done in the Barcelona Center for Subsurface 
!** Imaging (BCSI).
!** Goal: Set of tools to analyse the FWI results.
!** Authors: Daniel Dagnino.
!*************************************************************************/
%}

function [ sg_f ] = filt_butterworth_all( sg, dt, f_type, freq_min, freq_max, f_order, norm )
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Size.
  [nt, n_rec] = size(sg);
  
  % Filter.
  filt = filt_butterworth( nt, dt, f_type, freq_min, freq_max, f_order );
  
  % Filtering.
  sg_f = zeros(nt,n_rec);
  for ir=1:n_rec
    
    % Apply filter.
    sg_f(:,ir) = ifft(filt'.*fft(sg(:,ir)));
    
    % Normalization.
    if norm
      sg_f(:,ir) = sg_f(:,ir)/max(max(abs(sg_f(:,ir))));
    end
    
  end
  
end




