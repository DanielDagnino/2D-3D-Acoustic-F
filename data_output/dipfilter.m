%{
!*********************************************************************/
!** This code has been done in the Barcelona Center for Subsurface 
!** Imaging (BCSI).
!** Goal: Set of tools to analyse the FWI results.
!** Authors: Daniel Dagnino.
!*********************************************************************/
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% My Matlab script for the dipfilter of the SU.
function sg_filt = dipfilter( sg, dt, nt, slopes, amps, bias, dx )
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !rm -f *.su
  !rm -f *.segy
  !rm -f *.sgy
  !rm -f header
  !rm -f binary

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Write SGY.
  WriteSegy('sg_aux.sgy',sg','dt',dt,'ns',nt,'revision',1);

  % From SGY to SU.
  !segyread tape=sg_aux.sgy verbose=0 endian=0 | segyclean > sg_aux.su

  % !suximage perc=98 < sg_aux.su &
  % !sugethw key=dt < sg_aux.su

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Dip filter.
  command = [ 'sudipfilt dx=', num2str(dx) ];
  
  n_slopes = length(slopes);
  
  aux = [];
  for k=1:(n_slopes-1)
    aux = [ aux, num2str(slopes(k)),',' ];
  end
  aux = [ aux, num2str(slopes(n_slopes)) ];
  command = [ command, ' slopes=', num2str(aux) ];
  
  aux = [];
  for k=1:(n_slopes-1)
    aux = [ aux, num2str(amps(k)),',' ];
  end
  aux = [ aux, num2str(amps(n_slopes)) ];
  command = [ command, ' amps=', num2str(aux) ];
  
  command = [ command, ' bias=', num2str(bias) ];
  
  command = [ command, ' < sg_aux.su > sg_filt.su' ];
  [status,cmdout] = system(command);
  
  if status~=0
    disp(cmdout);
  end

  % % Plot depth-migrated model.
  % !suximage perc=98 < sg_filt.su &

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % SU to SEGY.
  !segyhdrs < sg_filt.su
  !segywrite < sg_filt.su tape=sg_filt.sgy

  % Read SEGY.
  [ sg_filt, ~, ~ ] = ReadSegy('sg_filt.sgy');
  sg_filt = sg_filt';
  
end




