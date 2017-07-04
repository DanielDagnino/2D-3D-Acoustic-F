%{
!*********************************************************************/
!** This code has been done in the Barcelona Center for Subsurface 
!** Imaging (BCSI).
!** Goal: Set of tools to analyse the FWI results.
!** Authors: Daniel Dagnino.
!*********************************************************************/
%}

function [ data, n_rec,nt, dt, check, time, data_min,data_max ] = read_data( file )
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 
  fid = fopen(file,'rb');
  
  %--------------------------------------------------------------%
  % 
  fread( fid, 1, 'int32' );
  ndim = fread( fid, 2, 'int32' );
  fread( fid, 1, 'int32' );
  n_rec = ndim(1);
  nt    = ndim(2);
%   [n_rec, nt]
  
  %--------------------------------------------------------------%
  % 
  fread( fid, 1, 'int32' );
  check = fread( fid, 1, 'int32' );
  fread( fid, 1, 'int32' );
  
  %--------------------------------------------------------------%
  % 
  fread( fid, 1, 'int32' );
  dt = double(fread( fid, 1, 'float32' ));
  fread( fid, 1, 'int32' );
  
  %--------------------------------------------------------------%
  % 
  data = zeros(n_rec,nt);
  for it=1:nt
    fread( fid, 1, 'int32' );
    data(:,it) = double(fread( fid, n_rec, 'float32' ));
    fread( fid, 1, 'int32' );
  end
  
  %--------------------------------------------------------------%
  % 
  fclose(fid);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 
  if nargout > 5
    time = dt*(0:nt-1);
  end
  
  % 
  if nargout > 6
    data_min = min(min(data));
    data_max = max(max(data));
  end
  
end

