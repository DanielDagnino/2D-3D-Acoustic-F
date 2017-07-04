%{
!*********************************************************************/
!** This code has been done in the Barcelona Center for Subsurface 
!** Imaging (BCSI).
!** Goal: Set of tools to analyse the FWI results.
!** Authors: Daniel Dagnino.
!*********************************************************************/
%}

function [ parameter, nx,ny,nz, dx, x,y,z, par_min,par_max ] = read_model( file )
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 
  fid = fopen(file,'rb');
  
  %--------------------------------------------------------------%
  % 
  fread( fid, 1, 'int32' );
  ndim = fread( fid, 3, 'int32' );
  fread( fid, 1, 'int32' );
  nx = ndim(1);
  ny = ndim(2);
  nz = ndim(3);
  
  %--------------------------------------------------------------%
  % 
  fread( fid, 1, 'int32' );
  dx = double(fread( fid, 1, 'float32' ));
  fread( fid, 1, 'int32' );
  
  %--------------------------------------------------------------%
  % 
  parameter = zeros(ny,nx);
  for iy=1:ny
    fread( fid, 1, 'int32' );
    parameter(iy,:) = double(fread( fid, nx, 'float32' ));
    fread( fid, 1, 'int32' );
  end
  
  %--------------------------------------------------------------%
  % 
  fclose(fid);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 
  if nargout > 5
    x = dx*(0:nx-1);
    y = dx*(0:ny-1);
    z = dx*(0:nz-1);
  end
  
  % 
  if nargout > 8
    par_min = min(min(parameter));
    par_max = max(max(parameter));
  end
  
end

