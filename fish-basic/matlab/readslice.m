
  function [date,time,dx,u,x0,y0,z0] = readslice(name,cycle,seg);

% function [date,time,dx,u,x0,y0,z0] = readslice(name,cycle,seg);
% date ... date of output
% time ... time of slice
% dx ... physical length of one zone width
% name ... model name
% cycle ... time step number
% seg ... segment number
% u(1:nu,x,y,z) ... array with all variables 1:nu at index x0+x, y0+y, z0+z
% x0, y0, z0 ... index offset to save memory in MATLAB

% convert time step number to string
  zstr = '00000';
  cycstr = int2str(cycle);
  cycstr = [zstr(1:5-size(cycstr,2)),cycstr];

%debug
  tmp = 0;

% check if file exists and open
  ip = 1; % number of file
  ia = 1; % number of data chunk
  filename = ['proc',zstr(1:4),'/',name,'_',cycstr,'_',zstr(1:4),'.dat'];
  while exist(filename,'file')==2
    disp(['reading ',filename]);
    fid = fopen(filename);

% read general information
    [dummy,count] = fread(fid,4,'uchar');
    [date,count] = fread(fid,8,'char');
    date = char(date');
    [dummy,count] = fread(fid,8,'uchar');
    [time,count] = fread(fid,1,'float');
    [dx,count] = fread(fid,1,'float');
    [dummy,count] = fread(fid,8,'uchar');
    [ns,count] = fread(fid,1,'int');
    [num,count] = fread(fid,1,'int');
    [mr,count] = fread(fid,1,'int');

% read segments and store in separate data chunks
    for is=1:min(seg,ns)
      [dummy,count] = fread(fid,8,'uchar');
      [nu(ia),count] = fread(fid,1,'int');
      [dummy,count] = fread(fid,8,'uchar');
      [n(1:6,ia),count] = fread(fid,6,'int');
      nx = max(n(2,ia)-n(1,ia)+1,0);
      ny = max(n(4,ia)-n(3,ia)+1,0);
      nz = max(n(6,ia)-n(5,ia)+1,0);
      [dummy,count] = fread(fid,8,'uchar');
      if nx*ny*nz>0
        if is<seg % skip data
          [dummy,count] = fread(fid,nu(ia)*nx*ny*nz,'float');
        else      % store data
          [tmp(1:nu(ia)*nx*ny*nz,ia),count] = fread(fid,nu(ia)*nx*ny*nz,'float');
        end
      end
    end

% read domain overview and store in separate data chunks
    if seg>ns
      [dummy,count] = fread(fid,8,'uchar');
      [nred,count] = fread(fid,1,'int');
      [i0,count] = fread(fid,1,'int');
      [j0,count] = fread(fid,1,'int');
      [k0,count] = fread(fid,1,'int');
      [dummy,count] = fread(fid,8,'uchar');
      [nu(ia),count] = fread(fid,1,'int');
      [dummy,count] = fread(fid,8,'uchar');
      [n(1:6,ia),count] = fread(fid,6,'int');
      nx = max(n(2,ia)-n(1,ia)+1,0);
      ny = max(n(4,ia)-n(3,ia)+1,0);
      nz = max(n(6,ia)-n(5,ia)+1,0);
      [dummy,count] = fread(fid,8,'uchar');
      if nx*ny*nz>0
        [tmp(1:nu(ia)*nx*ny*nz,ia),count] = fread(fid,nu(ia)*nx*ny*nz,'float');
      end
    end

% close file and prepare opening of next one
    fclose(fid);
    ip = ip+1;
    ipstr = int2str(ip-1);
    ipstr = [zstr(1:4-size(ipstr,2)),ipstr];
    filename = ['proc',ipstr,'/',name,'_',cycstr,'_',ipstr,'.dat'];
    na = ia;
    if nx*ny*nz>0
      ia = ia+1;
    end
  end

% subract unused index space to save memory in MATLAB
  x0 = min(n(1,:))-1;
  y0 = min(n(3,:))-1;
  z0 = min(n(5,:))-1;
  n(1,:) = n(1,:) - x0*ones(1,na);
  n(2,:) = n(2,:) - x0*ones(1,na);
  n(3,:) = n(3,:) - y0*ones(1,na);
  n(4,:) = n(4,:) - y0*ones(1,na);
  n(5,:) = n(5,:) - z0*ones(1,na);
  n(6,:) = n(6,:) - z0*ones(1,na);

% fill data chunks into one output array u(1:nu,ix,iy,iz)
  u = [];
  for ia=1:na
    nx = max(n(2,ia)-n(1,ia)+1,0);
    ny = max(n(4,ia)-n(3,ia)+1,0);
    nz = max(n(6,ia)-n(5,ia)+1,0);
    if nx*ny*nz>0
      u(:,n(1,ia):n(2,ia),n(3,ia):n(4,ia),n(5,ia):n(6,ia))...
        = reshape(tmp(1:nu(ia)*nx*ny*nz,ia),nu(ia),nx,ny,nz);
    end
  end

% tranform to larger dx in domain overview
  if seg>ns
    dx = dx*nred;
    u(10:11,:,:,:) = u(10:11,:,:,:)*nred^2;
  end
