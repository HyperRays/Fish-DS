
  function drawslice(name,cycle,seg,range)

% function drawslice(name,cycle,seg,range)
%
% name ... file name
% cycle ... time step number
% seg ... segment (usually plane) in output file to be displayed
% range = [xmin xmax ymin ymax zmin zmax] to be displayed
% example>> cd data;
% example>> drawslice('data',0,1,1e6*[-1 1 -1 1 -1 1]);
%
% characteristic constants:
% c=1 ... speed
% G=1 ... gravity
c = 1;
G = 1;

% read data
  [date,time,dx,u,x0,y0,z0] = readslice(name,cycle,seg);
% load timetable.dat;
% ibounce = find(timetable(:,5)==0);
% time = time-timetable(ibounce,4);
  for i=1:3
    u(i,:,:,:) = u(i,:,:,:)./u(7,:,:,:)*c;
  end
  for i=4:6
    u(i,:,:,:) = u(i,:,:,:)*c^2/dx*sqrt(4*pi/G);
  end
  for i=8:8
    u(i,:,:,:) = u(i,:,:,:)./u(7,:,:,:);
  end

% rotate data
  if seg==1
    xax = 'x [cm]';
    yax = 'y [cm]';
  elseif seg==2
    nu = size(u,1);
    nx = size(u,2);
    nz = size(u,4);
    u = reshape(u,nu,nx,nz);
    tmp = y0;
    y0 = z0;
    z0 = tmp;
    tmp = u(2,:,:,:);
    u(2,:,:,:) = u(3,:,:,:);
    u(3,:,:,:) = tmp;
    tmp = u(5,:,:,:);
    u(5,:,:,:) = u(6,:,:,:);
    u(6,:,:,:) = tmp;
    xax = 'x [cm]';
    yax = 'z [cm]';
    tmp = range(3:4);
    range(3:4) = range(5:6);
    range(5:6) = tmp;
  elseif seg==3
    nu = size(u,1);
    ny = size(u,3);
    nz = size(u,4);
    u = reshape(u,nu,ny,nz);
    tmp = x0;
    x0 = y0;
    y0 = z0;
    z0 = tmp;
    tmp = u(1,:,:,:);
    u(1,:,:,:) = u(2,:,:,:);
    u(2,:,:,:) = u(3,:,:,:);
    u(3,:,:,:) = tmp;
    tmp = u(4,:,:,:);
    u(4,:,:,:) = u(5,:,:,:);
    u(5,:,:,:) = u(6,:,:,:);
    u(6,:,:,:) = tmp;
    xax = 'y [cm]';
    yax = 'z [cm]';
    tmp = range(1:2);
    range(1:2) = range(3:4);
    range(3:4) = range(5:6);
    range(5:6) = tmp;
  end

% plot data
  nx = size(u,2);
  ny = size(u,3);
  nz = size(u,4);
  x = [x0-nx/2-0.5:x0+nx/2-1-0.5]'*dx;
  y = [y0-ny/2-0.5:y0+ny/2-1-0.5]'*dx;
  z = [z0-nz/2-0.5:z0+nz/2-1-0.5]'*dx;
  
% zoom into data
  [ans,imin] = min(abs(x-range(1)*ones(size(x))));
  [ans,imax] = min(abs(x-range(2)*ones(size(x))));
  [ans,jmin] = min(abs(y-range(3)*ones(size(y))));
  [ans,jmax] = min(abs(y-range(4)*ones(size(y))));
  [ans,kmin] = min(abs(z-range(5)*ones(size(z))));
  [ans,kmax] = min(abs(z-range(6)*ones(size(z))));
  nx = imax-imin+1;
  ny = jmax-jmin+1;
  nz = kmax-kmin+1;
  x = x(imin:imax);
  y = y(jmin:jmax);
  z = z(kmin:kmax);
  u = u(:,imin:imax,jmin:jmax,kmin:kmax);
  [X,Y,Z] = ndgrid(x,y,z);
  
  disp('drawing');

  figure(1);
  clf;
  hold on;
  box on;
  C = u(7,:,:,:)*(c/dx)^2/G;
  C = reshape(C,nx,ny,nz);
  surf(X,Y,Z,C);
  shading flat;
  colorbar;
  title(['density, t = ',num2str(time)]);
  xlabel(xax);
  ylabel(yax);
  pbaspect([1 1 1]);
  axis([x(1)-dx,x(nx)+dx,y(1)-dx,y(ny)+dx,z(1)-dx,z(nz)+dx]);
  if nz==1
    view(0,90)
  end

  figure(2);
  clf;
  hold on;
  box on;
  C = u(8,:,:,:);
  C = reshape(C,nx,ny,nz);
  surf(X,Y,Z,C);
  shading flat;
  colorbar;
  title(['specific energy, t = ',num2str(time)]);
  xlabel(xax);
  ylabel(yax);
  pbaspect([1 1 1]);
  axis([x(1)-dx,x(nx)+dx,y(1)-dx,y(ny)+dx,z(1)-dx,z(nz)+dx]);
  if nz==1
    view(0,90)
  end

  figure(3);
  clf;
  hold on;
  box on;
  v = sqrt(u(1,:,:,1).^2 + u(2,:,:,1).^2);
  v = reshape(v,nx,ny);
  pcolor(x,y,v');
  cax = caxis;
  shading flat;
  colorbar;
  title(['velocity, t = ',num2str(time)]);
  xlabel(xax);
  ylabel(yax);
  nv = 32;
  if (nx>nv) && (ny>nv)
    deltax = floor(nx/nv);
    dy = floor(ny/nv);
  else
    deltax = 1;
    dy = 1;
  end
  xv = x([deltax:deltax:nx]);
  xv = xv + 0.5*dx*ones(size(xv));
  yv = y([dy:dy:ny]);
  yv = yv + 0.5*dx*ones(size(yv));
  vx = u(1,[deltax:deltax:nx],[dy:dy:ny],1);
  vy = u(2,[deltax:deltax:nx],[dy:dy:ny],1);
  vx = reshape(vx,floor(nx/deltax),floor(ny/dy));
  vy = reshape(vy,floor(nx/deltax),floor(ny/dy));
  vx = vx./v([deltax:deltax:nx],[dy:dy:ny]);
  vy = vy./v([deltax:deltax:nx],[dy:dy:ny]);
  h = quiver(xv,yv,vx',vy',0.5,'w');
  pbaspect([1 1 1]);
  axis([x(1)-dx,x(nx)+dx,y(1)-dx,y(ny)+dx]);

  figure(4);
  clf;
  hold on;
  box on;
  map = [winter(32);[0:1/31:1]',ones(32,1),0.5*[1:-1/31:0]'];
  colormap(map);
  v = sqrt(u(4,:,:,1).^2 + u(5,:,:,1).^2);
  v = reshape(v,nx,ny);
  vmin = 10;
  vtmp = log10(max(v,vmin*ones(size(v))));
  pcolor(x,y,vtmp');
  caxis([vmin max(max(max(vtmp)),vmin+1)]);
  shading flat;
  colorbar;
  title(['log10(magnetic field), t = ',num2str(time)]);
  xlabel(xax);
  ylabel(yax);
  nv = 32;
  if (nx>nv) && (ny>nv)
    deltax = floor(nx/nv);
    dy = floor(ny/nv);
  else
    deltax = 1;
    dy = 1;
  end

% stream lines
  pbaspect([1 1 1]);
  axis([x(1)-dx,x(nx)+dx,y(1)-dx,y(ny)+dx]);
  vx = u(4,:,:,1);
  vy = u(5,:,:,1);
  vx = reshape(vx,nx,ny);
  vy = reshape(vy,nx,ny);
  [xx,yy] = meshgrid(x,y);
  str = streamslice(xx,yy,vx',vy');
  set(str,'Color','r');

% vector arrows
% ny = floor(ny/2);
% vx = u(7,[deltax:deltax:nx],[dy:dy:ny],1);
% vy = u(8,[deltax:deltax:nx],[dy:dy:ny],1);
% vx = reshape(vx,floor(nx/deltax),floor(ny/dy));
% vy = reshape(vy,floor(nx/deltax),floor(ny/dy));
% vx = vx./v([deltax:deltax:nx],[dy:dy:ny]);
% vy = vy./v([deltax:deltax:nx],[dy:dy:ny]);
% yv = y([dy:dy:ny]);
% yv = yv + 0.5*dx*ones(size(yv));
% h = quiver(xv,yv,vx',vy',0.5,'w');

  disp('done');

