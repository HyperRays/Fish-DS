
  function draw_basic(name,cycle,seg,range)

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

figure(3);
clf;
hold on;
box on;

for i=1:cycle

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

disp('done');

end



