% Solves the incompressible Navier-Stokes equations in a
% rectangular pipe with prescribed velocities along the
% boundary. The solution method is finite differencing on
% a staggered grid with implicit diffusion and a SIMPLE method 
% for the pressure-velocity coupling.
% Visualization is done by a colormap-isoline plot for
% pressure and normalized quiver and streamline plot for
% the velocity field.

% 01/2015 by Hua Zhou
% Take density variation into consideration
%-----------------------------------------------------------------------

clear all;
mu=0.01; % dynamic viscosity /(kg/(m*s))
Re=10;
dt = 1e-3; % time step
tf = 4e-1; % final time
lx = 1; % width of box
ly = 0.1; % height of box
nx = 100; % number of x-gridpoints
ny = 10; % number of y-gridpoints
nsteps = 10; % number of steps with graphic output
%-----------------------------------------------------------------------
nt = ceil(tf/dt); dt = tf/nt;
x = linspace(0,lx,nx+1); hx = lx/nx;
y = linspace(0,ly,ny+1); hy = ly/ny;
% [X,Y] = meshgrid(y,x);
%-----------------------------------------------------------------------
tic;
fprintf('start initialization\n');
% initial conditions
U = zeros(nx-1,ny); V = zeros(nx,ny-1);
P=zeros(nx,ny); pc=zeros(nx*ny,1);
% boundary conditions
uN = x*0; vN = avg(x)*0;
uS = x*0; vS = avg(x)*0;
uW = avg(y)*0+1; vW = y*0;
vE = y*0; uE = avg(y)*0; % this is not the actual bc for uE
%-----------------------------------------------------------------------
Ubc = [2*uS(2:end-1)' zeros(nx-1,ny-2) 2*uN(2:end-1)']/hx^2+...
[uW;zeros(nx-2,ny)]/hy^2;
Vbc = [vS' zeros(nx,ny-3) vN']/hx^2+...
[2*vW(2:end-1);zeros(nx-2,ny-1);2*vE(2:end-1)]/hy^2;

fprintf('Initialization complete, start prepocessing\n');
Lp = kron(eye(ny),Second_full(nx,hx,1,3))+kron(Second_full(ny,hy,1,1),eye (nx));
perp = symamd(Lp); mRp = chol(-1*Lp(perp,perp)); mRpt = mRp';

Lu = kron(eye(ny),Second_full(nx-1,hx,2,1))+...
kron(Second_full(ny,hy,3,3),eye(nx-1));
% peru = symamd(Lu); Ru = chol(Lu(peru,peru)); Rut = Ru';

Lv = kron(eye(ny-1),Second_full(nx,hx,3,3))+...
kron(Second_full(ny-1,hy,2,2),eye(nx));
% perv = symamd(Lv); Rv = chol(Lv(perv,perv)); Rvt = Rv';

fprintf('Prepocessing complete, start iteration\n');
for k = 1:nt
%%%% calculate u and v by substituting p into momentum equations
gamma = min(1.2*dt*max(max(max(abs(U)))/hx,max(max(abs(V)))/hy),1);
Ue = [uW;U;uE]; Ue = [2*uS'-Ue(:,1) Ue 2*uN'-Ue(:,end)];
Ve = [vS' V vN']; Ve = [2*vW-Ve(1,:);Ve;2*vE-Ve(end,:)];
Ua = avg(Ue')'; Ud = diff(Ue')'/2;
Va = avg(Ve); Vd = diff(Ve)/2;
UVx = diff(Ua.*Va-gamma*abs(Ua).*Vd)/hx;
UVy = diff((Ua.*Va-gamma*Ud.*abs(Va))')'/hy;
Ua = avg(Ue(:,2:end-1)); Ud = diff(Ue(:,2:end-1))/2;
Va = avg(Ve(2:end-1,:)')'; Vd = diff(Ve(2:end-1,:)')'/2;
U2x = diff(Ua.^2-gamma*abs(Ua).*Ud)/hx;
V2y = diff((Va.^2-gamma*abs(Va).*Vd)')'/hy;

T=300+500*linspace(0.5*hy,ly-0.5*hy,ny);
% T=300*ones(1,ny); % constant density 
rho_u=300./T;

T=300+500*linspace(0,ly,ny+1);
% T=300*ones(1,ny+1); % constant density 
rho_v=300./T;

Px=reshape(kron(eye(ny),First_back(nx,hx))*reshape(P,[],1),nx-1,ny);
Py=reshape(kron(First_back(ny,hy),eye(nx))*reshape(P,[],1),nx,ny-1);

Dx=mu*(reshape(Lu*reshape(U,[],1),nx-1,ny)+Ubc);
Dy=mu*(reshape(Lv*reshape(V,[],1),nx,ny-1)+Vbc);

rho_tmp=repmat(rho_u,nx-1,1);
U=U+dt*((Dx-Px)./rho_tmp-U2x-UVy(2:end-1,:));

rho_tmp=repmat(rho_v(2:end-1),nx,1);
V=V+dt*((Dy-Py)./rho_tmp-V2y-UVx(:,2:end-1));
uE=U(end,:);

%%%% calculate corrections for u, v, and p
Ue=[uW;U;uE];
rho_tmp=repmat(rho_u,nx+1,1);
RhoUex=diff(rho_tmp.*Ue)/hx;
Ve=[vS', V, vN'];
rho_tmp=repmat(rho_v,nx,1);
RhoVey=diff((rho_tmp.*Ve)')'/hy;

rhs=reshape(RhoUex+RhoVey,[],1)/dt;
% pc=Lp\rhs;
% pc=inv(Lp)*rhs;
pc(perp)=-1*mRp\(mRpt\rhs(perp));

Pc=reshape(pc,nx,ny);

Pcx=reshape(kron(eye(ny),First_back(nx,hx))*pc,nx-1,ny);
Pcy=reshape(kron(First_back(ny,hy),eye(nx))*pc,nx,ny-1);

rho_tmp=repmat(rho_u,nx-1,1);
Uc=-Pcx./rho_tmp*dt;

rho_tmp=repmat(rho_v(2:end-1),nx,1);
Vc=-Pcy./rho_tmp*dt;

P=P+Pc;
U=U+Uc;
V=V+Vc;
uE=U(end,:);

% visualization
if floor(25*k/nt)>floor(25*(k-1)/nt), fprintf('.'), end
if k==1||floor(nsteps*k/nt)>floor(nsteps*(k-1)/nt)
% stream function
%rhs = reshape(diff(U')'/hy-diff(V)/hx,[],1);
%q(perq) = Rq\(Rqt\rhs(perq));
%Q = zeros(nx+1,ny+1);
%Q(2:end-1,2:end-1) = reshape(q,nx-1,ny-1);
% clf, contourf(avg(x),avg(y),P',50,'LineColor','none'), hold on

figure(1)
contourf(avg(x),avg(y),P',50,'w-'), hold on
colorbar('location','southoutside');
% clf, contourf(avg(x),avg(y),P'), hold on
%contour(x,y,Q',20,'k-');
Ue = [uS' avg([uW;U;uE]')' uN'];
Ve = [vW;avg([vS' V vN']);vE];
Len = sqrt(Ue.^2+Ve.^2+eps);
quiver(x,y,(Ue./Len)',(Ve./Len)',.4,'k-')
hold off, axis equal, axis([0 lx 0 ly])
caxis([0 12]);
title(sprintf('Contour of static pressure /(Pa) for Re = %0.1g at t = %0.2g',Re,k*dt),'FontSize',12);
drawnow
filename=sprintf('Contour of static pressure at t = %0.2g.fig', k*dt);
savefig(filename);

figure(2)
y_tmp=linspace(hy/2,ly-hy/2,ny);
y_tmp=[0, y_tmp, ly];
U_tmp=[0, U(end, :),0];
plot(U_tmp, y_tmp,'--ko','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','r');
set(gca, 'FontSize',12);
title(sprintf('Velocity profile at pipe end /(m/s) for Re = %0.1g at t = %0.2g',Re,k*dt),'FontSize',12);
xlabel('Horizontal velocity /(m/s)','FontSize',14);
ylabel('Vertical position /(m)','FontSize',14);
drawnow
filename=sprintf('Velocity profile at pipe end at t = %0.2g.fig', k*dt);
savefig(filename);
end
end
toc;
fprintf('\n')

%%%% Plot figures invariant with time
% plot density 
figure(3)
rho_tmp=repmat(rho_u,nx-1,1);
x_tmp=linspace(hx,lx-hx,nx-1);
y_tmp=linspace(hy/2,ly-hy/2,ny);
contourf(x_tmp,y_tmp,rho_tmp',10,'w-'), hold on
colorbar('location','southoutside');
axis equal, axis([0 lx 0 ly]);
caxis([0.85 1]);
set(gca, 'FontSize',12);
title('Density distribution /(kg/m^3)','FontSize',14);

% plot outlet horizontal velocity given by this code and FLUENT calculation
figure(4)
y_tmp=linspace(hy/2,ly-hy/2,ny);
y_tmp=[0, y_tmp, ly];
U_tmp=[0, U(end, :),0];
y_FLUENT=[1.00E-01,9.00E-02,8.00E-02,7.00E-02,6.00E-02,5.00E-02,4.00E-02,3.00E-02,...
2.00E-02,1.00E-02,0.00E+00];
U_FLUENT=[0.00E+00,5.30E-01,9.42E-01,1.24E+00,1.41E+00,1.47E+00,1.41E+00,1.24E+00,...
9.42E-01,5.30E-01,0.00E+00];
% plot(U_tmp, y_tmp,'--ko','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','r'); hold on;
plot(U_tmp, y_tmp,'--ko','LineWidth',2,'MarkerSize',10); hold on;
plot(U_FLUENT, y_FLUENT,'--rx','LineWidth',2,'MarkerSize',10); hold off;
legend('This code','FLUENT')
set(gca, 'FontSize',12);
title(sprintf('Velocity profile at pipe end /(m/s) for Re = %0.1g at steady state',Re),'FontSize',12);
xlabel('Horizontal velocity /(m/s)','FontSize',16);
ylabel('Vertical position /(m)','FontSize',16);