function [uv,ue,H,xs,ys]=fastearth(Nin,tf,dtyear,Z,type,graphflag)
% FASTEARTH  Use the elastic-viscous earth model from
%     Bueler (2006) "Computation of a viscoelastic deformable earth
%     model for ice sheet simulation."
% to simulate the deformation of the earth under a disc load or under an
% artificial ice sheet load history.  Uses a modified "load response
% matrix" approach for the elastic component.  Uses a Fourier spectral
% method (Trefethen 2000) for the flat-earth, elastic plate over viscous
% half-space model.  The PDE for this latter model is
%     d/dt(2 eta Lap^{1/2} u) + rho g u + D Lap^2 u = sigma_zz.
% Time discretization of this PDE is Crank-Nicolson.  The source of the
% model is
%     Lingle & Clark (1985) "A numerical model of interactions between a
%     marine ice sheet and the solid earth: Application to a West Antarctic
%     ice stream," J. Geophysical Research 90 (C1) 1100--1114.
%
%     [UV,UE,HH,XX,YY]=FASTEARTH(N,TF,DTYEAR,Z)  computes and plots
%     deformation of the earth's surface under an ice disc load of radius
%     1000km and thickness 1000m.  Uses N grid points in each direction on
%     a square 4000 km by 4000 km region R of physical interest.  N must be
%     even; best choices are powers of 2.  Z is a small whole number
%     (Z=2 recommended) which determines the size of the computational
%     region Omega for the flat viscous model.  Namely, Omega is a square
%     with side 4000*Z km, discretized by N*Z points in each direction.  A
%     periodic boundary condition is applied at the edge of Omega.  TF is
%     # of years of the run with timestep DTYEAR.  Elastic deformation is
%     performed on R and is only done if a file "I{N}.mat" can be found;
%     this can be generated using GEFORCONV; see examples.  FASTEARTH
%     returns the gridded flat viscous deformation UV, the gridded elastic
%     deformation UE, the thickness HH--all at the final time---and the 
%     grid itself in XX,YY (as from MESHGRID); all outputs in meters.
%
%     [UV,UE,HH,XX,YY]=FASTEARTH(N,TF,DTYEAR,Z,'disc')  Same as above.
%
%     [UV,UE,HH,XX,YY]=FASTEARTH(N,TF,DTYEAR,Z,'testC')  is the same except
%     that the ice load is essentially test C, but with simple isostacy 
%     added, from
%         Bueler, et al. (2005) "Exact solutions and the verification of
%         numerical models for isothermal ice sheets," J. Glaciol. 51 (173)
%         291--306.
%     Thickness evolution is stopped at 40033 years and the last thickness
%     held.
%
%     [UV,UE,HH,XX,YY]=FASTEARTH(N,TF,DTYEAR,Z,'...',false)  No graphing.
%
% Example I; generate I for 16 x 16 case, store and use:
%     tic, I = geforconv(16,16,2000,2000); toc  % about a minute
%     save 'I16' I
%     fastearth(16,100000,500,2); % 2 secs
% Example II; to see ice cap on bed; needs "I128.mat":
%     [uv,ue,H,x,y]=fastearth(128,40000,500,2,'testC'); % 15 secs
%     figure, mesh(x/1000,y/1000,H+uv+ue)
%     hold on, mesh(x/1000,y/1000,uv+ue), hold off
% See also GEFORCONV, VISCDISC, TESTFAST.
% ELB 1/12/06

if nargin<5, Cflag=false; type='disc'; end
if nargin<6, graphflag=true; end
if strcmp(type,'testC'), Cflag=true;
else, Cflag=false; end
if ((floor(Z)~=ceil(Z)) || (Z<1)), error('Z must be positive integer'), end

g=9.81;  spera=3.1556926e7;  % seconds per year
rhoi=0.910e3;       % density of ice; kg/m^3
rhor=3.300e3;       % density of mantle; kg/m^3
D=5.0e24;           % flexural rigidity of lithosphere as an elastic plate
eta=1.0e21;         % viscosity of mantle

L0=2000e3;          % half-length of actual domain in each direction; m

% compute flat-earth, viscous half-space model on extended grid:
L=Z*2000e3;         % half-length of computational domain
N=Z*Nin;            % computational N
h=2*L/N;   x=-L+h:h:L;
[xx,yy]=meshgrid(x,x);  rr=sqrt(xx.^2+yy.^2);
M=ceil(tf/dtyear);   dt=(tf/M)*spera;

% Fourier coefficients of powers of Laplacian
cx=(pi/L)*[0:N/2 N/2-1:-1:1];
[ccxx,ccyy]=meshgrid(cx,cx);
cclap=ccxx.^2+ccyy.^2;
cchalf=sqrt(cclap);  % lap^{1/2} coeffs
ccbih=cclap.^2;      % biharmonic coeffs

% coefficients in transformed PDE
part1=2*eta*cchalf;
part2=(dt/2)*( rhor*g*ones(size(ccbih)) + D*ccbih );
right = part1 - part2;  % Fourier-transformed operator on right of eqn (19)
left = part1 + part2;   % ... on left of eqn

disp('  computing flat viscous deformation ...')
uun1=zeros(size(xx));
for n=0:M-1
    uun=uun1;
    if Cflag, H=getTestC(xx,yy,(n+1/2)*dt);
    else, H=1000*(rr<1e6); end
    % apply Fourier spectral method:
    sszz=-rhoi*g*H;
    frhs=right.*fft2(uun) + fft2(dt*sszz);
    uun1=real(ifft2( frhs./left ));
end
% tweak: find average value along "distant" bdry of [-ZL,ZL]X[-ZL,ZL],
% remove it and add back result for equivalent disc load
uun1=uun1-( sum(uun1(1,:))+sum(uun1(:,1)) )/(2*N);
if Cflag % assume equiv disc load has R0=1000km
    H0=h*h*sum(sum(H))/(pi*1e6^2);  % trapezoid rule
else, H0=1000; end
uun1=uun1+viscdisc(H0,1000,tf,L/1000);

sh=(N/2)-(N/(2*Z))+1:(N/2)+(N/(2*Z));  % extract central/actual part
xs=xx(sh,sh);   ys=yy(sh,sh);   uv=uun1(sh,sh);  H=H(sh,sh);

% get and use (w. convolution) elastic LRM if available; plot things
filename=['I' num2str(Nin) '.mat'];
doGE=(exist(filename)==2);
if doGE, disp('  [elastic load response matrix FOUND]')
else, disp('  [elastic load response matrix NOT found]'), end
if doGE
    disp(['  computing spherical elastic deformation using I' ...
        num2str(Nin) '.mat ...'])
    S=load(filename);
    ue=rhoi*conv2(H,S.I,'same');
    if graphflag
        figure(1),  clf,  subplot(2,1,2),  mesh(xs/1000,ys/1000,ue);
        xlabel('x (km)'),  ylabel('y (km)')
        zlabel('vertical displacement (m)')
        title('spherical, self-gravitating elastic model')
        subplot(2,1,1),  mesh(xs/1000,ys/1000,uv);
    end
else
    disp('  skipping spherical elastic deformation ...')
    ue=zeros(size(uv)); % for return value only
    if graphflag, figure(1),  clf,  mesh(xs/1000,ys/1000,uv); end
end
if graphflag
    xlabel('x (km)'),  ylabel('y (km)')
    zlabel('vertical displacement (m)')
    title('elastic plate over viscous half-space model')
end

function H=getTestC(x,y,t)
% compute exact H from test C (with simple isostacy added) in Bueler et al 
% "Exact solutions ... isothermal ice sheets".  Stops growth at time 
% t0=40033 years.
spera=3.1556926e7;  % seconds per year
H0=3600;  R0=750000;
% see "Exact solutions ..." eqn (9):
% t0=(2/9.0177e-13)*(7/((1-(910/3300))*4))^3*(R0^4/H0^7) is 40034 years
t0=40034*spera;
ts=t/t0;  if ts>1, ts=1.0; end
rscl=(ts^(-2))*(sqrt(x.^2+y.^2)/R0);
H=H0*ts*max(0, 1-rscl.^(4/3) ).^(3/7);
