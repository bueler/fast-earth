function II=geforconv(Nx,Ny,Lx,Ly);
% GEFORCONV  Computes matrix I=I(p,q), a form of the elastic "load response
% matrix" described in Lingle & Clark (1985) and in the technical report
% Bueler (2005).   Matrix I(p,q) is used by FASTEARTH.
%
%     I=GEFORCONV(NX,NY,LX,LY)  computes (2 NX + 1) by (2 NY + 1)
%     matrix I with entries
%             /dy/2   /dx/2
%     I(p,q)= |       |      G^E(sqrt((pdx-xi)^2+(qdy+eta^2))) dxi deta
%             /-dy/2  /-dx/2
%     where dx=2 LX/NX, dy=2 LY/NY.  Input LX,LY in km.  Assumes region 
%     extends from x=-LX to x=LX and y=-LY to y-LY.  Uses Matlab's 
%     dblquad to do integral.
%
% Example:
%     I=geforconv(16,16,2000,2000);
%
% See also FASTEARTH.
% ELB 12/11/05.

rm=[ 0.011  0.111  1.112  2.224  3.336  4.448  6.672  8.896  11.12 17.79 ...
     22.24  27.80  33.36  44.48  55.60  66.72  88.96  111.2  133.4 177.9 ...
     222.4  278.0  333.6  444.8  556.0  667.2  778.4  889.6 1001.0 1112.0 ...
    1334.0 1779.0 2224.0 2780.0 3336.0 4448.0 5560.0 6672.0 7784.0 8896.0 ...
   10008.0] * 1e3;  % converted to meters
% GE /(10^12 rm) is vertical displacement in meters
GE=[-33.64 -33.56 -32.75 -31.86 -30.98 -30.12 -28.44 -26.87 -25.41 ...
    -21.80 -20.02 -18.36 -17.18 -15.71 -14.91 -14.41 -13.69 -13.01 ...
    -12.31 -10.95 -9.757 -8.519 -7.533 -6.131 -5.237 -4.660 -4.272 ...
    -3.999 -3.798 -3.640 -3.392 -2.999 -2.619 -2.103 -1.530 -0.292 ...
     0.848  1.676  2.083  2.057  1.643];  
% linearly extrapolate GE to r=0;  GE(0) := -33.6488
GE=[interp1(rm(1:2),GE(1:2),0.0,'linear','extrap') GE];
rm=[0.0 rm];  % length(rm)=length(GE)=42

dx=2*Lx*1000/Nx;  dy=2*Ly*1000/Ny;  % dimensions of load element
II=zeros(2*Nx-1,2*Ny-1);
% compute entries of II by dblquad, using quad method and default TOL=1e-6 
for p=-Nx+1:Nx-1
    for q=-Ny+1:Ny-1
        II(p+Nx,q+Ny)=dblquad(@integrand,-dx/2,dx/2,-dy/2,dy/2);
    end
end

% nested function which is integrand in I(p,q)
    function z=integrand(xi,eta)
        r=sqrt((p*dx-xi).^2+(q*dy-eta).^2);
        z=zeros(length(r));
        for jj=1:length(r)
            if r(jj)==0.0  % treat normalization r as 11 m
                z(jj)=GE(1)/(rm(2)*1e12);
            elseif r(jj)>=rm(end), z(jj)=0;
            else  % linearly interpolate to get GE value; then normalize
                ii=find(rm>r(jj),1);  %  2 <= ii <= 42
                z(jj)=interp1(rm(ii-1:ii),GE(ii-1:ii),r(jj))/(r(jj)*1e12);
            end
        end
    end

end
