%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Elevation & Azimuth angle
%
%..........................................................................
%
% IN:   SAT= Position of SV in ECEF [m]
%       U=  Position of User in ECEF [m]
%       constanst= constant parameters
%
% OUT:  ELEV= Elevation Angle in fixed user point [deg]
%       AZ= Azimuth Angle in fixed user point [deg] 
%
%%%%%%%%%%i%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ELEV,AZ]=eleazi(SAT,U)

% D2R= constants.d2r;
% R2D= constants.r2d;
D2R=pi/180;
R2D=180/pi;

%L = xyz2LLH(U);         % WGS-84 coordinate (deg,deg,m)

%L(1:2)=L(1:2)*D2R;      %
%LSX=-sin(L(2))*SAT(1)+cos(L(2))*SAT(2)...
%    +U(1)*sin(L(2))-U(2)*cos(L(2));
%LSY=-sin(L(1))*cos(L(2))*SAT(1)-sin(L(1))*sin(L(2))*SAT(2)+cos(L(1))*SAT(3)...
%    +U(1)*sin(L(1))*cos(L(2))+U(2)*sin(L(1))*sin(L(2))-U(3)*cos(L(1));
%LSZ=cos(L(1))*cos(L(2))*SAT(1)+cos(L(1))*sin(L(2))*SAT(2)+sin(L(1))*SAT(3)...
%    -U(1)*cos(L(1))*cos(L(2))-U(2)*cos(L(1))*sin(L(2))-U(3)*sin(L(1));
%LSV = [LSX LSY LSZ];

LLH = xyz2llh(U);         % WGS-84 coordinate (rad,rad,m)

R= rot2(LLH(1,1),LLH(1,2));
rho_xyz= SAT-U(1:3);
LSV= R*rho_xyz'; %Local level position of SV
LSX= LSV(1); LSY= LSV(2); LSZ= LSV(3);

%ELEV = atan(LSZ/sqrt(LSX^2+LSY^2));
ELEV = asin(LSZ/sqrt(LSX^2+LSY^2+LSZ^2));
ELEV = ELEV*R2D;
AZ=atan2(LSX,LSY);
AZ=AZ*R2D;

if AZ < 0
   AZ= 180+(180-abs(AZ));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transformation Matrix (into NEU coordinate)
% ( r(spherical) <- R * r(xyz) )
% (x,y,z) -> (lambda,phi,r)
%
% by Kai Borre
%...........................................................................
% IN:  phi= latitude [rad]
%      ramda= longitude [rad]
%
% OUT:  rotm= rotation matirx [3x3]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function rotm=rot2(phi,ramda)

      cp= cos(phi);
      sp= sin(phi);
      cl= cos(ramda);
      sl= sin(ramda);
      rotm=[-sl cl 0 ; -sp*cl -sp*sl cp; cp*cl cp*sl sp];