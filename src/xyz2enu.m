function enu = xyz2enu(xyz,orgxyz)
%XYZ2ENU	Convert from WGS-84 ECEF cartesian coordinates to 
%               rectangular local-level-tangent ('East'-'North'-Up)
%               coordinates.
%
%	enu = XYZ2ENU(xyz,orgxyz)	
%
%    INPUTS
%	xyz(1) = ECEF x-coordinate in meters
%	xyz(2) = ECEF y-coordinate in meters
%	xyz(3) = ECEF z-coordinate in meters
%
%	orgxyz(1) = ECEF x-coordinate of local origin in meters
%	orgxyz(2) = ECEF y-coordinate of local origin in meters
%	orgxyz(3) = ECEF z-coordinate of local origin in meters
%
%    OUTPUTS
%       enu:  Column vector
%		enu(1,1) = 'East'-coordinate relative to local origin (meters)
%		enu(2,1) = 'North'-coordinate relative to local origin (meters)
%		enu(3,1) = Up-coordinate relative to local origin (meters)

%	Reference: Alfred Leick, GPS Satellite Surveying, 2nd ed.,
%	           Wiley-Interscience, John Wiley & Sons, 
%	           New York, 1995.
%
%	M. & S. Braasch 10-96
%	Copyright (c) 1996 by GPSoft
%	All Rights Reserved.

% 입력 검증 및 차벡터 계산 (xyz와 orgxyz 모두 필요)
if nargin<2, error('insufficient number of input arguments'), end

% xyz와 orgxyz의 벡터 크기가 다르면 orgxyz를 전치(transpose)하여 맞춤
tmpxyz = xyz;
tmporg = orgxyz;
if size(tmpxyz) ~= size(tmporg), tmporg = tmporg'; end

% 변환할 점과 기준점 간의 차벡터 (ECEF 기준 상대 위치)
difxyz = tmpxyz - tmporg;

% 차벡터가 행벡터면 열벡터로 변환
[m,n] = size(difxyz); 
if m<n, difxyz = difxyz'; end

% 기준점의 ECEF 좌표를 위도(latitude), 경도(longitude), 고도(height)로 변환
orgllh = xyz2llh(orgxyz);

phi = orgllh(1);   % 기준점 위도 (latitude)
lam = orgllh(2);   % 기준점 경도 (longitude)

sinphi = sin(phi); cosphi = cos(phi);
sinlam = sin(lam); coslam = cos(lam);

% ECEF → ENU 변환을 위한 회전 행렬 R 구성
R = [ -sinlam          coslam         0     ; ...
      -sinphi*coslam  -sinphi*sinlam  cosphi; ...
       cosphi*coslam   cosphi*sinlam  sinphi];

% 회전행렬 R을 곱해 ENU(동,북,천) 좌표로 변환
enu = R * difxyz;
