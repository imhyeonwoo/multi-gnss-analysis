function llh = xyz2llh(xyz)
%XYZ2LLH	Convert from ECEF cartesian coordinates to 
%               latitude, longitude and height.  WGS-84
%
%	llh = XYZ2LLH(xyz)	
%
%    INPUTS
%	xyz(1) = ECEF x-coordinate in meters
%	xyz(2) = ECEF y-coordinate in meters
%	xyz(3) = ECEF z-coordinate in meters
%
%    OUTPUTS
%	llh(1) = latitude in radians
%	llh(2) = longitude in radians
%	llh(3) = height above ellipsoid in meters

%	Reference: Understanding GPS: Principles and Applications,
%	           Elliott D. Kaplan, Editor, Artech House Publishers,
%	           Boston, 1996.
%
%	M. & S. Braasch 10-96
%	Copyright (c) 1996 by GPSoft
%	All Rights Reserved.
%
	x = xyz(1);
	y = xyz(2);
	z = xyz(3);
	x2 = x^2;
	y2 = y^2;
	z2 = z^2;

	a = 6378137.0000;	    % earth radius in meters
	b = 6356752.3142;       % earth semiminor in meters	
	e = sqrt (1-(b/a).^2);  % first eccentricity (1차 이심률) -> 지구가 얼마나 납작한지 나타내는 척도
	b2 = b*b;
	e2 = e^2;
	ep = e*(a/b);       % second eccentricity (보조 이심률) -> 위도 계산식에서 z축 성분 보정용
	r = sqrt(x2+y2);    % equatorial radius component (적도면 반경 성분) -> 적도면에서 지구 자전축까지의 거리
	r2 = r*r;
	E2 = a^2 - b^2;

    % 폐형식 변환 공식의 중간 계산 (Kaplan, 1996)
	F = 54*b2*z2;
	G = r2 + (1-e2)*z2 - e2*E2;
	c = (e2*e2*F*r2)/(G*G*G);
	s = ( 1 + c + sqrt(c*c + 2*c) )^(1/3);
	P = F / (3 * (s+1/s+1)^2 * G*G);
	Q = sqrt(1+2*e2*e2*P);
	ro = -(P*e2*r)/(1+Q) + sqrt((a*a/2)*(1+1/Q) ...
                                - (P*(1-e2)*z2)/(Q*(1+Q)) - P*r2/2);

	tmp = (r - e2*ro)^2;
	U = sqrt( tmp + z2 );           % distance from Earth center (지구 중심 거리)
	V = sqrt( tmp + (1-e2)*z2 );    % adjusted distance (이심률 보정 거리)
	zo = (b2*z)/(a*V);              % corrected z-component (보정된 z성분) -> 위도 계산식에 사용

	height = U*( 1 - b2/(a*V) );    % height above ellipsoid (Geodetric or Ellipsoidal Height)
	
	lat = atan( (z + ep*ep*zo)/r ); % geodetic latitude

    % longitude
	temp = atan(y/x);
	if x >=0	
		long = temp;
	elseif (x < 0) & (y >= 0)
		long = pi + temp;
	else
		long = temp - pi;
	end

	llh(1) = lat;
	llh(2) = long;
	llh(3) = height;

