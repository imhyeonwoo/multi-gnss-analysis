# Function Anaylsis of xyz2llh and xyz2enu

## 1. xyz2llh.m
This function converts coordinates from the ECEF (Earth-Centered, Earth-Fixed) coordinate system 
[𝑥,𝑦,𝑧] to the LLH (Latitude, Longitude, Height) coordinate system 
[$\phi$,$\lambda$,$h$] based on the WGS-84 reference ellipsoid. The units are [rad,rad,meter].

Here, the converted latitude $\phi$ and altitude $h$ are defined with respect to the geodetic reference frame (the ellipsoid’s surface normal), not the geocentric reference frame (the Earth’s center). In other words, the LLH coordinate system is based on the ellipsoid normal rather than the Earth’s center.

There are two main approaches to performing the ECEF ↔ LLH conversion:

- Iterative approach, such as using the Newton-Raphson method

- Closed-form approach, which does not require iteration

In this report, the `xyz2llh.m` function based on the closed-form method proposed by Kaplan (1996) was used to directly compute latitude, longitude, and altitude. This method has the advantages of being independent of initial values and having low computational cost. However, it also has disadvantages, such as limited extensibility to other ellipsoid models (e.g., GRS80) and potentially larger errors near the polar regions.

---

#### 0) Pre-Calculation
```matlab
a = 6378137.0000;	    % earth radius in meters
b = 6356752.3142;         % earth semiminor in meters	
e = sqrt (1-(b/a).^2);    % first eccentricity
b2 = b*b;
e2 = e^2;
ep = e*(a/b);       % second eccentricity
r = sqrt(x2+y2);    % equatorial radius component
r2 = r*r;
E2 = a^2 - b^2;

% Closed-Form Calculation Preprocess (Kaplan, 1996)
F = 54*b2*z2;
G = r2 + (1-e2)*z2 - e2*E2;
c = (e2*e2*F*r2)/(G*G*G);
s = ( 1 + c + sqrt(c*c + 2*c) )^(1/3);
P = F / (3 * (s+1/s+1)^2 * G*G);
Q = sqrt(1+2*e2*e2*P);
ro = -(P*e2*r)/(1+Q) + sqrt((a*a/2)*(1+1/Q) ...
                               - (P*(1-e2)*z2)/(Q*(1+Q)) - P*r2/2);

tmp = (r - e2*ro)^2;
U = sqrt( tmp + z2 );           % distance from Earth center
V = sqrt( tmp + (1-e2)*z2 );    % adjusted distance
zo = (b2*z)/(a*V);              % corrected z-component 
```

---

#### 1) Longitude $\lambda$ Calculation
Longitude is calculated using x, y components of ECEF coordinates.
$$\lambda = tan^-1(y/x)$$

To ensure that $\lambda \in [-\pi, \pi]$, use the following conditional statements based on the signs of $x$ and $y$:
```matlab
% longitude
temp = atan(y/x);
if x >=0	
	long = temp;
elseif (x < 0) & (y >= 0)
	long = pi + temp;
else
	long = temp - pi;
end
```
Longitude does not require considering the Earth's flattening or eccentricity,  
as it is simply calculated as an angle on the $xy$-plane.  
Moreover, the formula for longitude can be simplified in code by using the `atan2` function.

---

#### 2) Latitude $\phi$ Calculation
Since the Earth is not a perfect sphere but a flattened ellipsoid,  
latitude is calculated using the second eccentricity ($e'$) and the corrected $z$-component ($z_0$) as follows:

$$
\phi = \tan^{-1}\left(\frac{z + {e'}^2 z_0}{r}\right), \qquad 
z_0 = \frac{b^2 z}{aV}
$$

where

$
r = \sqrt{x^2 + y^2}, \qquad 
{e'}^2 = \frac{a^2 - b^2}{b^2} \; (\text{second eccentricity}), 
\qquad V = \text{adjusted distance (eccentricity correction)}.
$

```matlab
% geodetic latitude
lat = atan( (z + ep*ep*zo)/r );
```
$e$ is the first eccentricity, representing the Earth's flattening.  
$e'$ is the second (or auxiliary) eccentricity, used to correct the $z$-component when calculating latitude.  
$z_0$ is the corrected ECEF $z$-component adjusted to the curvature of the reference ellipsoid.  
$r$ is the horizontal distance computed from the ECEF $x$ and $y$ components.  

The latitude obtained through this process is the **geodetic latitude**.

---

#### 3) Height $h$ Calculation
Height $h$ is defined as the distance from the Earth's center minus the distance to the ellipsoid surface.  
Since the Earth is not a perfect sphere but a flattened ellipsoid, it is calculated using the semi-major axis ($a$) and semi-minor axis ($b$) of the ellipsoid as follows.  
For reference, in WGS-84: `a = 6378137.0000`, `b = 6356752.3142`.

$$
h = U \left(1 - \frac{b^2}{aV}\right)
$$

From `xyz2llh.m` [lines 57:58]:
```matlab
% height above ellipsoid (Geodetic or Ellipsoidal Height)
height = U*( 1 - b2/(a*V) );
```

---

## 2. xyz2enu.m
This function converts coordinates from the ECEF coordinate system \([x, y, z]\) based on WGS-84  
to the ENU coordinate system \([e, n, u]\) centered at a specified reference point (origin).

The ECEF (Earth-Centered, Earth-Fixed) coordinate system is a global coordinate system  
with its origin at the Earth's center.  
The ENU (East-North-Up) coordinate system is a local coordinate system  
with its origin at a specific point on the Earth's surface (e.g., user or receiver position).  
In the ENU system, the axes are defined as East, North, and Up from the perspective of the reference point.

This transformation is mainly used to express global positions measured by GPS/GNSS (in ECEF)  
as local positions relative to the user (in ENU).

The ECEF↔ENU conversion consists of the following main steps:  
1) Compute the relative vector from the reference point  
2) Obtain the latitude and longitude of the reference point  
3) Construct a rotation matrix using the latitude and longitude  
4) Apply the rotation matrix to obtain the ENU coordinates  

The provided function follows this procedure and internally calls the `xyz2llh` function  
to obtain the latitude and longitude of the reference point.

---

#### 1) Compute the relative vector from the reference point

The relative position vector $\Delta \mathbf{r}$ is obtained by subtracting the ECEF coordinates  
of the reference point (origin) $\mathbf{r_0}$ from the ECEF coordinates of the target point $\mathbf{r}$:

$$
\Delta \mathbf{r} = \mathbf{r} - \mathbf{r_0} 
= 
\begin{bmatrix}
x - x_0 \\[4pt]
y - y_0 \\[4pt]
z - z_0
\end{bmatrix}
$$

This vector represents the position of the target point relative to the reference point in the ECEF frame.  
It is used as the input for the ENU conversion, and will be projected by the rotation matrix  
to obtain the East, North, and Up components.

Here, $[x,y,z]$ are the ECEF coordinates to be converted,  
and $[x_0, y_0, z_0]$ are the ECEF coordinates of the origin (reference point).

---

**From `xyz2enu.m` [lines 31:44]:**
```matlab
% Validate inputs and compute difference vector (both xyz and orgxyz are required)
if nargin<2, error('insufficient number of input arguments'), end

% If xyz and orgxyz have different vector shapes, transpose orgxyz to match
tmpxyz = xyz;
tmporg = orgxyz;
if size(tmpxyz) ~= size(tmporg), tmporg = tmporg'; end

% Difference vector between target point and reference point (relative position in ECEF)
difxyz = tmpxyz - tmporg;

% If the difference vector is a row vector, convert it to a column vector
[m,n] = size(difxyz); 
if m<n, difxyz = difxyz'; end
```
As seen in the code, the function first validates the inputs and applies transposition if needed.
The orgxyz (a 3-element vector representing the local ENU origin in ECEF coordinates)
is a required input.
The function also automatically corrects row/column mismatches through transposition,
allowing various input formats.

---

#### 2) Obtain the latitude and longitude of the reference point

In the second step, to rotate the relative vector (`difxyz`) into the ENU frame,  
the geographic coordinates (geodetic latitude $\phi$ and longitude $\lambda$)  
of the reference point are required.  

This is because the ENU coordinate axes are always defined on the local tangent plane  
at the reference point.  
In other words, the orientation of the ENU axes depends on **where** the reference point is located  
(i.e., its latitude and longitude).  
Therefore, it is necessary to first convert the reference point’s ECEF coordinates (`orgxyz`)  
into LLH coordinates.

In the provided `xyz2enu.m` function, this is done by calling `xyz2llh(orgxyz)`  
to obtain the LLH coordinates of the reference point (origin):

**From `xyz2enu.m` [lines 46–50]:**
```matlab
% Convert the reference point's ECEF coordinates to latitude, longitude, and height
orgllh = xyz2llh(orgxyz);

phi = orgllh(1);   % latitude of the reference point
lam = orgllh(2);   % longitude of the reference point
```
Here, the latitude ($\phi$) and longitude ($\lambda$) of the reference point
are obtained from its LLH coordinates.
Because the latitude and longitude returned by xyz2llh are in radians,
they are immediately used to compute their sine and cosine values
for convenience in the ENU rotation calculation:
**From `xyz2enu.m` [lines 52-53]:**
```matlab
sinphi = sin(phi); cosphi = cos(phi);
sinlam = sin(lam); coslam = cos(lam);
```

---

#### 3) Construct the rotation matrix using latitude and longitude

The rotation matrix $R(\phi,\lambda)$ expresses the local axes  
(East, North, Up) defined at the reference point’s geographic coordinates  
(geodetic latitude $\phi$, longitude $\lambda$) in terms of the ECEF coordinate components.  

Multiplying this matrix by the relative vector $\Delta \mathbf{r}$
(from Step 1) gives the ENU components of that vector.

The local unit vectors at the reference point, written in ECEF components, are as follows  
(each depends on $\phi$ and $\lambda$):

- **East unit vector**
$$
\mathbf{e} =
\begin{bmatrix}
-\sin\lambda \\ \cos\lambda \\ 0
\end{bmatrix}
$$

- **North unit vector**
$$
\mathbf{n} =
\begin{bmatrix}
-\sin\phi\cos\lambda \\ -\sin\phi\sin\lambda \\ \cos\phi
\end{bmatrix}
$$

- **Up unit vector (ellipsoid normal)**
$$
\mathbf{u} =
\begin{bmatrix}
\cos\phi\cos\lambda \\ \cos\phi\sin\lambda \\ \sin\phi
\end{bmatrix}
$$

By stacking these row-wise, the ECEF → ENU rotation matrix is obtained:

$$
R(\phi,\lambda) =
\begin{bmatrix}
-\sin\lambda & \cos\lambda & 0 \\
-\sin\phi\cos\lambda & -\sin\phi\sin\lambda & \cos\phi \\
\cos\phi\cos\lambda & \cos\phi\sin\lambda & \sin\phi
\end{bmatrix}
$$

Therefore, for the relative vector  
$\Delta \mathbf{r} = [x - x_0,\, y - y_0,\, z - z_0]^T$ obtained in Step 1,  
the ENU coordinates are calculated as:

$$
\mathbf{enu} = R(\phi,\lambda) \, \Delta \mathbf{r}
$$

---

**From `xyz2enu.m` [lines 55–61]:**
```matlab
% Construct rotation matrix R for ECEF → ENU conversion
R = [ -sinlam          coslam         0     ; ...
      -sinphi*coslam  -sinphi*sinlam  cosphi; ...
       cosphi*coslam   cosphi*sinlam  sinphi];

% Apply the rotation matrix R to get ENU (East, North, Up) coordinates
enu = R * difxyz;
```