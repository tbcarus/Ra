function xyuvOut = xyuv(lambda_nm, spd)
[x,y,z] = ra.cmf1931_2deg(lambda_nm);

X = trapz(lambda_nm, spd.*x);
Y = trapz(lambda_nm, spd*y);
Z = trapz(lambda_nm, spd*z);

sumXYZ = max(X+Y+Z, eps);
xy = [X;Y]/sumXYZ;

den = max(X + 15*Y + 3*Z, eps);
up = 4*X/den;
vp = 9*Y/den;

xyuvOut.XYZ = [X;Y;Z];
xyuvOut.xy  = xy;
xyuvOut.uv  = [up; vp];

end

