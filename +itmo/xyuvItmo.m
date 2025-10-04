function [xyuvOut] = xyuvItmo(Lambda, I)

% XYZ
[xCurve,yCurve,zCurve] = ra.cie_xyz(Lambda);
kc = 100/trapz(Lambda, I.*yCurve);
X = kc * trapz(Lambda, I.*xCurve);
Y = kc * trapz(Lambda, I.*yCurve);
Z = kc * trapz(Lambda, I.*zCurve);
XYZ = [X; Y; Z];

x = X/(X+Y+Z);
y = Y/(X+Y+Z);

% Lu'v'
Luv = Y;
u = 4*X/(X+15*Y+3*Z);
v = 9*Y/(X+15*Y+3*Z);

xyuvOut.X = X;
xyuvOut.Y = Y;
xyuvOut.Z = Z;
xyuvOut.XYZ = XYZ;
xyuvOut.xy  = [x; y];
xyuvOut.uv  = [u; v];

end

