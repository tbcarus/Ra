function outCRI = criRaItmo(Lambda, I)

c = 3e8;
h = 6.6262e-34;
k = 1.38067e-23;
LL = Lambda*1e-9;
[xCurve,yCurve,zCurve] = data.cie_xyz(Lambda);
kc = 100/trapz(Lambda, I.*yCurve);
xyuv = itmo.xyuvItmo(Lambda, I);
X = xyuv.X;
Y = xyuv.Y;
Z = xyuv.Z;
u = xyuv.uv(1);
v = xyuv.uv(2);
TxcK = itmo.cctItmo(Lambda, I).CCT;

%% Ra около кривой планка с цветовым различием не более 0,01
plank = 2*h*c^2./(LL.^5.*(exp(h*c./(LL*k*TxcK))-1)); % эталонный источник
kc_ref = 100/trapz(LL, plank.*yCurve);
Xref = kc_ref * trapz(LL, plank.*xCurve);
Yref = kc_ref * trapz(LL, plank.*yCurve);
Zref = kc_ref * trapz(LL, plank.*zCurve);

Luv_ref = Yref;
u_ref = 4*Xref/(Xref+15*Yref+3*Zref);
v_ref = 9*Yref/(Xref+15*Yref+3*Zref);

R = data.rObjects(Lambda);
% R = R';

XRtest = zeros(1,14); YRtest = zeros(1,14); ZRtest = zeros(1,14);
LRtest = zeros(1,14); uRtest = zeros(1,14); vRtest = zeros(1,14);
LxRtest = zeros(1,14); uxRtest = zeros(1,14); vxRtest = zeros(1,14);

XRref = zeros(1,14); YRref = zeros(1,14); ZRref = zeros(1,14);
LRref = zeros(1,14); uRref = zeros(1,14); vRref = zeros(1,14);
LxRref = zeros(1,14); uxRref = zeros(1,14); vxRref = zeros(1,14);
dE = zeros(1, 14); CRI = zeros(1, 14);
for i=1:14
    % координаты цвета эталонных поврехностей, освещаемых исследуемым
    % (тестируемым) источником
    XRtest(1,i) = kc*trapz(Lambda, I.*R(:,i).*xCurve);
    YRtest(1,i) = kc*trapz(Lambda, I.*R(:,i).*yCurve);
    ZRtest(1,i) = kc*trapz(Lambda, I.*R(:,i).*zCurve);
    LRtest(1, i) = YRtest(1,i);
    uRtest(1, i) = 4*XRtest(1,i)/(XRtest(1,i)+15*YRtest(1,i)+3*ZRtest(1,i));
    vRtest(1, i) = 9*YRtest(1,i)/(XRtest(1,i)+15*YRtest(1,i)+3*ZRtest(1,i));
    LxRtest(1, i) = 116*(YRtest(1,i)/Y).^(1/3) - 16;
    uxRtest(1, i) = 13*LxRtest(1, i)*(uRtest(1, i)-u);
    vxRtest(1, i) = 13*LxRtest(1, i)*(vRtest(1, i)-v);

    % координаты цвета эталонных поврехностей, освещаемых эталонным
    % (референсным) источником
    XRref(1,i) = kc_ref*trapz(LL, plank.*R(:,i).*xCurve);
    YRref(1,i) = kc_ref*trapz(LL, plank.*R(:,i).*yCurve);
    ZRref(1,i) = kc_ref*trapz(LL, plank.*R(:,i).*zCurve);
    LRref(1, i) = YRref(1,i);
    uRref(1, i) = 4*XRref(1,i)/(XRref(1,i)+15*YRref(1,i)+3*ZRref(1,i));
    vRref(1, i) = 9*YRref(1,i)/(XRref(1,i)+15*YRref(1,i)+3*ZRref(1,i));
    LxRref(1, i) = 116*(YRref(1,i)/Y).^(1/3) - 16;
    uxRref(1, i) = 13*LxRref(1, i)*(uRref(1, i)-u);
    vxRref(1, i) = 13*LxRref(1, i)*(vRref(1, i)-v);
    
    % Цветовая разность
    dE(1, i) = sqrt(...
        (LxRtest(1, i)-LxRref(1, i))^2 + ...
        (uxRtest(1, i)-uxRref(1, i))^2 + ...
        (vxRtest(1, i)-vxRref(1, i))^2 ...
        );
    
    % Частные индексы цветопередачи
    CRI(1, i) = 100 - 4.6*dE(1, i);
end
    CRIo = mean(CRI);

%% Ra около кривой планка с цветовым различием более 0,01
plank2 = 2*h*c^2./(LL.^5.*(exp(h*c./(LL*k*TxcK))-1)); % эталонный источник
kc_ref2 = 100/trapz(LL, plank2.*yCurve);
Xref2 = kc_ref2 * trapz(LL, plank2.*xCurve);
Yref2 = kc_ref2 * trapz(LL, plank2.*yCurve);
Zref2 = kc_ref2 * trapz(LL, plank2.*zCurve);

Luv_ref2 = Yref2;
u_ref2 = 4*Xref2/(Xref2+15*Yref2+3*Zref2);
v_ref2 = 9*Yref2/(Xref2+15*Yref2+3*Zref2);

XRref2 = zeros(1,14); YRref2 = zeros(1,14); ZRref2 = zeros(1,14);
LRref2 = zeros(1,14); uRref2 = zeros(1,14); vRref2 = zeros(1,14);
LXRref2 = zeros(1,14); uXRref2 = zeros(1,14); vXRref2 = zeros(1,14);

rref2 = (4 - u_ref2 - 10 * v_ref2)/v_ref2;
fref2 = (1.708*v_ref2 + 0.404 - 1.481*u_ref2)/v_ref2;
rtest2 = (4 - u - 10 * v)/v;
ftest2 = (1.708*v + 0.404 - 1.481*u)/v;
rR2 = zeros(1,14); fR2 = zeros(1,14);
uxxRtest = zeros(1,14); vxxRtest = zeros(1,14);
LxxxRref = zeros(1,14); uxxxRref = zeros(1,14); vxxxRref = zeros(1,14);
LxxxRtest = zeros(1,14); uxxxRtest = zeros(1,14); vxxxRtest = zeros(1,14);
dE2 = zeros(1, 14); CRI2 = zeros(1, 14);
for i=1:14
    % координаты цвета эталонных поврехностей, освещаемых эталонным
    % (референсным) источником
    XRref2(1,i) = kc_ref2*trapz(LL, plank2.*R(:,i).*xCurve);
    YRref2(1,i) = kc_ref2*trapz(LL, plank2.*R(:,i).*yCurve);
    ZRref2(1,i) = kc_ref2*trapz(LL, plank2.*R(:,i).*zCurve);
    LRref2(1, i) = YRref2(1,i);
    uRref2(1, i) = 4*XRref2(1,i)/(XRref2(1,i)+15*YRref2(1,i)+3*ZRref2(1,i));
    vRref2(1, i) = 9*YRref2(1,i)/(XRref2(1,i)+15*YRref2(1,i)+3*ZRref2(1,i));
    LXRref2(1, i) = 116*(YRref2(1,i)/Y).^(1/3) - 16;
    uXRref2(1, i) = 13*LXRref2(1, i)*(uRref2(1, i)-u);
    vXRref2(1, i) = 13*LXRref2(1, i)*(vRref2(1, i)-v);

    rR2(1, i) = (4 - uRtest(1, i) - 10 * vRtest(1, i))/vRtest(1, i);
    fR2(1, i) = (1.708*vRtest(1, i) + 0.404 - 1.481*uRtest(1, i))/vRtest(1, i);

    uxxRtest(1, i) = (10.872+0.404*(rref2/rtest2)*rR2(1, i) - 4*(fref2/ftest2)*fR2(1, i))/(16.518+1.481*(rref2/rtest2)*rR2(1, i)-(fref2/ftest2)*fR2(1, i));
    vxxRtest(1, i) = 5.520/(16.518+1.481*(rref2/rtest2)*rR2(1, i)-(fref2/ftest2)*fR2(1, i));

    
    uxxtest = (10.872+0.404*rref2 - 4*fref2)/(16.518+1.481*rref2-fref2);
    vxxtest = 5.520/(16.518+1.481*rref2-fref2);

    LxxxRref(1, i) = (116*(YRref2(1,i)/Yref2)^(1/3)-16);
    uxxxRref(1, i) = 13*LxxxRref(1, i)*(uRref2(1, i)-u_ref2);
    vxxxRref(1, i) = 13*LxxxRref(1, i)*(vRref2(1, i)-v_ref2);

    LxxxRtest(1, i) = (116*(YRtest(1,i)/Y)^(1/3)-16);
    uxxxRtest(1, i) = 13*LxxxRtest(1, i)*(uxxRtest(1, i)-uxxtest);
    vxxxRtest(1, i) = 13*LxxxRtest(1, i)*(vxxRtest(1, i)-vxxtest);

% Цветовая разность
    dE2(1, i) = sqrt(...
        (LxxxRref(1, i)-LxxxRtest(1, i))^2 + ...
        (uxxxRref(1, i)-uxxxRtest(1, i))^2 + ...
        (vxxxRref(1, i)-vxxxRtest(1, i))^2 ...
        );
    
    % Частные индексы цветопередачи
    CRI2(1, i) = 100 - 4.6*dE2(1, i);
end
CRIo2 = mean(CRI2);

outCRI.nearPlank = CRIo;
outCRI.notNearPlnk = CRIo2;

end

