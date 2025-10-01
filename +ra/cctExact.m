function outCCT = cctExact(lambda_nm, spd)

% 1) u'v' для исследуемого SPD
m = ra.xyuv(lambda_nm, spd);
uv_target = m.uv;

% 2) целевая функция: Δuv(T)
fun = @(T) duv_to_blackbody(lambda_nm, uv_target, T);

% 3) минимизация по температуре
lb = 1000; ub = 25000;  % при желании расширь
opts = optimset('TolX',1e-6,'Display','off');
[Topt, fval] = fminbnd(fun, lb, ub, opts);

outCCT.CCT = Topt;
outCCT.duv = fval;
end


function duv = duv_to_blackbody(lambda_nm, uv_target, T)
% Δu'v' между SPD целевого и спектром чёрного тела при T
spd_bb = planck_spd(lambda_nm, T);
uv_bb  = ra.xyuv(lambda_nm, spd_bb).uv;
duv = hypot(uv_target(1)-uv_bb(1), uv_target(2)-uv_bb(2));
end

function spd = planck_spd(lambda_nm, T)
% Спектр чёрного тела (без конст. множителя), λ в нм
% B(λ,T) ~ (1/λ^5) * 1/(exp(c2/(λ*T))-1)
c2 = 1.438776877e7;   % нм·K (вторая константа излучения)
lam = lambda_nm(:);
spd = lam.^(-5) ./ (exp(c2./(lam*T)) - 1);
% Нормировать можно по максимуму, чтобы масштабы не мешали
spd = spd / max(spd);
end

