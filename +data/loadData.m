function D = loadData(L_interp)
% Функция для централизованной загрузки данных (пока не используется)
D.lambda = L_interp(:);
[D.x, D.y, D.z] = data.cieXyz(L_interp);
D.R     = data.rObjects(L_interp);
[D.S0, D.S1, D.S2] = data.s0s1s2(L_interp);
D.B     = data.blh(L_interp);
end

