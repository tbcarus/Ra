function [E, eItmo] = evaluateSpectrum(L_interp, spd, target)
%   L_interp : [n x 1] сетка длин волн (нм)
%   spd      : [n x 1] спектр смеси (например, fit)
%   target : [n x 1] эталонный спектр на той же сетке (опц.)

% E.color.XYZ — тристимульные значения (отн.).
% E.color.xy — координаты CIE 1931 xy
% E.color.uv — координаты CIE 1976 u'v'
% E.CCT.CCT — коррелированная цветовая температура (K).
% E.CCT.duv — модуль минимальной дистанции до планковского локуса в u′v′.
% E.CCT.uv_src — координаты u′v′ твоего источника.
% E.CCT.uv_bb — координаты u′v′ чёрного тела при E.CCT.CCT.
% E.CCT.T_bounds — интервал поиска температуры, [K].
% E.CCT.method — описание метода расчёта.
% E.compare.dxy — расстояние между xy fit и target (если target задан).
% E.compare.duv — расстояние между u′v′ fit и target.


% базовые метрики по текущему спектру xy/uv
col = ra.xyuv(L_interp, spd);
E.color = col;

% базовые метрики по текущему спектру xy/uv - альтернативный
colItmo = itmo.xyuvItmo(L_interp, spd);
eItmo.color = colItmo;

% Точная CCT и duv до локуса
cct = ra.cctExact(L_interp, spd);
E.CCT = cct;        % Кельвины

% Точная CCT - альтернативный
cctItmo = itmo.cctItmo(L_interp, spd);
eItmo.CCT = cctItmo;

% Ra
cri = ra.cri_ra(L_interp, spd);
E.CRI = cri;

% Сравнение с target (эталонный спектр)
E.compare.dxy = NaN;
E.compare.duv = NaN;
if ~isempty(target)
    col_t = ra.xyuv(L_interp, target);
    cct_t = ra.cctExact(L_interp, target);
    E.target.color = col_t;     % цветиметрия эталона
    E.target.CCT   = cct_t;     % CCT эталона (аналогичная структура)
    E.compare.dxy = norm(col.xy - col_t.xy, 2);
    E.compare.duv = norm(col.uv - col_t.uv, 2);
else
    % если таргета нет — возвращаем пустые структуры для единообразия
    E.target = struct('color',[],'CCT',[]);
end

% Ra альтернативный
eItmoCri = itmo.criRaItmo(L_interp, spd);
eItmo.CRI = eItmoCri;

end

