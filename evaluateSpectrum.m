function R = evaluate_spectrum(lambda_nm, spd, varargin)
% Единая “оценка спектра”: xy, uv, CCT (точный), Ra, BLH, и, опц., сравнение с target
p = inputParser;
p.addParameter('Target', [], @(v)isvector(v));
p.parse(varargin{:});
t = p.Results.Target;

% xy/uv
c = ra.xyuv(lambda_nm, spd);

% CCT (точный)
cct = ra.cctExact(lambda_nm, spd);

% CRI (если готовы TCS; иначе вернёт заглушку)
cri = ra.cri_ra(lambda_nm, spd);

% BHL индекс
blh = ra.blueHazard(lambda_nm, spd);

% Δuv / Δxy до target (если есть)
dxy = NaN; duv = NaN;
if ~isempty(t)
    ct = ra.xyuv(lambda_nm, t);
    dxy = norm(c.xy - ct.xy, 2);
    duv = norm(c.uv - ct.uv, 2);
end

% Собираем ответ
R.xy   = c.xy;
R.uv   = c.uv;
R.CCT  = cct.CCT;
R.duv  = cct.duv;
R.Ra   = cri.Ra;
R.BLH  = blh.index;
R.dxy  = dxy;
R.duv2 = duv;
end

