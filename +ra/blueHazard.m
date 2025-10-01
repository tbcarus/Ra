
function out = blueHazard(lambda_nm, spd, varargin)
% Индекс "blue-light hazard" как ∫ SPD(λ)*B(λ) dλ (нормированный).
% Name-Value:
%   'B' : [n x 1] весовая функция B(λ) на той же сетке (если нет — используем встроенную)

p = inputParser;
p.addParameter('B', [], @(v)isvector(v));
p.parse(varargin{:});
B = p.Results.B;

if isempty(B)
    B = default_B(lambda_nm); % вставь свою стандартную функцию
end
B = B(:);

num = trapz(lambda_nm, spd(:).*B);
den = trapz(lambda_nm, spd(:));         % нормируем на общую мощность
index = num / max(den, eps);

out.index = index;
out.weighted = num;      % невнормированный интеграл
out.B = B;
end

function B = default_B(lambda_nm)
% ЗАГЛУШКА: вставь стандартную кривую blue-light hazard B(λ) (≈ 300–700 нм)
B = zeros(numel(lambda_nm),1);
% Пример: можно задать плечо 400–500нм "колокольчиком"/таблицей и интерполировать
end


