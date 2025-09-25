%% === ВАШИ ДАННЫЕ ===
% lambda_led : [n x 1] вектор длин волн (нм), общая сетка для LED
% S          : [n x m] матрица спектров LED (каждый столбец нормирован к 1)
% lambda_t   : [nt x 1] сетка эталона
% t_raw      : [nt x 1] эталонный спектр (НЕ нормирован)

assert(exist('lambda_led','var')==1 && exist('S','var')==1, 'Нет S или lambda_led');
assert(exist('lambda_t','var')==1 && exist('t_raw','var')==1, 'Нет эталона');
[n, m] = size(S);

%% === Интерполяция эталона на сетку LED ===
t = interp1(lambda_t(:), t_raw(:), lambda_led(:), 'linear', 0);  % вне диапазона -> 0

% (Опционально) учитывать только перекрытие:
% in_rng = (lambda_led >= min(lambda_t)) & (lambda_led <= max(lambda_t));
% t(~in_rng) = 0;

%% === Нормировка эталона к 1 ===
t_peak = max(t);
if t_peak > 0
    t = t / t_peak;
else
    error('После интерполяции максимум эталона равен нулю — проверьте данные.');
end

%% === (Опц.) Весовая функция по длинам волн ===
u = ones(n,1);                 % сюда можно подставить V(?) или чувствительность сенсора
W = spdiags(u, 0, n, n);

%% === (Опц.) Ограничения на веса диодов ===
w_max = 1e6*ones(m,1);              % при необходимости задайте пределы
disabled_idx = []; % отключённые каналы
w_max(disabled_idx) = 0;

%% === CVX: L2-фит без масштаба ===
cvx_begin quiet
    variable w(m) nonnegative
    minimize( norm( W*(S*w - t), 2 ) + 1e-3*norm(w,1) )  % лёгкая L1-регуляризация
    subject to
        w >= 0;
        w <= w_max;
        % ВКЛЮЧИТЕ, если хотите смесь долей (при пик-нормировке столбцов S):
%         sum(w) == 1
cvx_end

w_opt = w;
fit   = S*w_opt;
resid = fit - t;

%% === Метрики и графики ===
rmse = sqrt(mean(resid.^2));
mae  = mean(abs(resid));
maxe = max(abs(resid));
fprintf('Status: %s\n', cvx_status);
fprintf('RMSE: %.4g, MAE: %.4g, Max|e|: %.4g\n', rmse, mae, maxe);
disp('Оптимальные веса:')
for k = 1:m
    fprintf('%-10s : %.4f\n', names{k}, w_opt(k));
end

figure;
plot(lambda_led, t, '-', 'LineWidth', 1.5); hold on;
plot(lambda_led, fit, '--', 'LineWidth', 1.5);
xlabel('\lambda, nm'); ylabel('Relative power (norm to 1)');
legend('Target (norm)', 'Mixture S*w', 'Location', 'best'); grid on;
title('Подбор смеси LED к нормированному эталону');

figure;
plot(lambda_led, resid, 'LineWidth', 1);
xlabel('\lambda, nm'); ylabel('Error'); grid on;
title('Остаток: S*w - t');
