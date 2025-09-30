clearvars;
close all;

%% === ДАННЫЕ ===
% names      : {1 x m} имена каналов в порядке столбцов S
% S          : [n x m] матрица спектров LED (каждый столбец нормирован к реальному потоку)
% lambda_led : [n x 1] вектор длин волн (нм), общая сетка для LED
% lambda_t   : [nt x 1] сетка эталона
% t_raw      : [nt x 1] эталонный спектр (нормирован к реальному потоку)

[Params, T, F_Sun_AM15] = led_config();
[names, name2idx] = ra.led_order();
[S, lambda_led, lambda_t, t_raw] = spectra(names, name2idx, Params, F_Sun_AM15);
assert(exist('lambda_led','var')==1 && exist('S','var')==1, 'Нет S или lambda_led');
assert(exist('lambda_t','var')==1 && exist('t_raw','var')==1, 'Нет эталона');
[n, m] = size(S);

%% === Интерполяция эталона на сетку LED ===
t = interp1(lambda_t(:), t_raw(:), lambda_led(:), 'linear', 0);  % вне диапазона -> 0
if max(t) <= 0
    error('После интерполяции максимум эталона <= 0. Проверьте данные.');
end

%% === (Опц.) Весовая функция по длинам волн ===
u = ones(n,1);                 % сюда можно подставить V(?) или чувствительность сенсора
W = spdiags(u, 0, n, n);

%% === (Опц.) Ограничения на веса диодов ===
disabled_names = {};
active_mask = ~ismember(names, disabled_names);
% Для контроля вывод отключённых каналов:
if any(~active_mask)
    fprintf('Отключены каналы: %s\n', strjoin(names(~active_mask), ', '));
else
    fprintf('Все каналы активны.\n');
end
BIG = 1e6;                  % “большая” граница вместо Inf
w_max = BIG * ones(m,1);
w_max(~active_mask) = 0;    % отключённым — строго ноль
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

%% === ВЫВОД ВЕСОВ С ПОДПИСЯМИ И ПОМЕТКОЙ [OFF] ===
fprintf('\nОптимальные веса по каналам:\n');
for k = 1:m
    tag = '';
    if ~active_mask(k), tag = ' [OFF]'; end
    fprintf('%-14s : %.6f%s\n', names{k}, w_opt(k), tag);
end

% Таблица (удобно копировать/сохранять)
% T = table(string(names(:)), w_opt, logical(active_mask(:)), ...
%           'VariableNames', {'LED','Weight','Active'});
% disp(T);

%% === ГРАФИКИ ===
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
