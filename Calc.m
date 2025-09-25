%% === ���� ������ ===
% lambda_led : [n x 1] ������ ���� ���� (��), ����� ����� ��� LED
% S          : [n x m] ������� �������� LED (������ ������� ���������� � 1)
% lambda_t   : [nt x 1] ����� �������
% t_raw      : [nt x 1] ��������� ������ (�� ����������)

assert(exist('lambda_led','var')==1 && exist('S','var')==1, '��� S ��� lambda_led');
assert(exist('lambda_t','var')==1 && exist('t_raw','var')==1, '��� �������');
[n, m] = size(S);

%% === ������������ ������� �� ����� LED ===
t = interp1(lambda_t(:), t_raw(:), lambda_led(:), 'linear', 0);  % ��� ��������� -> 0

% (�����������) ��������� ������ ����������:
% in_rng = (lambda_led >= min(lambda_t)) & (lambda_led <= max(lambda_t));
% t(~in_rng) = 0;

%% === ���������� ������� � 1 ===
t_peak = max(t);
if t_peak > 0
    t = t / t_peak;
else
    error('����� ������������ �������� ������� ����� ���� � ��������� ������.');
end

%% === (���.) ������� ������� �� ������ ���� ===
u = ones(n,1);                 % ���� ����� ���������� V(?) ��� ���������������� �������
W = spdiags(u, 0, n, n);

%% === (���.) ����������� �� ���� ������ ===
w_max = 1e6*ones(m,1);              % ��� ������������� ������� �������
disabled_idx = []; % ����������� ������
w_max(disabled_idx) = 0;

%% === CVX: L2-��� ��� �������� ===
cvx_begin quiet
    variable w(m) nonnegative
    minimize( norm( W*(S*w - t), 2 ) + 1e-3*norm(w,1) )  % ����� L1-�������������
    subject to
        w >= 0;
        w <= w_max;
        % ��������, ���� ������ ����� ����� (��� ���-���������� �������� S):
%         sum(w) == 1
cvx_end

w_opt = w;
fit   = S*w_opt;
resid = fit - t;

%% === ������� � ������� ===
rmse = sqrt(mean(resid.^2));
mae  = mean(abs(resid));
maxe = max(abs(resid));
fprintf('Status: %s\n', cvx_status);
fprintf('RMSE: %.4g, MAE: %.4g, Max|e|: %.4g\n', rmse, mae, maxe);
disp('����������� ����:')
for k = 1:m
    fprintf('%-10s : %.4f\n', names{k}, w_opt(k));
end

figure;
plot(lambda_led, t, '-', 'LineWidth', 1.5); hold on;
plot(lambda_led, fit, '--', 'LineWidth', 1.5);
xlabel('\lambda, nm'); ylabel('Relative power (norm to 1)');
legend('Target (norm)', 'Mixture S*w', 'Location', 'best'); grid on;
title('������ ����� LED � �������������� �������');

figure;
plot(lambda_led, resid, 'LineWidth', 1);
xlabel('\lambda, nm'); ylabel('Error'); grid on;
title('�������: S*w - t');
