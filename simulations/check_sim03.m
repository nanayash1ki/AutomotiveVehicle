%% CHECK_SIM03  Проверка БИНС и фильтра Калмана
%
%   Запускать после sim02_sensors (нужен results/sim02_sensors.mat).
%   Каждая проверка печатает «OK» или «!! ПРОВЕРИТЬ».
%
%   1. БИНС по идеальным датчикам повторяет траекторию
%      -> правильны уравнения счисления (ins_step).
%   2. Основной прогон: ошибки лежат в границах ±3σ, оценки погрешностей
%      датчиков сходятся к истинным -> фильтр согласован.
%   3. Большая ошибка начального курса устраняется
%      -> правильны знаки невязок и обратной связи (при ошибке знака
%         фильтр расходится).
%   4. Монте-Карло: 20 прогонов с новыми шумами датчиков. Средняя
%      нормированная квадратичная ошибка (e^2/σ^2) должна быть около 1:
%      много больше 1 — фильтр самоуверен, много меньше 1 — слишком осторожен.

clear; clc;
proj = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(proj, 'src')), fullfile(proj, 'config'));
load(fullfile(proj, 'results', 'sim02_sensors.mat'), 'truth', 'imu', 'gnss', 'odo', 'sp');
ep = ekf_params(sp);
I  = ep.I;
r2d = 180/pi;  d2r = pi/180;

%% 1. БИНС по идеальным датчикам
fprintf('\n=== 1. БИНС по идеальным датчикам (без коррекции) ===\n');
imu_i = imu;   imu_i.gyro.meas = imu.gyro.true;  imu_i.acc.meas = imu.acc.true;
gnss_i = gnss; gnss_i.pos = gnss.pos_true;       gnss_i.vel = gnss.vel_true;  % только для начальных условий
ep_free = ep;  ep_free.use_gnss = false;  ep_free.use_odo = false;
est = ekf_ins(imu_i, gnss_i, odo, ep_free, truth.psi(1));
tr  = truth_at(truth, est.t);
report('макс. ошибка координат за весь прогон, м', max(hypot(est.X - tr.X, est.Y - tr.Y)), 3);
report('макс. ошибка курса, °', max(abs(est.psi - tr.psi))*r2d, 0.1);
fprintf('    (остаток — погрешность численного интегрирования на пути %.0f м)\n', tr.s(end));

%% 2. Основной прогон
fprintf('\n=== 2. Основной прогон: согласованность и оценки ===\n');
est = ekf_ins(imu, gnss, odo, ep, truth.psi(1) + ep.psi0_err);
t = est.t;  m = t > 10;                   % после начальной сходимости
err = [est.X - tr.X, est.Y - tr.Y, est.psi - tr.psi, est.VX - tr.VX, est.VY - tr.VY];
idx = [I.pos, I.psi, I.vel];
names = {'X', 'Y', 'курс', 'VX', 'VY'};
for i = 1:5
    s3 = 3*sqrt(est.Pd(:, idx(i)));
    report(sprintf('%s: доля времени вне ±3σ, %%', names{i}), 100*mean(abs(err(m,i)) > s3(m)), 5);
end

% оценки погрешностей датчиков в конце прогона: |оценка - истина| <= 3σ
p_est  = [est.bg(end), est.ba(end,:), est.Sg(end), est.Sa(end,:), est.k(end)];
p_true = [imu.gyro.b0 + imu.gyro.bdrift(end), imu.acc.b0 + imu.acc.bdrift(end,:), ...
          imu.gyro.S, imu.acc.S, odo.k];
p_idx  = [I.bg, I.ba, I.Sg, I.Sa, I.k];
p_name = {'смещение гироскопа', 'смещение акселерометра X', 'смещение акселерометра Y', ...
          'масштаб гироскопа', 'масштаб акселерометра X', 'масштаб акселерометра Y', ...
          'масштаб одометра'};
fprintf('  Оценки погрешностей датчиков в конце прогона (|оценка - истина| / 3σ, должно быть <= 1):\n');
for i = 1:numel(p_est)
    ratio = abs(p_est(i) - p_true(i))/(3*sqrt(est.Pd(end, p_idx(i))));
    report(['  ' p_name{i}], ratio, 1);
end

%% 3. Сходимость при большой ошибке начального курса
fprintf('\n=== 3. Устранение ошибки начального курса ===\n');
for dpsi = [5 10 20]
    e3 = ekf_ins(imu, gnss, odo, ep, truth.psi(1) + dpsi*d2r);
    k20 = find(t >= 20, 1);
    report(sprintf('начальная ошибка %2d° -> ошибка курса через 20 с, °', dpsi), ...
           abs(e3.psi(k20) - tr.psi(k20))*r2d, 1);
end

%% 4. Монте-Карло
M = 20;
fprintf('\n=== 4. Монте-Карло, %d прогонов (около минуты) ===\n', M);
nees = zeros(M, 3);  out_max = zeros(M, 1);  frac = zeros(M, 2);
outage = sp.gnss.outages(1,:);
in_out = t >= outage(1) & t < outage(2);
for i = 1:M
    rng(1000 + i);
    imu_r  = imu_model(truth, sp.imu);
    gnss_r = gnss_model(truth, sp.gnss);
    odo_r  = odometer_model(truth, sp.odo);
    e4 = ekf_ins(imu_r, gnss_r, odo_r, ep, truth.psi(1) + ep.psi0_err);
    ex = e4.X - tr.X;  ey = e4.Y - tr.Y;  epsi = e4.psi - tr.psi;
    nees(i,:) = [mean(ex(m).^2./e4.Pd(m, I.pos(1))), mean(ey(m).^2./e4.Pd(m, I.pos(2))), ...
                 mean(epsi(m).^2./e4.Pd(m, I.psi))];
    frac(i,:) = [mean(abs(ex(m)) > 3*sqrt(e4.Pd(m, I.pos(1)))), ...
                 mean(abs(epsi(m)) > 3*sqrt(e4.Pd(m, I.psi)))]*100;
    out_max(i) = max(hypot(ex(in_out), ey(in_out)));
    fprintf('.');
end
fprintf('\n');
report_range('средняя e^2/σ^2 по X (норма около 1)',    mean(nees(:,1)), 0.5, 2);
report_range('средняя e^2/σ^2 по Y (норма около 1)',    mean(nees(:,2)), 0.5, 2);
report_range('средняя e^2/σ^2 по курсу (норма около 1)', mean(nees(:,3)), 0.5, 2);
report('X: доля времени вне ±3σ, среднее по прогонам, %', mean(frac(:,1)), 3);
report('курс: доля времени вне ±3σ, среднее по прогонам, %', mean(frac(:,2)), 3);
fprintf('  Макс. ошибка за время без ГНСС: среднее %.2f м, худший прогон %.2f м\n', ...
        mean(out_max), max(out_max));

fprintf('\nГотово.\n');

%% Локальные функции
function report(name, val, tol)
fprintf('  %-62s %9.3f (допуск %g)  %s\n', name, val, tol, passfail(abs(val) <= tol));
end

function report_range(name, val, lo, hi)
fprintf('  %-62s %9.3f (норма %g…%g)  %s\n', name, val, lo, hi, passfail(val >= lo && val <= hi));
end

function s = passfail(ok)
if ok, s = 'OK'; else, s = '!! ПРОВЕРИТЬ'; end
end
