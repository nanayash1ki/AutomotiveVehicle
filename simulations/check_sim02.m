%% CHECK_SIM02  Количественная проверка эталона и имитаторов датчиков
%
%   Запускать ПОСЛЕ sim02_sensors (берёт results/sim02_sensors.mat).
%   Каждая проверка печатает «OK» или «!! ПРОВЕРИТЬ».
%
%   1. Эталон: величины согласованы между собой (курс — интеграл угловой
%      скорости, ускорения — производные скоростей). Ловит ошибки систем
%      координат (мировая/связанная).
%   2. Счисление по ИДЕАЛЬНЫМ датчикам повторяет траекторию — значит,
%      ИИБ получает правильные сигналы. Для сравнения — счисление по
%      реальным (зашумлённым) показаниям ИИБ: так ведёт себя БИНС без коррекции.
%   3. Датчики: параметры, оценённые по данным, совпадают с заданными.

clear; clc;
proj = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(proj, 'src')), fullfile(proj, 'config'));
load(fullfile(proj, 'results', 'sim02_sensors.mat'), 'truth', 'imu', 'gnss', 'odo', 'sp');
r2d = 180/pi;
rel = @(err, sig) sqrt(mean(err.^2))/sqrt(mean(sig.^2));   % относительная СКО

%% 1. Эталон
fprintf('\n=== 1. Эталон: согласованность величин ===\n');
t = truth.t;
fprintf('  Курс в конце: %.1f° (по сценарию scn_drive: 90°)\n', truth.psi(end)*r2d);
fprintf('  Пройденный путь: %.0f м, макс. скорость %.1f м/с\n', truth.s(end), max(truth.V));

e_psi = truth.psi(1) + cumtrapz(t, truth.r) - truth.psi;
report('курс = интеграл угловой скорости, макс. ошибка', max(abs(e_psi))*r2d, 0.1, '°');
report('ax = dvx/dt - r*vy, относит. СКО', ...
       rel(gradient(truth.vx, t) - truth.r.*truth.vy - truth.ax, truth.ax), 0.03, '');
report('ay = dvy/dt + r*vx, относит. СКО', ...
       rel(gradient(truth.vy, t) + truth.r.*truth.vx - truth.ay, truth.ay), 0.03, '');
report('VX = dX/dt, относит. СКО', rel(gradient(truth.X, t) - truth.VX, truth.VX), 0.01, '');

%% 2. Счисление по идеальным и по реальным показаниям ИИБ
fprintf('\n=== 2. Счисление пути по ускорениям и угловой скорости ===\n');
[Xi, Yi] = dead_reckon(t, truth.r, truth.ax, truth.ay, truth);
err_ideal = hypot(Xi - truth.X, Yi - truth.Y);
report('по идеальным датчикам: макс. ошибка координат', max(err_ideal), 3, 'м');
fprintf('    (остаток — погрешность численного интегрирования, не датчиков)\n');

[Xm, Ym] = dead_reckon(imu.t, imu.gyro.meas, imu.acc.meas(:,1), imu.acc.meas(:,2), truth);
tr_imu = truth_at(truth, imu.t);
err_meas = hypot(Xm - tr_imu.X, Ym - tr_imu.Y);
fprintf('  по реальным показаниям ИИБ (БИНС без коррекции):\n');
for tt = [30 60 100 imu.t(end)]
    [~, k] = min(abs(imu.t - tt));
    fprintf('    ошибка через %5.0f с: %8.1f м\n', imu.t(k), err_meas(k));
end
fprintf('    (именно этот уход должен устранить фильтр Калмана)\n');

%% 3. Датчики
fprintf('\n=== 3. Датчики: оценка по данным против заданного ===\n');

% --- ИИБ: шум по первым разностям (смещение и медленный дрейф в разностях сокращаются)
eg = imu.gyro.meas - imu.gyro.true;
ea = imu.acc.meas  - imu.acc.true;
report_ratio('гироскоп: СКО шума, °/с', std(diff(eg))/sqrt(2)*r2d, imu.sigma_gyro*r2d, 0.05);
report_ratio('акселерометр X: СКО шума, м/с^2', std(diff(ea(:,1)))/sqrt(2), imu.sigma_acc, 0.05);
report_ratio('акселерометр Y: СКО шума, м/с^2', std(diff(ea(:,2)))/sqrt(2), imu.sigma_acc, 0.05);

% --- ИИБ: смещение и масштаб по регрессии «погрешность = b0 + S*истина»
n = numel(imu.t);
c = [ones(n,1) imu.gyro.true] \ eg;
report('гироскоп: смещение, оценка минус заданное, °/с', ...
       (c(1) - imu.gyro.b0 - mean(imu.gyro.bdrift))*r2d, 0.005, '°/с');
report('гироскоп: масштаб, оценка минус заданное, %', (c(2) - imu.gyro.S)*100, 0.2, '%');
lbl = 'XY';
for ax = 1:2
    c = [ones(n,1) imu.acc.true(:,ax)] \ ea(:,ax);
    report(sprintf('акселерометр %s: смещение, оценка минус заданное', lbl(ax)), ...
           c(1) - imu.acc.b0(ax) - mean(imu.acc.bdrift(:,ax)), 1e-3, 'м/с^2');
end

% --- ГНСС: 100 прогонов (погрешность коррелирована, один прогон мало что говорит)
M = 100;
rng(100);
he_all = [];  ex_all = [];  ev_all = [];
for i = 1:M
    g = gnss_model(truth, sp.gnss);
    he_all = [he_all; hypot(g.pos_err(:,1), g.pos_err(:,2))]; %#ok<AGROW>
    ex_all = [ex_all; g.pos_err(:)];                          %#ok<AGROW>
    ev_all = [ev_all; g.vel_err(:)];                          %#ok<AGROW>
end
report_ratio(sprintf('ГНСС: CEP по %d прогонам, м', M), median(he_all), sp.gnss.CEP, 0.05);
report_ratio('ГНСС: СКО координаты по оси, м', std(ex_all), sp.gnss.pos_sigma, 0.05);
report_ratio('ГНСС: СКО скорости по оси, м/с', std(ev_all), sp.gnss.vel_sigma, 0.05);
report('ГНСС: доля решений вне ±3σ, %', 100*mean(abs(ex_all) > 3*sp.gnss.pos_sigma), 0.6, '%');
dur = sum(diff(sp.gnss.outages, 1, 2));
fprintf('  ГНСС: решений нет %.1f с (задано %.1f с)\n', sum(~gnss.valid)/sp.gnss.fs, dur);

% --- Одометр: v_изм = ds_nom/T = v/(1+k), поэтому средняя относительная
%     погрешность должна быть около -k/(1+k) (знак минус!)
mv = odo.v_true > 1;
rel_odo = mean((odo.v(mv) - odo.v_true(mv))./odo.v_true(mv));
report('одометр: средняя относит. погрешность минус ожидаемая, %', ...
       (rel_odo + odo.k/(1 + odo.k))*100, 0.2, '%');
fprintf('    (k = %.2f %%, ожидаемая погрешность %.2f %%, получено %.2f %%)\n', ...
        odo.k*100, -odo.k/(1 + odo.k)*100, rel_odo*100);
stop = odo.v_true < 1e-3;
if any(stop)
    report('одометр: скорость на стоянке, м/с', max(abs(odo.v(stop))), 1e-9, 'м/с');
end

fprintf('\nГотово. Все строки «OK» — эталон и имитаторы работают как задумано.\n');

%% Локальные функции
function [X, Y] = dead_reckon(t, w, ax, ay, truth)
%DEAD_RECKON  Счисление на плоскости: курс, скорость, координаты (метод трапеций).
s0  = truth_at(truth, t(1));
psi = s0.psi + cumtrapz(t, w);
aN  = ax.*cos(psi) - ay.*sin(psi);        % ускорение в навигационной СК
aE  = ax.*sin(psi) + ay.*cos(psi);
VX  = s0.VX + cumtrapz(t, aN);
VY  = s0.VY + cumtrapz(t, aE);
X   = s0.X  + cumtrapz(t, VX);
Y   = s0.Y  + cumtrapz(t, VY);
end

function report(name, val, tol, unit)
fprintf('  %-58s %10.4f %-5s (допуск ±%g)  %s\n', name, val, unit, tol, passfail(abs(val) <= tol));
end

function report_ratio(name, got, want, reltol)
fprintf('  %-58s %10.4f (задано %.4f, ±%g%%)  %s\n', name, got, want, reltol*100, ...
        passfail(abs(got/want - 1) <= reltol));
end

function s = passfail(ok)
if ok, s = 'OK'; else, s = '!! ПРОВЕРИТЬ'; end
end
