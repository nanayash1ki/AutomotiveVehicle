%% SIM02_SENSORS  Моделирование датчиков комплексированной системы
%
%   1) строит эталонную траекторию (модель автомобиля или тестовая);
%   2) формирует показания ИИБ, приёмника ГНСС и датчика скорости колеса;
%   3) строит графики показаний и погрешностей;
%   4) проверяет, что статистика имитаторов соответствует заданной;
%   5) сохраняет данные и рисунки в results/.
%
%   Результат sim02 — входные данные для БИНС и фильтра Калмана (sim03).

clear; close all; clc;

%% 0. Пути и параметры
proj = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(proj, 'src')), fullfile(proj, 'config'));

sp = sensor_params();
rng(sp.seed);                       % одинаковые случайные числа при каждом запуске

out_dir = fullfile(proj, 'results');
fig_dir = fullfile(out_dir, 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

set(groot, 'defaultAxesFontName', 'Times New Roman', ...
           'defaultTextFontName', 'Times New Roman', ...
           'defaultAxesFontSize', 12, 'defaultLineLineWidth', 1);

%% 1. Эталонная траектория
% 'vehicle' — модель автомобиля (vehicle_model + run_truth) по сценарию scn_drive;
% 'test'    — упрощённая кинематическая траектория (для отладки имитаторов).
TRUTH_SOURCE = 'vehicle';
T_END   = 150;          % длительность моделирования, с
DT_TRUE = 0.005;        % шаг интегрирования модели автомобиля, с

switch TRUTH_SOURCE
    case 'vehicle'
        vp = vehicle_params();
        tr = run_truth(vp, @(t, x) scn_drive(t, x, vp), T_END, DT_TRUE);
        truth = truth_from_run_truth(tr);
    case 'test'
        [t, x] = make_test_truth(DT_TRUE, T_END);
        truth = truth_from_state(t, x);
end

%% 2. Датчики
imu  = imu_model(truth, sp.imu);
gnss = gnss_model(truth, sp.gnss);
odo  = odometer_model(truth, sp.odo);

r2d = 180/pi;

%% 3. Графики
% --- Рисунок 1. Показания ИИБ ---
f1 = figure('Name', 'ИИБ: показания', 'Position', [50 50 800 700]);
tl = tiledlayout(3, 1, 'TileSpacing', 'compact');
nexttile;
plot(imu.t, imu.gyro.meas*r2d, 'Color', [0.65 0.65 0.65]); hold on;
plot(imu.t, imu.gyro.true*r2d, 'k', 'LineWidth', 1.3);
ylabel('\omega_z, °/с'); legend('измерение', 'истина', 'Location', 'best'); grid on;
title('Гироскоп (ось Z)');
nexttile;
plot(imu.t, imu.acc.meas(:,1), 'Color', [0.65 0.65 0.65]); hold on;
plot(imu.t, imu.acc.true(:,1), 'k', 'LineWidth', 1.3);
ylabel('a_x, м/с^2'); grid on; title('Акселерометр (продольная ось X)');
nexttile;
plot(imu.t, imu.acc.meas(:,2), 'Color', [0.65 0.65 0.65]); hold on;
plot(imu.t, imu.acc.true(:,2), 'k', 'LineWidth', 1.3);
ylabel('a_y, м/с^2'); grid on; title('Акселерометр (поперечная ось Y)');
xlabel(tl, 'Время, с', 'FontName', 'Times New Roman');

% --- Рисунок 2. Погрешности ИИБ ---
% систематическая часть: S*истина + b0 + дрейф (то, что будет оценивать фильтр)
g_err   = imu.gyro.meas - imu.gyro.true;
g_sys   = imu.gyro.S*imu.gyro.true + imu.gyro.b0 + imu.gyro.bdrift;
a_err   = imu.acc.meas - imu.acc.true;
a_sys   = imu.acc.S.*imu.acc.true + imu.acc.b0 + imu.acc.bdrift;

f2 = figure('Name', 'ИИБ: погрешности', 'Position', [100 50 800 700]);
tl = tiledlayout(3, 1, 'TileSpacing', 'compact');
nexttile;
plot(imu.t, g_err*r2d, 'Color', [0.7 0.7 0.7]); hold on;
plot(imu.t, g_sys*r2d, 'k', 'LineWidth', 1.3);
ylabel('\Delta\omega_z, °/с'); grid on; title('Гироскоп');
legend('полная погрешность', 'смещение и масштаб', 'Location', 'best');
for c = 1:2
    nexttile;
    plot(imu.t, a_err(:,c), 'Color', [0.7 0.7 0.7]); hold on;
    plot(imu.t, a_sys(:,c), 'k', 'LineWidth', 1.3);
    grid on;
    if c == 1, ylabel('\Deltaa_x, м/с^2'); title('Акселерометр X');
    else,      ylabel('\Deltaa_y, м/с^2'); title('Акселерометр Y'); end
end
xlabel(tl, 'Время, с', 'FontName', 'Times New Roman');

% --- Рисунок 3. Траектория и решения ГНСС ---
f3 = figure('Name', 'ГНСС: траектория', 'Position', [150 50 800 650]);
plot(truth.X, truth.Y, 'k', 'LineWidth', 1.5); hold on;
plot(gnss.pos(:,1), gnss.pos(:,2), '.', 'Color', [0.45 0.45 0.45], 'MarkerSize', 5);
in_out = ~gnss.valid;
plot(gnss.pos_true(in_out,1), gnss.pos_true(in_out,2), '--', ...
     'Color', [0.55 0.55 0.55], 'LineWidth', 3);
plot(truth.X(1), truth.Y(1), 'ko', 'MarkerFaceColor', 'k');
axis equal; grid on; xlabel('X, м'); ylabel('Y, м');
legend('истинная траектория', 'решения ГНСС', 'участок без ГНСС', 'старт', ...
       'Location', 'best');

% --- Рисунок 4. Погрешности ГНСС ---
f4 = figure('Name', 'ГНСС: погрешности', 'Position', [200 50 800 650]);
tl = tiledlayout(2, 1, 'TileSpacing', 'compact');
pe = gnss.pos_err;  pe(~gnss.valid,:) = NaN;
ve = gnss.vel_err;  ve(~gnss.valid,:) = NaN;
nexttile;
plot(gnss.t, pe(:,1), 'k'); hold on;
plot(gnss.t, pe(:,2), 'Color', [0.55 0.55 0.55]);
yline([3 -3]*sp.gnss.pos_sigma, '--k');
shade_outages(sp.gnss.outages);
ylabel('Погрешность, м'); grid on; title('Координаты');
legend('по X', 'по Y', '\pm3\sigma', 'Location', 'best');
nexttile;
plot(gnss.t, ve(:,1), 'k'); hold on;
plot(gnss.t, ve(:,2), 'Color', [0.55 0.55 0.55]);
yline([3 -3]*sp.gnss.vel_sigma, '--k');
shade_outages(sp.gnss.outages);
ylabel('Погрешность, м/с'); grid on; title('Скорость');
xlabel(tl, 'Время, с', 'FontName', 'Times New Roman');

% --- Рисунок 5. Датчик скорости колеса ---
f5 = figure('Name', 'Одометр', 'Position', [250 50 800 600]);
tl = tiledlayout(2, 1, 'TileSpacing', 'compact');
nexttile;
plot(odo.t, odo.v, 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5); hold on;
plot(odo.t, odo.v_true, 'k');
ylabel('V, м/с'); grid on; title('Скорость');
legend('измерение', 'истина', 'Location', 'best');
nexttile;
plot(odo.t, odo.v - odo.v_true, 'k');
ylabel('\DeltaV, м/с'); grid on; title('Погрешность');
xlabel(tl, 'Время, с', 'FontName', 'Times New Roman');

%% 4. Проверка имитаторов
fprintf('\n=== Проверка имитаторов датчиков ===\n');

% ИИБ: шум по участку стоянки (если он есть в начале траектории)
still = imu.t < 5 & abs(imu.gyro.true) < 1e-4 & all(abs(imu.acc.true) < 1e-3, 2);
if nnz(still) > 100
    fprintf('ИИБ (стоянка, %d отсчётов):\n', nnz(still));
    fprintf('  СКО шума гироскопа:      задано %.4f °/с,  получено %.4f °/с\n', ...
            imu.sigma_gyro*r2d, std(imu.gyro.meas(still))*r2d);
    fprintf('  СКО шума акселерометра:  задано %.4f м/с^2, получено %.4f / %.4f м/с^2\n', ...
            imu.sigma_acc, std(imu.acc.meas(still,1)), std(imu.acc.meas(still,2)));
    fprintf('  Смещение гироскопа:      %.4f °/с (среднее на стоянке %.4f °/с)\n', ...
            imu.gyro.b0*r2d, mean(imu.gyro.meas(still))*r2d);
else
    fprintf('ИИБ: в начале траектории нет стоянки — проверка шума пропущена\n');
end
fprintf('  Масштаб гироскопа S = %.2f %%, акселерометров S = %.2f / %.2f %%\n', ...
        imu.gyro.S*100, imu.acc.S*100);

% ГНСС
he = hypot(gnss.pos_err(gnss.valid,1), gnss.pos_err(gnss.valid,2));
fprintf('ГНСС (%d решений):\n', nnz(gnss.valid));
fprintf('  CEP: задано %.2f м, получено %.2f м\n', sp.gnss.CEP, median(he));
fprintf('  СКО координат: задано %.2f м, получено %.2f / %.2f м\n', sp.gnss.pos_sigma, ...
        std(gnss.pos_err(gnss.valid,1)), std(gnss.pos_err(gnss.valid,2)));
fprintf('  СКО скорости:  задано %.3f м/с, получено %.3f / %.3f м/с\n', sp.gnss.vel_sigma, ...
        std(gnss.vel_err(gnss.valid,1)), std(gnss.vel_err(gnss.valid,2)));
fprintf('  (СКО координат по одному прогону отличается от заданного:\n');
fprintf('   погрешность коррелирована с tau = %g с)\n', sp.gnss.pos_tau);

% одометр
mv = odo.v_true > 1;
rel = (odo.v(mv) - odo.v_true(mv))./odo.v_true(mv);
fprintf('Одометр:\n');
fprintf('  Погрешность радиуса k = %.2f %%: ожидается погрешность скорости -k/(1+k) = %.2f %%, получено %.2f %%\n', ...
        odo.k*100, -odo.k/(1 + odo.k)*100, mean(rel)*100);
fprintf('  СКО погрешности скорости (V > 1 м/с): %.3f м/с\n', ...
        std(odo.v(mv) - odo.v_true(mv)));

%% 5. Сохранение
save(fullfile(out_dir, 'sim02_sensors.mat'), 'truth', 'imu', 'gnss', 'odo', 'sp');
figs  = [f1 f2 f3 f4 f5];
names = {'sim02_imu', 'sim02_imu_err', 'sim02_gnss_track', 'sim02_gnss_err', 'sim02_odo'};
for i = 1:numel(figs)
    exportgraphics(figs(i), fullfile(fig_dir, [names{i} '.png']), 'Resolution', 300);
end
fprintf('\nДанные: %s\nРисунки: %s\n', fullfile(out_dir, 'sim02_sensors.mat'), fig_dir);

%% Локальные функции
function shade_outages(outages)
%SHADE_OUTAGES  Серая заливка интервалов без ГНСС на текущих осях.
yl = ylim;  xl = xlim;
outages = min(max(outages, xl(1)), xl(2));   % обрезать по длительности моделирования
for i = 1:size(outages, 1)
    if outages(i,2) <= outages(i,1), continue; end
    h = patch([outages(i,1) outages(i,2) outages(i,2) outages(i,1)], ...
              [yl(1) yl(1) yl(2) yl(2)], [0.9 0.9 0.9], ...
              'EdgeColor', 'none', 'HandleVisibility', 'off');
    uistack(h, 'bottom');
end
ylim(yl);
end
