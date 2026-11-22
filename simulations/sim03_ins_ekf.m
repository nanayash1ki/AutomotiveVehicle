%% SIM03_INS_EKF  Комплексированная навигационная система: БИНС + ГНСС + одометр
%
%   Берёт показания датчиков из sim02 (results/sim02_sensors.mat) и
%   сравнивает три варианта обработки:
%     1) автономная БИНС (без коррекции);
%     2) БИНС + ГНСС (фильтр Калмана, без одометра);
%     3) БИНС + ГНСС + одометр (полная система).
%   Строит графики ошибок, оценок погрешностей датчиков, печатает таблицу
%   точности и сохраняет результаты в results/.

clear; close all; clc;

%% 0. Пути и данные
proj = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(proj, 'src')), fullfile(proj, 'config'));
f_in = fullfile(proj, 'results', 'sim02_sensors.mat');
if ~exist(f_in, 'file'), error('Нет %s — сначала запустите sim02_sensors', f_in); end
load(f_in, 'truth', 'imu', 'gnss', 'odo', 'sp');

fig_dir = fullfile(proj, 'results', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
set(groot, 'defaultAxesFontName', 'Times New Roman', ...
           'defaultTextFontName', 'Times New Roman', ...
           'defaultAxesFontSize', 12, 'defaultLineLineWidth', 1);
r2d = 180/pi;

%% 1. Три варианта обработки
ep = ekf_params(sp);
psi0 = truth.psi(1) + ep.psi0_err;           % начальный курс с ошибкой

ep_free = ep;  ep_free.use_gnss = false;  ep_free.use_odo = false;
ep_noodo = ep; ep_noodo.use_odo = false;

tic;
est_free = ekf_ins(imu, gnss, odo, ep_free,  psi0);
est_gnss = ekf_ins(imu, gnss, odo, ep_noodo, psi0);
est      = ekf_ins(imu, gnss, odo, ep,       psi0);
fprintf('Расчёт: %.1f с\n', toc);

tr = truth_at(truth, est.t);                 % истина на сетке БИНС
t  = est.t;
herr = @(e) hypot(e.X - tr.X, e.Y - tr.Y);   % горизонтальная ошибка, м
e_free = herr(est_free);  e_gnss = herr(est_gnss);  e_full = herr(est);
g_err  = hypot(gnss.pos_err(:,1), gnss.pos_err(:,2));
g_err(~gnss.valid) = NaN;

out = sp.gnss.outages(1,:);                  % первый интервал без ГНСС

%% 2. Таблица точности
in_out  = t >= out(1) & t < out(2);
with_g  = t > 10 & ~in_out & ~(t >= out(2) & t < out(2) + 10);  % ГНСС есть, после сходимости
k_end   = find(t < out(2), 1, 'last');
gv      = gnss.valid & gnss.t > 10;
fprintf('\n=== Точность определения координат (горизонтальная ошибка) ===\n');
fprintf('%-30s %14s %16s %18s\n', 'Вариант', 'СКО с ГНСС, м', 'макс. без ГНСС, м', 'в конце пропадания, м');
fprintf('%-30s %14.2f %16s %18s\n', 'Только ГНСС', sqrt(mean(g_err(gv).^2)), '—', '—');
rows = {'Автономная БИНС', e_free; 'БИНС + ГНСС', e_gnss; 'БИНС + ГНСС + одометр', e_full};
for i = 1:size(rows, 1)
    e = rows{i,2};
    fprintf('%-30s %14.2f %16.2f %18.2f\n', rows{i,1}, sqrt(mean(e(with_g).^2)), ...
            max(e(in_out)), e(k_end));
end
fprintf('Автономная БИНС: ошибка через 30 / 60 / 100 с = %.0f / %.0f / %.0f м\n', ...
        e_free(find(t >= 30, 1)), e_free(find(t >= 60, 1)), e_free(find(t >= 100, 1)));
fprintf('Длительность пропадания ГНСС %.0f с, путь за это время %.0f м\n', ...
        diff(out), tr.s(k_end) - tr.s(find(t >= out(1), 1)));

% согласованность: доля времени, когда ошибка вне ±3σ
ex = est.X - tr.X;  sx = 3*sqrt(est.Pd(:, ep.I.pos(1)));
epsi = est.psi - tr.psi;  spsi = 3*sqrt(est.Pd(:, ep.I.psi));
m = t > 10;
fprintf('Вне границ ±3σ: координата X — %.1f %% времени, курс — %.1f %%\n', ...
        100*mean(abs(ex(m)) > sx(m)), 100*mean(abs(epsi(m)) > spsi(m)));

%% 3. Графики
gray = [0.6 0.6 0.6];

% --- Рисунок 1. Траектории ---
f1 = figure('Name', 'Траектории', 'Position', [50 50 800 650]);
plot(tr.X, tr.Y, 'k', 'LineWidth', 2); hold on;
plot(gnss.pos(:,1), gnss.pos(:,2), '.', 'Color', gray, 'MarkerSize', 5);
plot(est.X, est.Y, '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 1.3);
plot(est_free.X, est_free.Y, ':k', 'LineWidth', 1.3);
plot(tr.X(in_out), tr.Y(in_out), 'k', 'LineWidth', 4);
axis equal; grid on;
pad = 0.1*max(max(tr.X) - min(tr.X), max(tr.Y) - min(tr.Y));
xlim([min(tr.X) - pad, max(tr.X) + pad]);  ylim([min(tr.Y) - pad, max(tr.Y) + pad]);
xlabel('X, м'); ylabel('Y, м');
legend('истинная', 'ГНСС', 'БИНС + ГНСС + одометр', 'автономная БИНС', ...
       'участок без ГНСС', 'Location', 'best');

% --- Рисунок 2. Ошибки координат: сравнение вариантов ---
f2 = figure('Name', 'Ошибки координат', 'Position', [100 50 800 700]);
tl = tiledlayout(2, 1, 'TileSpacing', 'compact');
nexttile;
plot(t, e_free, ':k', 'LineWidth', 1.5); hold on;
plot(t, e_full, 'k', 'LineWidth', 1.3);
shade(out); grid on; ylabel('Ошибка, м');
title('Автономная БИНС и комплексированная система');
legend('автономная БИНС', 'БИНС + ГНСС + одометр', 'Location', 'northwest');
nexttile;
plot(gnss.t, g_err, '.', 'Color', gray, 'MarkerSize', 4); hold on;
plot(t, e_gnss, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.3);
plot(t, e_full, 'k', 'LineWidth', 1.3);
shade(out); grid on; ylabel('Ошибка, м');
ylim([0, max([e_gnss(t > 10); 8])*1.1]);
title('Комплексированная система (увеличенный масштаб)');
legend('только ГНСС', 'БИНС + ГНСС', 'БИНС + ГНСС + одометр', 'Location', 'northwest');
xlabel(tl, 'Время, с', 'FontName', 'Times New Roman');

% --- Рисунок 3. Ошибки по осям и границы ±3σ ---
f3 = figure('Name', 'Ошибки и 3 сигма', 'Position', [150 50 800 700]);
tl = tiledlayout(3, 1, 'TileSpacing', 'compact');
ey = est.Y - tr.Y;  sy = 3*sqrt(est.Pd(:, ep.I.pos(2)));
nexttile; plot_3s(t, ex, sx, out); ylabel('\DeltaX, м'); title('Координата X');
nexttile; plot_3s(t, ey, sy, out); ylabel('\DeltaY, м'); title('Координата Y');
nexttile; plot_3s(t, epsi*r2d, spsi*r2d, out); ylabel('\Delta\psi, °'); title('Курс');
xlabel(tl, 'Время, с', 'FontName', 'Times New Roman');

% --- Рисунок 4. Ошибки скорости ---
f4 = figure('Name', 'Ошибки скорости', 'Position', [200 50 800 500]);
evx = est.VX - tr.VX;  evy = est.VY - tr.VY;
plot(t, evx, 'k'); hold on; plot(t, evy, 'Color', gray);
shade(out); grid on; xlabel('Время, с'); ylabel('Ошибка, м/с');
ylim([-1 1]*max(0.2, 1.2*max(abs([evx(t > 10); evy(t > 10)]))));
legend('V_X', 'V_Y', 'Location', 'best');

% --- Рисунок 5. Оценки погрешностей датчиков ---
f5 = figure('Name', 'Оценки погрешностей датчиков', 'Position', [250 50 900 750]);
tl = tiledlayout(3, 2, 'TileSpacing', 'compact');
b_true = interp1(imu.t, imu.gyro.b0 + imu.gyro.bdrift, t);
nexttile; plot_est(t, est.bg*r2d, b_true*r2d, out); ylabel('°/с'); title('Смещение гироскопа');
nexttile; plot_est(t, est.Sg*100, imu.gyro.S*100*ones(size(t)), out); ylabel('%'); title('Масштаб гироскопа');
ba_true = interp1(imu.t, imu.acc.b0 + imu.acc.bdrift, t);
nexttile; plot_est(t, est.ba(:,1), ba_true(:,1), out); ylabel('м/с^2'); title('Смещение акселерометра X');
nexttile; plot_est(t, est.ba(:,2), ba_true(:,2), out); ylabel('м/с^2'); title('Смещение акселерометра Y');
nexttile; plot_est(t, est.Sa*100, imu.acc.S*100.*ones(size(t)), out); ylabel('%');
title('Масштабы акселерометров X, Y');
nexttile; plot_est(t, est.k*100, odo.k*100*ones(size(t)), out); ylabel('%'); title('Масштаб одометра');
lg = legend('оценка', 'истинное значение'); lg.Layout.Tile = 'south'; lg.Orientation = 'horizontal';
xlabel(tl, 'Время, с', 'FontName', 'Times New Roman');

%% 4. Сохранение
save(fullfile(proj, 'results', 'sim03_ins_ekf.mat'), 'est', 'est_gnss', 'est_free', 'ep', 'tr');
figs  = [f1 f2 f3 f4 f5];
names = {'sim03_track', 'sim03_pos_err', 'sim03_err_3sigma', 'sim03_vel_err', 'sim03_sensor_est'};
for i = 1:numel(figs)
    exportgraphics(figs(i), fullfile(fig_dir, [names{i} '.png']), 'Resolution', 300);
end
fprintf('\nРисунки: %s\n', fig_dir);

%% Локальные функции
function shade(out)
%SHADE  Серая полоса на интервале без ГНСС (на текущих осях, под графиками).
yl = ylim;
h = patch([out(1) out(2) out(2) out(1)], [yl(1) yl(1) yl(2) yl(2)], [0.9 0.9 0.9], ...
          'EdgeColor', 'none', 'HandleVisibility', 'off');
uistack(h, 'bottom');  ylim(yl);
end

function plot_3s(t, e, s3, out)
%PLOT_3S  Ошибка оценки и границы ±3σ из ковариационной матрицы фильтра.
plot(t, e, 'k'); hold on;
plot(t, s3, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2);
plot(t, -s3, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, 'HandleVisibility', 'off');
m = t > 10;  ylim([-1.3 1.3]*max(s3(m)));
shade(out); grid on;
legend('ошибка', '\pm3\sigma', 'Location', 'northwest');
end

function plot_est(t, x_est, x_true, out)
%PLOT_EST  Оценка параметра и его истинное значение.
plot(t, x_est, 'k', 'LineWidth', 1.2); hold on;
plot(t, x_true, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2);
v = [reshape(x_est(t > 5, :), [], 1); x_true(:)];
yl = [min(v), max(v)];
ylim(yl + [-0.2 0.2]*max(diff(yl), eps));
shade(out); grid on;
end
