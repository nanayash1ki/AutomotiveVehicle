function odo = odometer_model(truth, p)
%ODOMETER_MODEL  Имитатор датчика скорости колеса (зубчатое кольцо + датчик Холла).
%
%   odo = odometer_model(truth, p)
%
%   truth — эталон (truth_from_state); p — sp.odo из sensor_params.
%
%   Принцип: датчик выдаёт импульс при прохождении каждого зуба.
%   Колесо проходит за один зуб путь ds = 2*pi*r/N. Действительный радиус
%   отличается от номинального на долю k (давление, износ), поэтому
%       ds_true = 2*pi*r_nom*(1 + k)/N  — фактический путь между фронтами,
%       ds_nom  = 2*pi*r_nom/N          — значение, заложенное в вычислитель.
%   Скорость вычисляется по периоду между двумя последними фронтами
%   (а не подсчётом импульсов за 0,01 с — при 0,04 м на зуб это дало бы
%   ступеньку ~4 м/с):
%       v = ds_nom / (t_last - t_prev).
%   Время фронтов округляется до разрешения таймера. Если фронта нет дольше
%   p.timeout, скорость считается нулевой (малые скорости не измеряются).
%   Пассивный учёт направления: движение считается только вперёд.
%   Проскальзывание колеса не моделируется.
%
%   Выход:
%       odo.t      — время выдачи, с
%       odo.v      — измеренная скорость, м/с
%       odo.v_true — эталонная продольная скорость vx, м/с
%       odo.k      — фактическая погрешность радиуса (доля)
%       odo.t_edge — моменты фронтов, с

N = p.N_teeth;
odo.k = p.k_3s/3*randn;
ds_true = 2*pi*p.r_nom*(1 + odo.k)/N;
ds_nom  = 2*pi*p.r_nom/N;

% моменты фронтов: путь s(t) достигает k*ds_true
[s_u, iu] = unique(truth.s);          % s не убывает; убираем участки стоянки
n_edge = floor(s_u(end)/ds_true);
t_edge = interp1(s_u, truth.t(iu), (1:n_edge).'*ds_true, 'linear');
t_edge = round(t_edge/p.timer_res)*p.timer_res;
odo.t_edge = t_edge;

% выдача скорости с частотой fs
dt = 1/p.fs;
odo.t = (ceil(truth.t(1)*p.fs)*dt : dt : truth.t(end)).';
n = numel(odo.t);
odo.v = zeros(n, 1);
j = 0;                                % число фронтов, прошедших к моменту t(k)
for k = 1:n
    while j < n_edge && t_edge(j+1) <= odo.t(k)
        j = j + 1;
    end
    if j >= 2 && (odo.t(k) - t_edge(j)) < p.timeout
        odo.v(k) = ds_nom/(t_edge(j) - t_edge(j-1));
    end
end

tr = truth_at(truth, odo.t);
odo.v_true = tr.vx;

end
