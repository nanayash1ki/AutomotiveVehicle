function gnss = gnss_model(truth, p)
%GNSS_MODEL  Имитатор приёмника ГНСС (координаты и скорость, 10 Гц).
%
%   gnss = gnss_model(truth, p)
%
%   truth — эталон (truth_from_state); p — sp.gnss из sensor_params.
%   Слабосвязанная схема: приёмник выдаёт готовое решение — координаты
%   X, Y и составляющие скорости VX, VY в навигационной СК.
%
%   Погрешность координат = марковская составляющая (медленная) + белый шум;
%   погрешность скорости — белый шум.
%   На интервалах p.outages решения нет: valid = false, значения NaN.
%   Задержкой решения относительно данных ИИБ пренебрегаем.
%
%   Выход:
%       gnss.t            — время, с
%       gnss.pos          — [X Y], м (NaN при отсутствии решения)
%       gnss.vel          — [VX VY], м/с
%       gnss.valid        — признак наличия решения
%       gnss.pos_true, gnss.vel_true — эталон на той же сетке
%       gnss.pos_err,  gnss.vel_err  — погрешности (есть и при пропадании)

dt = 1/p.fs;
gnss.t = (ceil(truth.t(1)*p.fs)*dt : dt : truth.t(end)).';
n = numel(gnss.t);
tr = truth_at(truth, gnss.t);

gnss.pos_true = [tr.X  tr.Y];
gnss.vel_true = [tr.VX tr.VY];

gnss.pos_err = gauss_markov(n, dt, p.pos_sigma_gm, p.pos_tau, 2) ...
               + p.pos_sigma_w*randn(n, 2);
gnss.vel_err = p.vel_sigma*randn(n, 2);

gnss.valid = true(n, 1);
for i = 1:size(p.outages, 1)
    gnss.valid(gnss.t >= p.outages(i,1) & gnss.t < p.outages(i,2)) = false;
end

gnss.pos = gnss.pos_true + gnss.pos_err;
gnss.vel = gnss.vel_true + gnss.vel_err;
gnss.pos(~gnss.valid, :) = NaN;
gnss.vel(~gnss.valid, :) = NaN;

end
