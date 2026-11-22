function imu = imu_model(truth, p)
%IMU_MODEL  Имитатор ИИБ (гироскоп Z, акселерометры X и Y).
%
%   imu = imu_model(truth, p)
%
%   truth — эталон (truth_from_state); p — sp.imu из sensor_params.
%   Движение плоское, поэтому моделируются три канала, нужные БИНС 2D:
%   угловая скорость рыскания и продольное/поперечное ускорения.
%
%   Модель каждого канала:
%       meas = (1 + S)*true + b0 + bdrift(t) + noise
%   S, b0      — случайные постоянные (свои в каждом запуске);
%   bdrift(t)  — марковский дрейф нуля (gauss_markov);
%   noise      — белый шум, СКО = N*sqrt(fs/2).
%
%   Выход (все в СИ):
%       imu.t                     — время, с
%       imu.gyro.meas, .true      — угловая скорость, рад/с (n×1)
%       imu.acc.meas,  .true      — ускорения [ax ay], м/с^2 (n×2)
%       imu.gyro.S, .b0, .bdrift, .noise   — составляющие погрешности
%       imu.acc.S,  .b0, .bdrift, .noise   (по столбцам X, Y)
%       imu.sigma_gyro, imu.sigma_acc      — СКО шума

dt = 1/p.fs;
imu.t = (ceil(truth.t(1)*p.fs)*dt : dt : truth.t(end)).';
n = numel(imu.t);
tr = truth_at(truth, imu.t);

imu.sigma_gyro = p.gyro.N*sqrt(p.fs/2);
imu.sigma_acc  = p.acc.N *sqrt(p.fs/2);

%% гироскоп
g = p.gyro;
if p.calibrated, b0_3s = g.b0cal_3s; else, b0_3s = g.b0_3s; end
imu.gyro.true   = tr.r;
imu.gyro.S      = g.S_3s/3*randn;
imu.gyro.b0     = b0_3s/3*randn;
imu.gyro.bdrift = gauss_markov(n, dt, g.bi_sigma, g.bi_tau);
imu.gyro.noise  = imu.sigma_gyro*randn(n, 1);
imu.gyro.meas   = (1 + imu.gyro.S)*imu.gyro.true + imu.gyro.b0 ...
                  + imu.gyro.bdrift + imu.gyro.noise;

%% акселерометры X, Y
a = p.acc;
if p.calibrated, b0_3s = a.b0cal_3s; else, b0_3s = a.b0_3s; end
imu.acc.true   = [tr.ax tr.ay];
imu.acc.S      = a.S_3s/3*randn(1, 2);
imu.acc.b0     = b0_3s/3*randn(1, 2);
imu.acc.bdrift = gauss_markov(n, dt, a.bi_sigma, a.bi_tau, 2);
imu.acc.noise  = imu.sigma_acc*randn(n, 2);
imu.acc.meas   = (1 + imu.acc.S).*imu.acc.true + imu.acc.b0 ...
                 + imu.acc.bdrift + imu.acc.noise;

end
