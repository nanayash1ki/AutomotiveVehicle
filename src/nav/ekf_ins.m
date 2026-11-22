function est = ekf_ins(imu, gnss, odo, ep, psi0)
%EKF_INS  БИНС 2D + расширенный фильтр Калмана по ошибкам с обратной связью.
%
%   est = ekf_ins(imu, gnss, odo, ep, psi0)
%
%   imu, gnss, odo — показания датчиков (sim02); ep — ekf_params;
%   psi0 — начальный курс, рад.
%   Режимы задаются флагами ep.use_gnss, ep.use_gnss_vel, ep.use_odo.
%   При всех флагах false получается автономная БИНС (без коррекции).
%
%   Цикл на каждом такте ИИБ (200 Гц):
%     1) шаг БИНС и матрица модели ошибок F            (ins_step)
%     2) прогноз ковариации P = Phi*P*Phi' + Q,  Phi = I + F*dt
%     3) если пришло измерение ГНСС / одометра — невязка z, коррекция:
%        K = P*H'/(H*P*H' + R),  dx = K*z,  P — в форме Джозефа
%     4) обратная связь: оценки ошибок вычитаются из БИНС, dx обнуляется.
%
%   Выход est: t, X, Y, psi, VX, VY, оценки bg, ba, Sg, Sa, k, eG
%   и Pd — диагональ ковариационной матрицы на каждом такте (n×14).

I  = ep.I;
dt = imu.t(2) - imu.t(1);              % шаг БИНС = период ИИБ
n  = numel(imu.t);

%% начальное состояние: координаты и скорость — по первому решению ГНСС
j0 = find(gnss.valid, 1);
nav.X  = gnss.pos(j0,1);   nav.Y  = gnss.pos(j0,2);
nav.VX = gnss.vel(j0,1);   nav.VY = gnss.vel(j0,2);
nav.psi = psi0;
nav.bg = 0;  nav.ba = [0; 0];  nav.Sg = 0;  nav.Sa = [0; 0];
nav.k  = 0;  nav.eG = [0; 0];
P = ep.P0;

%% память под результаты
est.t = imu.t;
[est.X, est.Y, est.psi, est.VX, est.VY, est.bg, est.Sg, est.k] = deal(zeros(n,1));
[est.ba, est.Sa, est.eG] = deal(zeros(n,2));
est.Pd = zeros(n, ep.n);

jg = 1;  jo = 1;                       % указатели на следующие измерения
ng = numel(gnss.t);  no = numel(odo.t);
tol = 1e-9;

for k = 1:n
    %% 1–2. БИНС и прогноз
    if k > 1
        [nav, F] = ins_step(nav, imu.gyro.meas(k), imu.acc.meas(k,:).', dt);
        Phi = eye(ep.n) + F*dt;
        Phi(I.eg, I.eg) = ep.phi_e*eye(2);
        nav.eG = ep.phi_e*nav.eG;
        P = Phi*P*Phi.' + ep.Q;
    end

    %% 3. коррекция по ГНСС
    while jg <= ng && gnss.t(jg) <= imu.t(k) + tol
        if ep.use_gnss && gnss.valid(jg)
            % координаты: ГНСС = истина + медленная погрешность + шум
            H = zeros(2, ep.n);  H(:, I.pos) = eye(2);  H(:, I.eg) = eye(2);
            z = [nav.X; nav.Y] + nav.eG - gnss.pos(jg,:).';
            [dx, P] = kf_update(P, z, H, ep.R_pos);
            nav = feedback(nav, dx, I);
            if ep.use_gnss_vel
                H = zeros(2, ep.n);  H(:, I.vel) = eye(2);
                z = [nav.VX; nav.VY] - gnss.vel(jg,:).';
                [dx, P] = kf_update(P, z, H, ep.R_vel);
                nav = feedback(nav, dx, I);
            end
        end
        jg = jg + 1;
    end

    %% 3. коррекция по одометру (продольная скорость в связанной СК)
    while jo <= no && odo.t(jo) <= imu.t(k) + tol
        if ep.use_odo
            v_odo = odo.v(jo);
            c = cos(nav.psi);  s = sin(nav.psi);
            vxb =  c*nav.VX + s*nav.VY;          % скорость БИНС вдоль оси автомобиля
            vyb = -s*nav.VX + c*nav.VY;
            H = zeros(1, ep.n);
            H(I.psi) = vyb;  H(I.vel) = [c s];  H(I.k) = -v_odo;
            z = vxb - v_odo*(1 + nav.k);         % одометр с учётом оценки масштаба
            [dx, P] = kf_update(P, z, H, ep.R_odo);
            nav = feedback(nav, dx, I);
        end
        jo = jo + 1;
    end

    %% запись
    est.X(k) = nav.X;   est.Y(k) = nav.Y;   est.psi(k) = nav.psi;
    est.VX(k) = nav.VX; est.VY(k) = nav.VY;
    est.bg(k) = nav.bg; est.ba(k,:) = nav.ba.';
    est.Sg(k) = nav.Sg; est.Sa(k,:) = nav.Sa.';
    est.k(k)  = nav.k;  est.eG(k,:) = nav.eG.';
    est.Pd(k,:) = diag(P).';
end

end

%% ---------------------------------------------------------------
function [dx, P] = kf_update(P, z, H, R)
%KF_UPDATE  Коррекция фильтра Калмана (ковариация — в форме Джозефа).
S  = H*P*H.' + R;
K  = P*H.'/S;
dx = K*z;
A  = eye(size(P)) - K*H;
P  = A*P*A.' + K*R*K.';
end

function nav = feedback(nav, dx, I)
%FEEDBACK  Обратная связь: вычесть оценки ошибок из БИНС и оценок датчиков.
%   Ошибка = оценка минус истина, поэтому  истина = оценка - ошибка.
%   После вызова оценка ошибки считается нулевой (dx в фильтре не хранится).
nav.X   = nav.X   - dx(I.pos(1));
nav.Y   = nav.Y   - dx(I.pos(2));
nav.psi = nav.psi - dx(I.psi);
nav.VX  = nav.VX  - dx(I.vel(1));
nav.VY  = nav.VY  - dx(I.vel(2));
nav.bg  = nav.bg  - dx(I.bg);
nav.ba  = nav.ba  - dx(I.ba);
nav.Sg  = nav.Sg  - dx(I.Sg);
nav.Sa  = nav.Sa  - dx(I.Sa);
nav.k   = nav.k   - dx(I.k);
nav.eG  = nav.eG  - dx(I.eg);
end
