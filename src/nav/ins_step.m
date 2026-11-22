function [nav, F] = ins_step(nav, w_m, a_m, dt)
%INS_STEP  Один шаг БИНС на плоскости и матрица модели ошибок.
%
%   [nav, F] = ins_step(nav, w_m, a_m, dt)
%
%   nav — состояние БИНС (структура): X, Y, psi, VX, VY и текущие оценки
%         погрешностей датчиков bg, ba(2×1), Sg, Sa(2×1);
%   w_m — измеренная угловая скорость, рад/с;
%   a_m — измеренные ускорения [ax; ay] в связанной СК, м/с^2;
%   F   — матрица непрерывной модели ошибок (14×14), вычисленная
%         по данным до шага.
%
%   1) компенсация оценённых погрешностей датчиков:
%        w = w_m - bg - Sg*w_m,   a = a_m - ba - Sa.*a_m
%   2) перевод ускорения в навигационную СК:  [aN; aE] = C(psi)*a
%   3) интегрирование: курс, скорость, координаты.
%
%   Модель ошибок (ошибка = БИНС минус истина):
%     d(dpsi)/dt = -dbg - w_m*dSg
%     d(dV)/dt   = [-aE; aN]*dpsi - C(psi)*(dba + diag(a_m)*dSa)
%     d(dX)/dt   = dV
%   смещения и масштабы — постоянные (их медленное изменение — в Q).

I = struct('pos',1:2,'psi',3,'vel',4:5,'bg',6,'ba',7:8,'Sg',9,'Sa',10:11,'k',12,'eg',13:14);

w = w_m - nav.bg - nav.Sg*w_m;
a = a_m - nav.ba - nav.Sa.*a_m;

% поворот ускорения на средний за шаг курс (точнее, чем на курс в начале шага)
psi_m = nav.psi + 0.5*w*dt;
c = cos(psi_m);  s = sin(psi_m);
aN = c*a(1) - s*a(2);
aE = s*a(1) + c*a(2);

% --- матрица модели ошибок ---
F = zeros(14);
F(I.pos, I.vel) = eye(2);
F(I.psi, I.bg)  = -1;
F(I.psi, I.Sg)  = -w_m;
F(I.vel, I.psi) = [-aE; aN];
F(I.vel, I.ba)  = -[c -s; s c];
F(I.vel, I.Sa)  = -[c*a_m(1) -s*a_m(2); s*a_m(1) c*a_m(2)];

% --- интегрирование ---
nav.X   = nav.X  + nav.VX*dt + 0.5*aN*dt^2;
nav.Y   = nav.Y  + nav.VY*dt + 0.5*aE*dt^2;
nav.VX  = nav.VX + aN*dt;
nav.VY  = nav.VY + aE*dt;
nav.psi = nav.psi + w*dt;

end
