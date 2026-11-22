function u = scn_drive(t, x, vp)
%SCN_DRIVE  Сценарий движения для run_truth: управление u = [F_drive; delta].
%
%   u = scn_drive(t, x, vp)
%   Использование:  truth = run_truth(vp, @(t,x) scn_drive(t,x,vp), 150, 0.005);
%
%   Водитель задаётся двумя регуляторами:
%   - тяга F_drive отслеживает заданную скорость v_ref(t)
%     (упреждение по ускорению + П-регулятор + компенсация сопротивлений);
%   - угол руля delta отслеживает заданную угловую скорость рыскания r_ref(t)
%     (кинематическое упреждение L*r/V + П-регулятор).
%
%   Профиль: старт со скоростью vp.v0 (при v0 = 0 — стоянка 5 с, по ней
%   проверяется шум ИИБ), разгон до 15 м/с, повороты влево и вправо на 90°,
%   змейка, торможение до 8 м/с, поворот влево на 90°. Часть манёвров
%   приходится на интервал пропадания ГНСС 45–85 с.
%
%   x = [X Y phi vx vy dphi], vx, vy — в мировой СК (как в vehicle_model).

d2r = pi/180;
g = 9.81;
prm = vp.prm;
f_roll = getf(prm, 'f_roll', 0.015);
rho = getf(prm, 'rho', 1.225);  Cd = getf(prm, 'Cd', 0.32);  A = getf(prm, 'A', 2.2);
L_wb = vp.h_a + vp.h_b;

% --- заданная скорость, м/с ---
tv = [0 5 15 105 115 1e4];
vv = [vp.v0 vp.v0 15 15 8 8];
v_ref = interp1(tv, vv, t);
a_ref = (interp1(tv, vv, t + 0.05) - interp1(tv, vv, max(t - 0.05, 0)))/0.1;

% --- заданная угловая скорость рыскания: [начало конец °/с] ---
turns = [ 30  40   9.0;
          55  65  -9.0;
          75  85   4.5;
          85  95  -4.5;
         120 135   6.0];
ramp = 1.0;                                 % нарастание, с
r_ref = 0;
for i = 1:size(turns, 1)
    up   = (t - turns(i,1) + ramp/2)/ramp;
    down = (turns(i,2) + ramp/2 - t)/ramp;
    r_ref = r_ref + turns(i,3)*d2r*min(max(min(up, down), 0), 1);
end

% --- продольная скорость в связанной СК ---
phi = x(3);
Vxb = x(4)*cos(phi) + x(5)*sin(phi);

% --- тяга ---
K_v = 1.0;                                  % 1/с
F = vp.m*(a_ref + K_v*(v_ref - Vxb)) + 0.5*rho*Cd*A*Vxb*abs(Vxb) + f_roll*vp.m*g;
F = min(max(F, -0.8*vp.m*g), 0.5*vp.m*g);   % ограничения тяги и торможения
if v_ref < 0.1 && Vxb < 0.2
    F = 0;                                  % стоянка
end

% --- руль ---
K_r = 0.1;                                  % с
delta = L_wb*r_ref/max(Vxb, 1.0) + K_r*(r_ref - x(6));
delta = min(max(delta, -0.5), 0.5);         % ±28,6°

u = [F; delta];
end

function v = getf(s, name, def)
if isfield(s, name), v = s.(name); else, v = def; end
end
