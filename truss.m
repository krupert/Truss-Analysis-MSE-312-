%variable declaration
ab = 1;
ao = 1;
bc = 1;
bo = 1;
bn = 1;
cd = 1;
cn = 1;
cm = 1;
de = 1;
dm = 1;
dl = 1;
ef = 1;
ek = 1;
el = 1;
fg = 1;
fk = 1;
fj = 1;
gh = 1;
gj = 1;
gi = 1;
hi = 1;
ij = 1;
jk = 1;
kl = 1;
lm = 1;
mn = 1;
no = 1;
F = 1; %down
phi = 0.2;
%angles:
angle_bao = getAngle(bo, ab, ao);
angle_cbn = getAngle(cn, bc, bn);
angle_dcm = getAngle(dm, cd, cm);

clear Nax Nay Nby Tab Tao Tbc Tbo Tbn Tcn
syms Nax Nay Nby Tab Tao Tbc Tbo Tbn Tcn

eqns = [
0 == Nax + Tab + cos(angle_bao)*Tao, %sum of Fx_a
0 == Nay + sin(angle_bao)*Tao, %sum of Fy_a
0 == Tbn*cos(angle_cbn) + Tbc - cos(angle_bao)*Tbo - Tab, %sum of Fx_b
0 == Nby + Tbn*cos(angle_cbn) + Tbc - cos(angle_bao),

];

S = solve(eqns);
S.Tab

