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

%angles:
angle_bao = getAngle(bo, ab, ao);

clear Tab Tao
syms Tab Tao

eqns = [
Tab == 1, 
0 == Tab + cos(angle_bao)*Tao
];

S = solve(eqns);
S.Tab

