%variable declaration
ab = 1;
ao = 1;
bc = 1;
bn = 1;
bp = 0.5;
op = 0.5;
on = 0.5;
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
Fapp = 1; %down
phi = 0.2;
%angles:
angle_bao = getAngle((bp+op), ab, ao);
angle_cbn = getAngle(cn, bc, bn);
angle_abp = getAngle(ao, ab, (bp+op));
angle_dcm = getAngle(dm, cd, cm);
angle_bcn = getAngle(bn, bc, cn);
angle_edl = getAngle(el, de, dl);
angle_cdm = getAngle(cm, cd, dm);
angle_del = getAngle(dl, de, el);
angle_fek = getAngle(ek, ef, fk);
angle_efk = getAngle(ek, ek, fk);
angle_gfj = getAngle(gj, fg, fj);
angle_fgi = getAngle(fi, fg, gi);
angle_hgi = getAngle(hi, gh, hi);
angle_fgj = getAngle(fj, fg, fj);
angle_aop = getAngle(ab, ao, (bp+op));

clear Nax Nay Nby Tab Tao Tbc Tbo Tbp Tbn Tcd Tcm Tcn Tde Tdl Tdm Tef Tek Tel Tfg Tfk Tfj 
syms Nax Nay Nby Tab Tao Tbc Tbo Tbp Tbn Tcd Tcm Tcn Tde Tdl Tdm Tef Tek Tel Tfg Tfk Tfj 

eqns = [
0 == Nax + Tab + cos(angle_bao)*Tao,... %sum of Fx_a
0 == Nay + sin(angle_bao)*Tao,... %sum of Fy_a
0 == Tbc - Tab + Tbn*cos(angle_cbn) - cos(angle_abp)*Tbp,... %sum of Fx_b
0 == Nby + Tbn*cos(angle_cbn) + Tbc - cos(angle_abp)*Tbp,...
0 == Tcd*cos(phi) + Tcm*cos(phi+angle_dcm) - Tcn*cos(angle_bcn) - Tbc,...
0 == Tcm*sin(phi+angle_dcm) + Tcn*sin(angle_bcn) + Tcd*cos(phi),...
0 == Tde*cos(phi) - Tcd*cos(phi) + Tdl*cos(phi+angle_edl) - Tdm*cos(phi+angle_cdm),...
0 == Tde*sin(phi) - Tcd*sin(phi) + Tdl*sin(phi+angle_edl) + Tdm *sin(phi+angle_cdm),...
0 == Tef*cos(phi) - Tde*cos(phi) + Tek*cos(phi+angle_fek) - Tel*cos(phi + angle_del),...
0 == Tef*sin(phi) - Tde*sin(phi) + Tek*sin(phi+angle_fek) + Tel*sin(phi+ angle_del),...
0 == Tfg - cos(phi)*Tef - Tfk*cos(phi+angle_efk)+Tfj*cos(angle_gfj),...
%missing eq
0 == Tgh - Tfg + Tgi * cos(angle_hgi) - Tgj*cos(angle_fgi),...
0 == Tgi*sin(angle_hgi) + Tgj*sin(angle_fgj),...
0 == - Tgh - cos(angle_ghi)*Thi,...
0 == Thi*sin(angle_ghi) - F,...
0 == -Tij - Tgi * cos(angle_gij) + Thi * cos(angle_ghi),...
0 == - Tgi*sin(angle_gij) - Thi * sin(angle_ghi),...
0 == Tij - Tjk + Tgj*cos(angle_gji) - Tfj*cos(angle_fjk),...
0 == Tgj*sin(angle_gji) - Tfj*sin(angle_fjk),...

0 == Tjk - Tkl + Tfk*cos(angle_fkj) - Tek*cos(angle_ekl),...
0 == Tfk*sin(angle_fkj) - Tek*sin(angle_ekl),...

0 == Tkl - Tlm + Tel*cos(angle_elk) - Tdl*cos(angle_dlm),...
0 == Tel*sin(angle_elk) - Tdl*sin(angle_dlm),...

0 == Tlm - Tmn + Tdm*cos(angle_dml) - Tcm*cos(angle_cmn),...
0 == Tdm*sin(angle_dml) - Tcm*sin(angle_cmn),...

0 == Tmn - Tno - Tnp*cos(angle_onp) + Tcn*cos(angle_cnm) - Tbn*cos(angle_cbn),...
0 == - Tnp*sin(angle_onp) - Tcn*sin(angle_cnm) - Tbn*sin(angle_cbn),...

0 == Tno + Top*cos(angle_nop) - Tao*cos(angle_aob)
0 == - Top*sin(angle_nop) - Tao*sin(angle_aob)
];

S = solve(eqns);
S.Tab

