%parameters
Fapp = -1; %down
phi = 0.1;

ab = 0.05;
ao = 0.0559;
bc = 0.05;
bn = 0.0559;
bp = 0.0559/2;
op = 0.0559/2;
cd = 0.05/cos(phi);
cn = 0.0559;
cm = 0.0559;
de = 0.05/cos(phi);
dm = 0.05;%approx
dl = 0.047;%approx
ef = 0.05/cos(phi);
ek = 0.043;%approx
el = 0.040;%approx
fg = 0.05;
fk = 0.038;%approx
fj = 0.0375;
gh = 0.05;
gj = 0.0375;
gi = 0.0375;
hi = 0.0375;
ij = 0.05;
jk = 0.05;
kl = 0.05;
lm = 0.05;
mn = 0.05;
no = 0.05;
np = 0.0659;%approx


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
angle_fgj = getAngle(fj, fg, gi);
angle_hgi = getAngle(hi, gh, hi);
angle_fgj = getAngle(fj, fg, fj);
angle_ghi = getAngle(gi, gh, hi);
angle_gij = getAngle(gj, ij, gi);
angle_gji = getAngle(gi, gj, ij);
angle_fjk = getAngle(fk, jk, fk);
angle_fkj = getAngle(fj, fk, jk);
angle_ekl = getAngle(el, ek, kl);
angle_elk = getAngle(ek, el, kl);
angle_dlm = getAngle(dm, dl, lm);
angle_dml = getAngle(dl, dm, lm);
angle_cmn = getAngle(cn, mn, cm);
angle_cnm = getAngle(cm, cn, mn);
angle_nop = getAngle(np, op, no);
angle_onp = getAngle(op, no, np);
angle_aob = getAngle(ab, (bp + op), ao);
angle_aop = getAngle(ab, ao, (bp+op));


clear Nax Nay Nby Tab Tao Tbc Tbn Tbp Tcd Tcm Tcn Tde Tdl Tdm Tef Tek Tel Tfg Tfk Tfj Tgj Tgi Tgh Thi Tij Tjk Tkl Tlm Tmn Tno Tnp Top
syms Nax Nay Nby Tab Tao Tbc Tbn Tbp Tcd Tcm Tcn Tde Tdl Tdm Tef Tek Tel Tfg Tfk Tfj Tgj Tgi Tgh Thi Tij Tjk Tkl Tlm Tmn Tno Tnp Top

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
0 ==  - sin(phi)*Tef + Tfk*sin(phi+angle_efk)+Tfj*sin(angle_gfj),...
0 == Tgh - Tfg + Tgi * cos(angle_hgi) - Tgj*cos(angle_fgj),...
0 == Tgi*sin(angle_hgi) + Tgj*sin(angle_fgj),...
0 == - Tgh - cos(angle_ghi)*Thi,...
0 == Thi*sin(angle_ghi) - Fapp,...
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
0 == Tno + Top*cos(angle_nop) - Tao*cos(angle_aob),...
0 == - Top*sin(angle_nop) - Tao*sin(angle_aob),...
0 == - Fapp*(bc+cd+de+ef+fg+gh) - Nay*ab,...
0 == - Fapp*(ab+bc+cd+de+ef+fg+gh) + Nby*ab

];

vars = [Nax; Nay; Nby; Tab; Tao; Tbc; Tbn; Tbp; Tcd; Tcm; Tcn; Tde; Tdl; Tdm; Tef; Tek; Tel; Tfg; Tfk; Tfj; Tgj; Tgi; Tgh; Thi; Tij; Tjk; Tkl; Tlm; Tmn; Tno; Tnp; Top];
lengths = [0; 0; 0; ab; ao; bc; bn; bp; cd; cm; cn; de; dl; dm; ef; ek; el; fg; fk; fj; gj; gi; gh; hi; ij; jk; kl; lm; mn; no; np; op];


S = solve(eqns);
solu = [double(S.Nax); double(S.Nay); double(S.Nby); double(S.Tab); double(S.Tao); double(S.Tbc); double(S.Tbn); double(S.Tbp); double(S.Tcd); double(S.Tcm); double(S.Tcn); double(S.Tde); double(S.Tdl); double(S.Tdm); double(S.Tef); double(S.Tek); double(S.Tel); double(S.Tfg); double(S.Tfk); double(S.Tfj); double(S.Tgj); double(S.Tgi); double(S.Tgh); double(S.Thi); double(S.Tij); double(S.Tjk); double(S.Tkl); double(S.Tlm); double(S.Tmn); double(S.Tno); double(S.Tnp); double(S.Top)];
