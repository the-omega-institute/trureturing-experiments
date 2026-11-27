"""两个质量标尺 (k,d) 上的变号 Pauli-Villars 正则：三个调节质量、只消两阶发散的一参数族。

说明见“三个不同调节质量的指数和型PV正则仅在对数发散系数为负时有驻点且两种整数质量比下驻点倍率夹在一点零三九与四点四三之间而截断须把该系数调进两个窄区间才进窗.md”。

沿用 ch39-mass-scale-fixed-regulator-truncations.py 与 ch39-mass-scale-mixed-supertrace-stationarity.py 的
场表与记号：98 行有质量重场，奇荷 32 行 (g<0)、偶荷 66 行；锚 Lambda 为与 (k,d) 无关的常数，
x_r=y_s h_r^2，k=sqrt(y_e/y_o)，t=log y。每行谱函数 f'(x)=int_0^inf s^{-2} R(s) e^{-sx} ds。

一、族的结构：R(s)=1+sum_{j=1}^3 C_j e^{-a_j s}，a_j>0 两两不同，只要求 sum C_j=-1、sum C_j a_j=0。
    解集为 R_nu=F+nu p，F、p 是 S=Exp(a1)+Exp(a2)+Exp(a3) 的分布函数与密度，nu=sum C_j a_j^2/(a1 a2 a3)。
    f'_nu=A+nu B，A=int s^{-2}F e^{-sx}，B=int s^{-2}p e^{-sx}，Q=A/B 严格递减（对数凹 + 似然比序）。
二、截断：sigma(nu)=2(m1+nu m2)exp(c-1-nu m1)；进窗的 nu 恰为两个区间，端点由朗伯 W 给出。
三、精确势：驻点要求 nu 在 (-Q0,0) 内；对 a=(1,2,3) 与 (1,4,9)，在每个这样的 nu 上都有驻点，且全部驻点的
    k 在 (1.039835,4.429476) 内（逐格 Arb 认证 + 两端解析尾界），所以不进两窗必要区间 (29.9841147,35.0010644)。
"""
from pathlib import Path
from runpy import run_path
import bisect
import math
import time

from flint import arb, ctx

T0 = time.time()


def inside(x, lo, hi):
    assert arb(lo) < x < arb(hi), (x, lo, hi)
    return x


def iv(lo, hi):
    """包含 [lo,hi] 的球（lo、hi 为球）。"""
    return lo.lower().union(hi.upper())


def rmax(xs):
    m = None
    for x in xs:
        m = x if m is None else m.max(x)                     # 逐个归约，不比较可能重叠的球
    return m


def load_rows(root):
    prev = run_path(str(root/'ch39-mass-scale-spectral-stationarity.py'))
    physical = run_path(str(root/'ch39-first-three-mode-observables.py'))
    model, _ = physical['model_inputs']()
    ctx.prec, ctx.cap = 256, 3
    joint = physical['joint']
    layers = prev['layers_at'](joint['chiral'], model, prev['FLAGSHIP'])
    rows = prev['field_rows'](joint['vacuum'], layers)       # 逐行核对原 minimal['fields']
    sec = {'o': [(arb(g), h) for g, _, e, h in rows if not e], 'e': [(arb(g), h) for g, _, e, h in rows if e]}
    assert len(sec['o']) == 32 and len(sec['e']) == 66 and all(g < 0 for g, _ in sec['o'])
    return sec


# ---------------------------------------------------------------- 一、族的结构
class Shape:
    """质量平方比 a=(a1<a2<a3)。C 为最少 PV（三阶全消，R=F），v 为 sum v a^m=(0,0,1) 的解，P=a1 a2 a3。"""

    def __init__(self, al):
        a = [arb(x) for x in al]
        assert a[0] < a[1] < a[2] and a[0] > 0
        self.al, self.name = a, ','.join(str(x) for x in al)
        self.P = a[0]*a[1]*a[2]
        self.C = [-(a[(j+1) % 3]*a[(j+2) % 3])/((a[(j+1) % 3]-a[j])*(a[(j+2) % 3]-a[j])) for j in range(3)]
        self.v = [1/((a[j]-a[(j+1) % 3])*(a[j]-a[(j+2) % 3])) for j in range(3)]
        for m, want in ((0, -1), (1, 0), (2, 0)):
            assert (sum(self.C[j]*a[j]**m for j in range(3))-want).contains(0)
        for m, want in ((0, 0), (1, 0), (2, 1)):
            assert (sum(self.v[j]*a[j]**m for j in range(3))-want).contains(0)
        C, v = self.C, self.v
        self.m1 = sum(C[j]*a[j]*a[j].log() for j in range(3))          # E[1/S]=M0(F)
        self.m2 = self.P*sum(v[j]*a[j]*a[j].log() for j in range(3))   # E[1/S^2]=B(0)
        self.c = sum(C[j]*a[j].log() for j in range(3))               # E log S + gamma
        self.Q0 = self.m1/self.m2
        self.a1 = self.c-1                                              # a'(0)，a(x)=sum C(x+a)log(x+a)
        self.b1 = self.P*sum(v[j]*a[j].log() for j in range(3))        # B'(0)
        self.e1 = self.Q0*self.b1-self.a1
        assert self.m1 > 0 and self.m2 > 0

    def AB(self, x):
        a, C, v = self.al, self.C, self.v
        L = [(x+a[j])*(x+a[j]).log() for j in range(3)]
        return x*x.log()+sum(C[j]*L[j] for j in range(3)), self.P*sum(v[j]*L[j] for j in range(3))

    def ABp(self, x):
        a, C, v = self.al, self.C, self.v
        L = [(x+a[j]).log() for j in range(3)]
        return x.log()+sum(C[j]*L[j] for j in range(3)), self.P*sum(v[j]*L[j] for j in range(3))

    def dd(self, xball):
        """e''=Q0 B''-a'' 与 B'' 在 xball 上的包围。"""
        a = self.al
        app = sum(self.C[j]/(xball+a[j]) for j in range(3))
        bpp = self.P*sum(self.v[j]/(xball+a[j]) for j in range(3))
        return self.Q0*bpp-app, bpp

    def sigma(self, nu):
        return 2*(self.m1+nu*self.m2)*(self.c-1-nu*self.m1).exp()


def structure_checks(sh):
    """独立数值积分核对闭式（仅作核对）：M0、C、A、B、R=F+nu p，以及 F/p 的单调性抽样。"""
    import mpmath as mp
    mp.mp.dps = 50
    a = [mp.mpf(x.mid().str(40, radius=False)) for x in sh.al]
    C = [-(a[(j+1) % 3]*a[(j+2) % 3])/((a[(j+1) % 3]-a[j])*(a[(j+2) % 3]-a[j])) for j in range(3)]
    v = [1/((a[j]-a[(j+1) % 3])*(a[j]-a[(j+2) % 3])) for j in range(3)]
    P = a[0]*a[1]*a[2]
    def expsum(c0, cs, s):                                     # c0+sum c_j e^{-a_j s}，小 s 用泰勒级数避免相消
        if s > mp.mpf('0.05'):
            return c0+sum(cs[j]*mp.e**(-a[j]*s) for j in range(3))
        tot, n, fact = mp.mpf(0), 1, mp.mpf(1)                  # n=0,1 两阶由相消条件恰为零
        while True:
            n += 1
            fact *= n
            term = sum(cs[j]*(-a[j])**n for j in range(3))*s**n/fact
            tot += term
            if n > 6 and abs(term) < mp.mpf(10)**-45*(abs(tot)+mp.mpf(10)**-300):
                return tot
    F = lambda s: expsum(1, C, s)
    p = lambda s: expsum(0, [P*x for x in v], s)
    dens = lambda s: mp.diff(F, s)
    for s in (mp.mpf('0.3'), mp.mpf(2), mp.mpf(7)):
        assert abs(dens(s)-p(s)) < mp.mpf(10)**-20 and 0 < F(s) < 1 and p(s) > 0
    for nu in (mp.mpf('-0.3'), mp.mpf('0.7')):
        Cn = [C[j]+nu*P*v[j] for j in range(3)]
        assert abs(sum(Cn)+1) < mp.mpf(10)**-25 and abs(sum(Cn[j]*a[j] for j in range(3))) < mp.mpf(10)**-25
        assert abs(sum(Cn[j]*a[j]**2 for j in range(3))-nu*P) < mp.mpf(10)**-25
        R = lambda s: F(s)+nu*p(s)
        M0 = mp.quad(lambda s: R(s)/s**2, [0, 1, mp.inf])
        Cc = 1-mp.euler+mp.quad(lambda s: R(s)/s, [0, 1])+mp.quad(lambda s: (R(s)-1)/s, [1, mp.inf])
        sig = 2*M0*mp.e**(-Cc)
        want = sh.sigma(arb(str(nu)))
        assert abs(sig-mp.mpf(want.mid().str(30, radius=False))) < mp.mpf(10)**-15
    for x in (mp.mpf('0.05'), mp.mpf(1), mp.mpf(20)):
        A = mp.quad(lambda s: F(s)/s**2*mp.e**(-s*x), [0, 1, mp.inf])
        B = mp.quad(lambda s: p(s)/s**2*mp.e**(-s*x), [0, 1, mp.inf])
        Aa, Bb = sh.AB(arb(str(x)))
        assert abs(A-mp.mpf(Aa.mid().str(30, radius=False))) < mp.mpf(10)**-15
        assert abs(B-mp.mpf(Bb.mid().str(30, radius=False))) < mp.mpf(10)**-15
        assert P/(6*(x+a[2])**2) < A < P/(6*x**2) and P/(2*(x+a[2])) < B < P/(2*(x+a[0]))
    rat = [F(s)/p(s) for s in (mp.mpf(i)/4 for i in range(1, 80))]
    assert all(r2 > r1 for r1, r2 in zip(rat, rat[1:]))       # F/p 递增（书面证明：对数凹）


# ---------------------------------------------------------------- 三、精确势：逐格包围与尾界
class Exact:
    def __init__(self, sh, sec, Tm=-14, Tp=8, ht=1/512, hu=1/512):
        self.sh, self.sec, self.Tm, self.Tp, self.ht, self.hu = sh, sec, Tm, Tp, ht, hu
        self.Tfar = Tp+4
        self.S = {(s, p): sum(g*h**(2*p) for g, h in v) for s, v in sec.items() for p in (-1, 0, 1, 2)}
        self.L4 = {s: sum(g*h**4*(h*h).log() for g, h in v) for s, v in sec.items()}
        self.lh = {s: [((h*h).log().lower(), (h*h).log().upper()) for g, h in v] for s, v in sec.items()}
        self.hmid = {s: [arb(h.mid()) for g, h in sec[s]] for s in sec}
        self.rho = {s: [abs((h*h).log()-(hm*hm).log()).upper() for (g, h), hm in zip(sec[s], self.hmid[s])] for s in sec}
        lmin = min(float(l[0]) for s in self.lh for l in self.lh[s])-1
        lmax = max(float(l[1]) for s in self.lh for l in self.lh[s])+1
        self.u0 = math.floor((Tm+lmin)/hu)*hu
        self.nU = int(math.ceil((self.Tfar+lmax-self.u0)/hu))+2
        G = {k: [] for k in ('A', 'B', 'Ap', 'Bp', 'X')}
        for j in range(self.nU):
            x = (arb(self.u0)+arb(j)*arb(hu)).exp()
            A, B = sh.AB(x)
            Ap, Bp = sh.ABp(x)
            assert A > 0 and B > 0 and Ap < 0 and Bp < 0       # A、B 正且递减（书面证明）
            for k_, val in zip(('A', 'B', 'Ap', 'Bp', 'X'), (A, B, Ap, Bp, x)):
                G[k_].append(val)
        self.G = G

    def sector(self, s, ta, tb):
        """[ta,tb] 上：P=sum g x A、Q=sum g x B 的格包围 (Pc,Qc)、导数包围 (Pd,Qd) 与中点值 (Pm,Qm)。"""
        G, hu, u0 = self.G, self.hu, self.u0
        Em = arb((ta+tb)/2).exp()
        Pc = Qc = Pd = Qd = Pm = Qm = arb(0)
        for (g, h), (llo, lhi), hm, rr in zip(self.sec[s], self.lh[s], self.hmid[s], self.rho[s]):
            ia = int(math.floor((ta+float(llo)-u0)/hu))-1
            ib = int(math.ceil((tb+float(lhi)-u0)/hu))+1
            assert 0 <= ia < ib < self.nU
            xa, xb = G['X'][ia], G['X'][ib]
            xA = iv(xa*G['A'][ib], xb*G['A'][ia])               # A 递减、x 递增
            xB = iv(xa*G['B'][ib], xb*G['B'][ia])
            x2Ap = iv(xb*xb*G['Ap'][ia], xa*xa*G['Ap'][ib])     # A'<0 递增（A''>0）
            x2Bp = iv(xb*xb*G['Bp'][ia], xa*xa*G['Bp'][ib])
            dA, dB = xA+x2Ap, xB+x2Bp                            # d(xA)/dlog x 与 d(xB)/dlog x
            Pc += g*xA; Qc += g*xB; Pd += g*dA; Qd += g*dB
            xm = Em*hm*hm                                        # 中点 h 的精确值，h 的球半径另加 rho*|d/du|
            A, B = self.sh.AB(xm)
            Pm += g*(xm*A+rr*arb(0, abs(dA).upper()))
            Qm += g*(xm*B+rr*arb(0, abs(dB).upper()))
        return Pm, Qm, Pc, Qc, Pd, Qd

    def build(self):
        sh, Q0 = self.sh, self.sh.Q0
        half = arb(-self.ht/2).union(arb(self.ht/2))
        cells = {'o': [], 'e': []}
        nT = int(round((self.Tp-self.Tm)/self.ht))
        for s in ('o', 'e'):
            for i in range(nT):
                ta, tb = self.Tm+i*self.ht, self.Tm+(i+1)*self.ht
                Pm, Qm, Pc, Qc, Pd, Qd = self.sector(s, ta, tb)
                assert (Qc > 0) if s == 'e' else (Qc < 0)         # 偶荷 sum g x B>0：N_e 无极点
                Nd = -(Pd*Qc-Pc*Qd)/(Qc*Qc)
                N = (-Pm/Qm+half*Nd).intersection(-Pc/Qc)       # 中值定理，与单调包围取交
                cells[s].append((ta, tb, N.lower(), N.upper(), N.mid()))
        # 深尾 t<=Tm：M=N+Q0=e^t kap(-t+w+rho)/(1+delta)，|rho|<=rN，|delta|<=dD。
        Y = arb(self.Tm).exp()
        deep = {}
        for s in ('o', 'e'):
            X = arb((Y*rmax((h*h).upper() for g, h in self.sec[s])).upper())
            assert X < sh.al[0]-arb('0.5')                         # 收敛半径内
            epp, bpp = sh.dd(arb(0).union(X))
            S6abs = sum(abs(g)*h**6 for g, h in self.sec[s])
            S4, S2 = self.S[(s, 2)], self.S[(s, 1)]
            rN = (Y/2*S6abs*abs(epp).upper()/abs(S4)).upper()
            dD = ((abs(sh.b1)*Y*abs(S4)+Y*Y/2*S6abs*abs(bpp).upper())/(sh.m2*abs(S2))).upper()
            kap = S4/(sh.m2*S2)
            w = sh.e1-self.L4[s]/S4
            assert kap > 0 and dD < arb('0.1')
            assert -arb(self.Tm)+w-rN > 1                          # 两条包络在 (-inf,Tm] 上递增且为正
            deep[s] = dict(kap=kap, w=w, rN=rN, dD=dD)
        self.deep = deep
        assert -Q0+self.env('o', arb(self.Tm), 1) < 0 and -Q0+self.env('e', arb(self.Tm), 1) < 0
        # 远尾 t>=Tp：N=-e^{-t} Gamma，Gamma 为与 t 无关的包围（xA、xB 的两侧界）。
        Gam = {}
        for s in ('o', 'e'):
            numA = denB = arb(0)
            for g, h in self.sec[s]:
                X = arb(self.Tp).exp()*(h*h).lower()
                th = X/(X+sh.al[2])
                assert 0 < th < 1
                numA += g*(h**-2)*(th*th).lower().union(arb(1))
                denB += g*th.lower().union(arb(1))
            assert (denB > 0) if s == 'e' else (denB < 0)
            Gam[s] = numA/denB/3
            assert Gam[s] > 0 and (-arb(self.Tp)).exp()*Gam[s] < Q0
        self.Gam = Gam
        for s in ('o', 'e'):
            k = 0
            while self.Tp+(k+1)*self.ht <= self.Tfar+1e-12:
                a, b = self.Tp+k*self.ht, self.Tp+(k+1)*self.ht
                lo, hi = (-(-arb(a)).exp()*Gam[s]).lower(), (-(-arb(b)).exp()*Gam[s]).upper()
                cells[s].append((a, b, lo, hi, (lo+hi)/2))
                k += 1
        for s in ('o', 'e'):
            for c in cells[s]:
                assert c[2] > -Q0 and c[3] < 0                   # N_s 处处在 (-Q0,0) 内
        self.cells = cells
        return cells

    def env(self, s, t, side):
        d = self.deep[s]
        if side < 0:
            return t.exp()*d['kap']*(-t+d['w']-d['rN'])/(1+d['dD'])
        return t.exp()*d['kap']*(-t+d['w']+d['rN'])/(1-d['dD'])

    def certify(self, Dup, Dlo):
        """断言：没有驻点满足 t_e-t_o>=Dup，也没有驻点满足 t_e-t_o<=Dlo。Dup、Dlo 为二进有理数。"""
        Q0, O, E, Tm, Tp, Tfar = self.sh.Q0, self.cells['o'], self.cells['e'], self.Tm, self.Tp, self.Tfar
        Oleft, Oright = [c[0] for c in O], [c[1] for c in O]
        pref, m = [], None
        for c in O:
            m = c[3] if m is None else m.max(c[3])
            pref.append(m)
        suf, m = [None]*len(O), None
        for i in range(len(O)-1, -1, -1):
            m = O[i][2] if m is None else m.min(O[i][2])
            suf[i] = m
        Gam, do, de = self.Gam, self.deep['o'], self.deep['e']
        farO_lo = (-(-arb(Tfar)).exp()*Gam['o']).lower()         # t_o>=Tfar 时 N_o 的下界
        DU, DL, U = arb(Dup), arb(Dlo), -arb(Tm)
        gaps = []
        # 上侧：t_o<=t_e-Dup 时 N_e(t_e)>N_o(t_o)。
        for ta, tb, lo, hi, _ in E:
            tau = tb-Dup
            j = bisect.bisect_right(Oleft, tau)-1
            bound = -Q0+self.env('o', arb(min(tau, Tm)), 1)
            if j >= 0:
                bound = bound.max(pref[j])
            assert lo > bound, ('up', ta)
            gaps.append(lo-bound)
        assert DU.exp()*Gam['o'] > Gam['e']                       # 两根都在远尾
        jT = bisect.bisect_left(Oleft, Tp)-1
        assert (-(-arb(Tfar)).exp()*Gam['e']).lower() > pref[jT].max(-Q0+self.env('o', arb(Tm), 1))
        ls, rs = de['kap']*(1-do['dD'])*DU.exp(), do['kap']*(1+de['dD'])   # 两根都在深尾
        assert ls > rs and ls*(U+de['w']-de['rN']) > rs*(U+DU+do['w']+do['rN'])
        assert de['kap']*(U+de['w']-de['rN'])*(1-do['dD'])*DU.exp() > rs
        # 下侧：t_o>=t_e-Dlo 时 N_e(t_e)<N_o(t_o)。
        for ta, tb, lo, hi, _ in E:
            j = bisect.bisect_left(Oright, ta-Dlo)
            bound = suf[j].min(farO_lo) if j < len(O) else farO_lo
            if ta-Dlo < Tm:
                bound = bound.min(-Q0+self.env('o', arb(ta-Dlo), -1))
            assert hi < bound, ('lo', ta)
            gaps.append(bound-hi)
        assert -Q0+self.env('e', arb(Tm), 1) < suf[0].min(farO_lo)
        ls, rs = de['kap']*(1+do['dD'])*DL.exp(), do['kap']*(1-de['dD'])
        assert ls < rs and de['kap']*(U+de['w']+de['rN'])*(1+do['dD'])*DL.exp() < rs*(U+DL+do['w']-do['rN'])
        assert Gam['e'] > DL.exp()*Gam['o']
        return rmax(-g for g in gaps)                             # = -最小间隙

    def psi(self, s, t, nu):
        Pm, Qm, *_ = self.sector(s, t, t)
        return Pm+nu*Qm

    def witness(self, pick):
        """取格中点的近似最值位置，构造 nu 与两扇区根的严格括号，返回 (nu, k 的包围)。"""
        O, E = self.cells['o'], self.cells['e']
        om = [float(c[4]) for c in O]
        best = None
        for ta, tb, lo, hi, mid in E[:int(round((self.Tp-self.Tm)/self.ht))]:
            j = bisect.bisect_left(om, float(mid))
            if 0 < j < len(O):
                to = O[j][0]
                val = (ta+tb)/2-to
                if best is None or pick*val > pick*best[0]:
                    best = (val, (ta+tb)/2)
        te = best[1]
        Pm, Qm, *_ = self.sector('e', te, te)
        nu = arb((-Pm/Qm).mid())
        d = 1/4096
        a, b = self.psi('e', te-d, nu), self.psi('e', te+d, nu)
        assert (a < 0 < b) or (b < 0 < a)
        lo, hi = te-best[0]-0.05, te-best[0]+0.05
        fl, fh = self.psi('o', lo, nu), self.psi('o', hi, nu)
        assert (fl < 0 < fh) or (fh < 0 < fl)
        while hi-lo > d:
            md = (lo+hi)/2
            fm = self.psi('o', md, nu)
            if fm.contains(0):
                break
            if (fm < 0) == (fl < 0):
                lo, fl = md, fm
            else:
                hi = md
        assert hi-lo <= d, (lo, hi)                              # 括号长度不超过 1/4096
        assert (fl < 0 < self.psi('o', hi, nu)) or (self.psi('o', hi, nu) < 0 < fl)
        return nu, ((arb(te-d)-arb(hi))/2).exp(), ((arb(te+d)-arb(lo))/2).exp()


# ---------------------------------------------------------------- 二、截断与两窗
def window_sigma(S, L4):
    off = {s: L4[s]/S[(s, 2)] for s in ('o', 'e')}
    A = {s: S[(s, 1)]/S[(s, 2)]*off[s].exp() for s in ('o', 'e')}
    inside(A['o'], '1.45097670', '1.45097671')
    inside(A['e'], '1.1008010', '1.1008012')

    def logk(sig):                                               # (W_{-1} 奇, W_0 偶) 分支，随 sigma 严格递减
        return ((-sig*A['e']/2).lambertw()-(-sig*A['o']/2).lambertw(branch=-1)-off['e']+off['o'])/2

    def pre(K):
        a, b = arb('0.02'), arb('0.06')
        for _ in range(60):
            m = arb(((a+b)/2).mid())
            if logk(m).mid() > K.log().mid():
                a = m
            else:
                b = m
        step = arb('1e-12')
        while not (logk(a) > K.log() > logk(b)):
            a, b, step = a-step, b+step, 2*step
            assert step < arb('1e-8')
        assert logk(a) > K.log() > logk(b)                     # 严格括号（单调性见常数锚混合超迹泛函定理 2）
        return a, b
    lo = pre(arb('35.0010644'))
    hi = pre(arb('29.9841147'))
    inside(lo[0], '0.0345056616', '0.0345056620'); inside(lo[1], '0.0345056616', '0.0345056620')
    inside(hi[0], '0.0441123193', '0.0441123197'); inside(hi[1], '0.0441123193', '0.0441123197')
    return lo, hi


def truncation(sh, sig_lo, sig_hi):
    """sigma(nu)=sigma* 的两支解 nu_b=(u_b-m1)/m2，u_b=-W_b(z)/Q0，z=-Q0 sigma* e^{1-c-Q0 m1}/2。"""
    def nus(sig, br):
        z = -sh.Q0*sig*(1-sh.c-sh.Q0*sh.m1).exp()/2
        u = -z.lambertw(branch=br)/sh.Q0
        nu = (u-sh.m1)/sh.m2
        assert (sh.sigma(nu)-sig).contains(0) or abs(sh.sigma(nu)-sig) < arb('1e-40')
        return nu
    nstar = 1/sh.m1-sh.m1/sh.m2                                  # sigma 的唯一极大点
    smax = sh.sigma(nstar)
    assert smax > sh.sigma(arb(0)) > arb('0.0441123197')          # 两支都非空，nu=0 在两区间之间
    I0 = (nus(sig_lo[1], 0).lower(), nus(sig_hi[0], 0).upper())  # sigma 在 (-Q0,nstar) 上递增
    I1 = (nus(sig_hi[0], -1).lower(), nus(sig_lo[1], -1).upper())  # 在 (nstar,inf) 上递减
    # 外端：用 sigma 区间的外端点求像，得到包含进窗集合的区间；内端点同理给出被进窗集合包含的区间
    I0in = (nus(sig_lo[0], 0).upper(), nus(sig_hi[1], 0).lower())
    I1in = (nus(sig_hi[1], -1).upper(), nus(sig_lo[0], -1).lower())
    assert -sh.Q0 < I0[0] and I0[1] < 0 < I1[0]
    return dict(nstar=nstar, smax=smax, I0=I0, I1=I1, I0in=I0in, I1in=I1in)


def scan_numeric(sh, sec):
    """仅数值：点值网格上 N_o、N_e 的曲线与 k 的范围（不作证明）。"""
    ex = Exact.__new__(Exact)
    ex.sh, ex.sec = sh, sec
    rows = {s: [(g, arb(h.mid())**2) for g, h in sec[s]] for s in sec}

    def N(s, t):
        y = arb(t).exp()
        P = Q = arb(0)
        for g, h2 in rows[s]:
            A, B = sh.AB(y*h2)
            P += g*y*h2*A; Q += g*y*h2*B
        return float((-P/Q).mid())
    ts = [-12+i/16 for i in range(16*18+1)]
    No = [N('o', t) for t in ts]
    assert all(b > a for a, b in zip(No, No[1:]))
    ks = []
    for t in ts:
        ne = N('e', t)
        j = bisect.bisect_left(No, ne)
        if 0 < j < len(No):
            to = ts[j-1]+(ts[j]-ts[j-1])*(ne-No[j-1])/(No[j]-No[j-1])
            ks.append(math.exp((t-to)/2))
    return min(ks), max(ks)


def main():
    root = Path(__file__).resolve().parent
    sec = load_rows(root)
    ctx.prec = 256
    S = {(s, p): sum(g*h**(2*p) for g, h in v) for s, v in sec.items() for p in (-1, 0, 1, 2)}
    L4 = {s: sum(g*h**4*(h*h).log() for g, h in v) for s, v in sec.items()}
    assert S[('o', 1)] < 0 and S[('o', 2)] < 0 and S[('e', 1)] > 0 and S[('e', 2)] > 0
    assert S[('o', -1)] < 0 and S[('o', 0)] < 0 and S[('e', -1)] > 0 and S[('e', 0)] > 0
    assert (S[('e', 0)]-200).contains(0) and (S[('o', 0)]+448).contains(0)
    kinf = inside((S[('e', 1)]*S[('o', 2)]/(S[('e', 2)]*S[('o', 1)])).sqrt(), '2.1021534145', '2.1021534146')
    kflat = inside((S[('e', -1)]*S[('o', 0)]/(S[('e', 0)]*S[('o', -1)])).sqrt(), '4.2267177', '4.2267179')
    print('nu->-Q0 极限倍率', kinf.str(11), '  nu->0^- 极限倍率', kflat.str(8))
    sig_lo, sig_hi = window_sigma(S, L4)

    Dup, Dlo = 381/128, 10/128
    kup = inside(arb(Dup/2).exp(), '4.429475', '4.429476')
    klow = inside(arb(Dlo/2).exp(), '1.039835', '1.039836')
    expect = {
        '1,2,3': dict(m1=('0.86304621735', '0.86304621736'), m2=('1.5697444312', '1.5697444313'),
                      Q0=('0.5498004644', '0.5498004645'), s0=('1.6933171207', '1.6933171208'),
                      smax=('2.10998', '2.10999'), nstar=('0.6088860', '0.6088861'),
                      I0=('-0.5427874', '-0.5427873', '-0.5408196', '-0.5408195'),
                      I1=('7.30836', '7.30838', '7.64099', '7.64101'),
                      w0=('0.0035790', '0.0035791'), d0=('0.0127557', '0.0127558'), r1=('1.04551', '1.04552'),
                      kmax=('4.40736', '4.40887'), kmin=('1.05083', '1.05120'),
                      nuhi=('-0.0155174', '-0.0155172'), nulo=('-0.0475488', '-0.0475486')),
        '1,4,9': dict(m1=('1.34960434', '1.34960435'), m2=('4.48909320', '4.48909321'), Q0=('0.30064074', '0.30064075'),
                      s0=('1.8312756', '1.8312757'), smax=('2.49125', '2.49127'), nstar=('0.440316', '0.440318'),
                      I0=('-0.2968459', '-0.2968458', '-0.2957825', '-0.2957823'),
                      I1=('4.86851', '4.86853', '5.08025', '5.08027'),
                      w0=('0.003536', '0.003538'), d0=('0.01261', '0.01263'), r1=('1.0434', '1.0436'),
                      kmax=('4.39017', '4.39169'), kmin=('1.06426', '1.06463'),
                      nuhi=('-0.0071455', '-0.0071453'), nulo=('-0.0253611', '-0.0253609')),
    }
    for al in ((1, 2, 3), (1, 4, 9)):
        t1 = time.time()
        sh = Shape(al)
        ex = expect[sh.name]
        structure_checks(sh)
        inside(sh.m1, *ex['m1']); inside(sh.m2, *ex['m2']); inside(sh.Q0, *ex['Q0'])
        s0 = inside(sh.sigma(arb(0)), *ex['s0'])
        assert s0 > 2*(arb.const_euler()-1).exp()                 # nu=0 即最少 PV，落在已证的 sigma>=2e^{gamma-1}
        tr = truncation(sh, sig_lo, sig_hi)
        inside(tr['smax'], *ex['smax']); inside(tr['nstar'], *ex['nstar'])
        a0, b0, a1, b1 = ex['I0']
        assert arb(a0) < tr['I0'][0] < arb(b0) and arb(a1) < tr['I0'][1] < arb(b1)
        a0, b0, a1, b1 = ex['I1']
        assert arb(a0) < tr['I1'][0] < arb(b0) and arb(a1) < tr['I1'][1] < arb(b1)
        w0 = (tr['I0'][1]-tr['I0'][0])/sh.Q0
        w0in = (tr['I0in'][1]-tr['I0in'][0])/sh.Q0
        r1 = tr['I1'][1]/tr['I1'][0]
        d0 = (tr['I0'][0]+sh.Q0)/sh.Q0
        print(f'[{sh.name}] m1={sh.m1.str(11)} m2={sh.m2.str(11)} Q0={sh.Q0.str(10)} sigma(0)={s0.str(10)} '
              f'sigma_max={tr["smax"].str(6)} at nu*={tr["nstar"].str(6)}')
        print(f'  截断进窗 nu 区间 I0 包含于 ({tr["I0"][0].str(8)}, {tr["I0"][1].str(8)})，包含 ({tr["I0in"][0].str(8)}, '
              f'{tr["I0in"][1].str(8)})；相对 Q0 宽 {w0.str(4)}（内 {w0in.str(4)}），离 -Q0 的相对距离 {d0.str(4)}')
        print(f'  I1 包含于 ({tr["I1"][0].str(6)}, {tr["I1"][1].str(6)})，端点比 {r1.str(5)}')
        inside(w0, *ex['w0']); inside(w0in, *ex['w0']); inside(d0, *ex['d0']); inside(r1, *ex['r1'])
        if sh.name == '1,2,3':
            assert (sh.m1-(arb(64)/27).log()).contains(0)
        E_ = Exact(sh, sec)
        E_.build()
        worst = E_.certify(Dup, Dlo)
        assert worst < 0
        nu_hi, khlo, khhi = E_.witness(+1)
        nu_lo, kllo, klhi = E_.witness(-1)
        inside(khlo, *ex['kmax']); inside(khhi, *ex['kmax'])
        inside(kllo, *ex['kmin']); inside(klhi, *ex['kmin'])
        assert -sh.Q0 < nu_hi < 0 and -sh.Q0 < nu_lo < 0
        inside(nu_hi, *ex['nuhi']); inside(nu_lo, *ex['nulo'])
        print(f'  精确势：nu 在 (-Q0,0) 内每处都有驻点，全部驻点 {klow.str(6)} < k < {kup.str(5)}（最小间隙 {(-worst).str(3)}）')
        print(f'  见证：nu={nu_hi.str(8)} 处 k 在 ({khlo.str(7)},{khhi.str(7)})；nu={nu_lo.str(8)} 处 k 在 '
              f'({kllo.str(7)},{klhi.str(7)})；用时 {time.time()-t1:.0f} s')
    assert kup < arb('29.9841147')
    # 仅数值：其余质量比
    for al, rng in (((1, 1.01, 1.02), ('1.0457', '1.0459', '4.4149', '4.4151')),
                    ((1, 10, 100), ('1.1092', '1.1094', '4.3569', '4.3571')),
                    ((1, 100, 10000), ('1.2081', '1.2083', '4.2717', '4.2719'))):
        sh = Shape([arb(str(a)) for a in al])
        kmn, kmx = scan_numeric(sh, sec)
        assert float(rng[0]) < kmn < float(rng[1]) and float(rng[2]) < kmx < float(rng[3]), (al, kmn, kmx)
        print(f'  仅数值 a={al}: 网格上 k 的范围约 ({kmn:.4f}, {kmx:.4f})')
    print(f'完成，用时 {time.time()-T0:.0f} s')


if __name__ == '__main__':
    main()
