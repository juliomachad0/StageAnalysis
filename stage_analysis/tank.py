from numpy import pi
"""     For a spherical tank     """


def ts(P_int, r, Sw):
    """Spherical wall thickness"""
    return (P_int * r) / (2 * Sw)


def Vst(r):
    """For a spherical tank"""
    return 4 * pi * (r ^ 3) / 3


def calcular_raio_esfera(volume):
    return ((3 * volume) / (4 * pi)) ** (1 / 3)


def Ws(r, ts, rho):
    return 4 * pi * r ^ 2 * ts * rho


"""     For a tank with ellipsoidal and spherical ends     """


def b(V, r):
    """ V : Volume de propelente,
        r : raio interno do estágio"""
    return 3 * V / 2 * pi * r ^ 2


def K(r, b):
    return r / b


def Wse(r, ts, rho):
    """Spherical tank end weight"""
    return 2 * pi * (r ^ 2) * ts * rho


def V_ee(r, b):
    """For ellipsoidal end"""
    return 2 * (r ^ 2) * b / 3


def V_se(r):
    """For spherical end"""
    return 2 * pi * (r ^ 3) / 3


def tc(P_int, r, Sw):
    return P_int * r / Sw


def lc(Vp, r, V_ee, V_se, tc):
    a = r - tc
    Vc = Vp - V_ee - V_se
    return Vc / a


def Wc(r, lc, tc, rho):
    return 2 * pi * r * lc * tc * rho


def tse(P_int, r, K, Sw):
    return (P_int * r * K + 0.5) / 2 * Sw


def tk(P_int, K, a, Sw):
    return K * P_int * a / Sw


def tcr(P_int, R, Sw):
    return P_int * R / 2 * Sw


def tee(tk, tcr):
    return (tk + tcr) / 2


def Wee(a, te, E, rho, K):
    return pi * a ^ 2 * te * E * rho / 2 * K


def Propellant_tank(r, Vp, P_int, rho, Sw, Fs, E):
    """ r           :   raio externo
        Vp          :   Volume de propelente
        P_int       :   pressão de trabalho interna
        rho         :   densidade do material utilizado
        Sw          :   limite de tensão do material
        Fs          :   Fator de segurança adotado
        E           :   

        A função retorna uma lista 9 itens, onde:
        
        Wtotal      ->  é a massa do tank completo
        a           ->  é o raio principal do tank
        r_b         ->  é a altura da sessão elíptica do tank
        K_ratio     ->  é a relação entre a e r_b
        h_lc        ->  é o comprimento da sessão cilíndrica do tank
        h_total     ->  é o comprimento total do tank
        tcilindro   ->  é a espessura da sessão cilíndrica do tank
        t_s         ->  é a espessura da sessão esférica do tank
        t_e         ->  é a espessura da sessão elíptica do tank"""

    Sw = Sw / Fs
    Ve = Vst(r)
    if Ve >= Vp:
        a = calcular_raio_esfera(Vp)
        tsphe = ts(P_int, rho, Sw)
        Wsphe = Ws(r, tsphe, rho)
        h_total = 2 * a + 2 * tsphe

        return [Wsphe, a, 0, 0, 0, h_total, 0, tsphe, 0]
    else:
        r_b = b(Vp, r)
        K_ratio = K(r, r_b)
        tcilindro = tc(P_int, r, Sw)
        a = r - tcilindro
        R = K_ratio * a
        t_k = tk(P_int, K, a, Sw)
        t_cr = tcr(P_int, R, Sw)
        t_e = tee(t_k, t_cr)
        t_s = tse(P_int, a, K, Sw)
        t_c = tc(P_int, a, Sw)
        Vee = V_ee(a, b)
        Vse = V_se(a)
        h_lc = lc(Vp, a, Vee, Vse, t_c)
        Wcilindro = Wc(a, h_lc, t_c, rho)
        Wsphe = Wse(a, t_s, rho)
        Wellip = Wee(a, t_e, E, rho, K)
        Wtotal = Wcilindro + Wsphe + Wellip
        h_total = a + r_b + h_lc + t_s + t_e

        return [Wtotal, a, r_b, K_ratio, h_lc, h_total, tcilindro, t_s, t_e]
