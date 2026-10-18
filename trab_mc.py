import math
import random as rd
import numpy as np
import matplotlib.pyplot as plt


#Modelo do tanque 
def f(r, h):
    return 150 * math.pi * r**2 + 100 * 2 * math.pi * r * h + 200 * math.pi * r**2

def g(r, h):
    return 500 - math.pi * r**2 * h

def P(r, h, mu):
    viol = max(0.0, g(r, h))
    return f(r, h) + mu * viol * viol

def gradiente(r, h, mu, delta=1e-5):
    """Derivadas parciais por diferenças finitas adiantadas (eqs. 4 e 5)."""
    p0 = P(r, h, mu)
    return (P(r + delta, h, mu) - p0) / delta, (P(r, h + delta, mu) - p0) / delta


#Método do gradiente (eqs. 6 e 7): mu e alpha fixos

def metodo_gradiente(mu, alpha, r0, h0, tol=1.0, max_iter=3000000):
    r, h = r0, h0
    hist = []
    for k in range(max_iter):
        if k % 100 == 0:
            hist.append(P(r, h, mu))
        dr, dh = gradiente(r, h, mu)
        if math.hypot(dr, dh) < tol:                 
            return r, h, k, hist, True
        r, h = r - alpha * dr, h - alpha * dh        
        if not (math.isfinite(r) and math.isfinite(h)) or abs(r) > 1e6 or abs(h) > 1e6:
            return r, h, k, hist, False              
    return r, h, max_iter, hist, False

def resumo(mu, alpha, r0, h0):
    r, h, it, hist, ok = metodo_gradiente(mu, alpha, r0, h0)
    return dict(mu=mu, alpha=alpha, r=r, h=h, it=it, hist=hist, ok=ok,
                custo=f(r, h), V=math.pi * r**2 * h)


#Experimentos

def multiplas_execucoes(n, mu, alpha, tol_viab=0.5):
    res = [resumo(mu, alpha, rd.uniform(0.1, 10), rd.uniform(0.1, 20)) for _ in range(n)]
    for i, s in enumerate(res, 1):
        viavel = s["ok"] and 500 - s["V"] <= tol_viab
        s["viavel"] = viavel
        print(f"Exec {i:02d} | it: {s['it']:6d} | r: {s['r']:.4f} | h: {s['h']:.4f} | "
              f"V: {s['V']:.3f} | Custo: R$ {s['custo']:,.2f} | {'viável' if viavel else 'INVIÁVEL/não conv.'}")
    ok = [s for s in res if s["viavel"]]
    if not ok:
        print("Nenhuma execução viável."); return
    c = np.array([s["custo"] for s in ok])
    best = min(ok, key=lambda s: s["custo"])
    print(f"\nViáveis: {len(ok)}/{n} | Custo mín: R$ {c.min():,.2f} | médio: R$ {c.mean():,.2f} | "
          f"máx: R$ {c.max():,.2f}")
    print(f"Ótimo: r = {best['r']:.4f} m, h = {best['h']:.4f} m, V = {best['V']:.3f} m³")

def estudo(casos, r0, h0, titulo, ax, rotulo):
    print(f"\n{titulo}\n{'mu':>8} {'alpha':>9} {'iter':>7} {'r':>8} {'h':>8} {'V':>9} {'custo':>12}  status")
    for mu, alpha in casos:
        s = resumo(mu, alpha, r0, h0)
        print(f"{mu:8.0e} {alpha:9.1e} {s['it']:7d} {s['r']:8.4f} {s['h']:8.4f} {s['V']:9.3f} "
              f"{s['custo']:12,.2f}  {'convergiu' if s['ok'] else 'NÃO convergiu/divergiu'}")
        if all(math.isfinite(v) for v in s["hist"]):
            ax.plot(np.arange(len(s["hist"])) * 100, s["hist"], label=rotulo(mu, alpha))
    ax.set(xscale="symlog", yscale="log", xlabel="Iteração", ylabel="P(r, h)", title=titulo)
    ax.grid(True, ls=":"); ax.legend()

def main():
    try:
        n = int(input("Digite o número de execuções: "))
    except ValueError:
        n = 10; print("Entrada inválida. Usando 10 execuções.")
    MU, ALPHA = 100.0, 8e-8     

    print("=" * 70 + "\n MÚLTIPLAS EXECUÇÕES (chutes iniciais aleatórios) \n" + "=" * 70)
    multiplas_execucoes(n, MU, ALPHA)

    fig, axs = plt.subplots(1, 2, figsize=(14, 5))
    estudo([(m, 8e-6 / m) for m in (10, 30, 100, 300)], 5.0, 10.0,
           "Efeito de mu (alpha ~ 1/mu)", axs[0], lambda m, a: f"mu={m:.0e}, alpha={a:.1e}")
    estudo([(MU, a) for a in (3e-8, 8e-8, 1.1e-7, 3e-7)], 5.0, 10.0,
           f"Efeito de alpha (mu={MU:.0f})", axs[1], lambda m, a: f"alpha={a:.1e}")
    plt.tight_layout(); plt.show()

if __name__ == "__main__":
    main()
