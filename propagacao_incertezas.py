#!/usr/bin/env python3
"""Frequência de oscilação de um fio esticado com propagação de incertezas.

Modelo:  f = (n / 2L) * sqrt(T / mu),  T = m_obj * g,  mu = m_fio / l_fio

Hipóteses: fio ideal, sem atrito na polia, tensão = peso do objeto,
n exato e grandezas medidas independentes (sem correlação).

Interface gráfica (padrão):  python propagacao_incertezas.py
Terminal interativo:         python propagacao_incertezas.py --cli
Uso por argumentos (unidades: L, l_fio em cm; massas em g):
    python propagacao_incertezas.py --n 3 --L 80 --sL 0.2 --m-obj 200 --sm-obj 1 \
        --m-fio 1.2 --sm-fio 0.05 --l-fio 200 --sl-fio 0.5 --g 9.78 --mc
"""
from __future__ import annotations

import argparse
import math
import random
import statistics
from dataclasses import dataclass

UNIDADES_COMP = {"m": 1.0, "cm": 1e-2, "mm": 1e-3}
UNIDADES_MASSA = {"kg": 1.0, "g": 1e-3, "mg": 1e-6}


@dataclass
class Entrada:
    """Todas as grandezas em unidades SI (m, kg, m/s²)."""
    n: float
    L: float
    sL: float
    m_obj: float
    sm_obj: float
    m_fio: float
    sm_fio: float
    l_fio: float
    sl_fio: float
    g: float = 9.78      # g local aproximado; ajuste para o seu laboratório
    sg: float = 0.0

    def validar(self) -> None:
        for nome in ("n", "L", "m_obj", "m_fio", "l_fio", "g"):
            if not getattr(self, nome) > 0:
                raise ValueError(f"'{nome}' deve ser > 0.")
        for nome in ("sL", "sm_obj", "sm_fio", "sl_fio", "sg"):
            if getattr(self, nome) < 0:
                raise ValueError(f"'{nome}' deve ser >= 0.")


@dataclass
class Resultado:
    f: float
    sigma_f: float
    T: float
    sigma_T: float
    mu: float
    sigma_mu: float
    contribuicoes: dict[str, float]   # fração da variância de f (soma = 1)


def calcular(e: Entrada) -> Resultado:
    """Propagação linear (derivadas parciais) para f(L, T, mu)."""
    e.validar()
    T = e.m_obj * e.g
    sigma_T = T * math.hypot(e.sm_obj / e.m_obj, e.sg / e.g)
    mu = e.m_fio / e.l_fio
    sigma_mu = mu * math.hypot(e.sm_fio / e.m_fio, e.sl_fio / e.l_fio)

    f = e.n / (2 * e.L) * math.sqrt(T / mu)

    df_dL = -f / e.L
    df_dT = f / (2 * T)
    df_dmu = -f / (2 * mu)

    termos = {
        "L": (df_dL * e.sL) ** 2,
        "T (m_obj e g)": (df_dT * sigma_T) ** 2,
        "mu (m_fio e l_fio)": (df_dmu * sigma_mu) ** 2,
    }
    var = sum(termos.values())
    contrib = {k: (v / var if var > 0 else 0.0) for k, v in termos.items()}
    return Resultado(f, math.sqrt(var), T, sigma_T, mu, sigma_mu, contrib)


def monte_carlo(e: Entrada, n_amostras: int = 100_000, semente: int | None = 42):
    """Verificação independente: sorteia as grandezas (normais) e mede o desvio de f.
    Usa só a biblioteca padrão (sem numpy)."""
    rng = random.Random(semente)
    fs = []
    for _ in range(n_amostras):
        L = rng.gauss(e.L, e.sL)
        m_obj = rng.gauss(e.m_obj, e.sm_obj)
        m_fio = rng.gauss(e.m_fio, e.sm_fio)
        l_fio = rng.gauss(e.l_fio, e.sl_fio)
        g = rng.gauss(e.g, e.sg)
        fs.append(e.n / (2 * L) * math.sqrt((m_obj * g) / (m_fio / l_fio)))
    return statistics.fmean(fs), statistics.stdev(fs)


def arredonda(valor: float, sigma: float, alg: int = 2) -> tuple[str, str]:
    """σ com `alg` algarismos significativos; valor no mesmo decimal."""
    if sigma <= 0 or not math.isfinite(sigma):
        return f"{valor:.4g}", "0"
    casas = max(0, alg - 1 - math.floor(math.log10(sigma)))
    return f"{valor:.{casas}f}", f"{sigma:.{casas}f}"


def erro_normalizado(f: float, sf: float, f_exp: float, sf_exp: float = 0.0) -> float:
    return abs(f - f_exp) / math.hypot(sf, sf_exp)


def relatorio(e: Entrada, r: Resultado, mc: tuple[float, float] | None = None,
              f_exp: float | None = None, sf_exp: float = 0.0) -> str:
    v, s = arredonda(r.f, r.sigma_f)
    linhas = [
        f"Frequência f = ({v} ± {s}) Hz   [{100 * r.sigma_f / r.f:.2f} %]",
        f"Tensão T = {r.T:.4f} ± {r.sigma_T:.4f} N",
        f"Densidade linear μ = {r.mu * 1e3:.4f} ± {r.sigma_mu * 1e3:.4f} g/m",
        "",
        "Contribuição para a variância de f:",
    ]
    for nome, frac in sorted(r.contribuicoes.items(), key=lambda kv: -kv[1]):
        linhas.append(f"  {nome:<22} {100 * frac:6.1f} %")
    if mc:
        linhas += ["", f"Monte Carlo: f = {mc[0]:.4f} Hz, σ_f = {mc[1]:.4f} Hz "
                       f"(analítico: {r.sigma_f:.4f} Hz)"]
    if f_exp is not None:
        en = erro_normalizado(r.f, r.sigma_f, f_exp, sf_exp)
        veredito = "compatível" if en <= 1 else "incompatível"
        linhas += ["", f"Comparação com f_exp = {f_exp} Hz: E_n = {en:.2f} → {veredito} (|E_n| ≤ 1)"]
    return "\n".join(linhas)


def _ler(prompt: str) -> float:
    while True:
        try:
            return float(input(prompt).strip().replace(",", "."))
        except ValueError:
            print("  Valor inválido, tente de novo.")


def _ler_unid(nome: str, tabela: dict[str, float], padrao: str) -> float:
    while True:
        u = input(f"Unidade de {nome} {list(tabela)} [{padrao}]: ").strip() or padrao
        if u in tabela:
            return tabela[u]
        print("  Unidade desconhecida.")


def modo_interativo() -> Entrada:
    print("Preencha as medidas (vírgula ou ponto decimal).\n")
    n = _ler("n (nº de ventres): ")
    uc = _ler_unid("comprimento (L e l_fio)", UNIDADES_COMP, "cm")
    um = _ler_unid("massa (m_obj e m_fio)", UNIDADES_MASSA, "g")
    L, sL = _ler("L: ") * uc, _ler("σ_L: ") * uc
    mo, smo = _ler("m_obj: ") * um, _ler("σ_m_obj: ") * um
    mf, smf = _ler("m_fio: ") * um, _ler("σ_m_fio: ") * um
    lf, slf = _ler("l_fio: ") * uc, _ler("σ_l_fio: ") * uc
    g = _ler("g local (m/s²) [ex.: 9.78]: ")
    return Entrada(n, L, sL, mo, smo, mf, smf, lf, slf, g)


def abrir_gui() -> None:
    """Interface gráfica (Tkinter, já vem com o Python)."""
    import tkinter as tk
    from tkinter import ttk

    raiz = tk.Tk()
    raiz.title("Propagação de incertezas · fio esticado")
    raiz.minsize(880, 560)
    est = ttk.Style(raiz)
    if "clam" in est.theme_names():
        est.theme_use("clam")
    est.configure("Err.TEntry", fieldbackground="#fde2e2")
    est.configure("Horizontal.TProgressbar", background="#2f5bea", troughcolor="#e5e7eb",
                  bordercolor="#c9ced6", lightcolor="#2f5bea", darkcolor="#2f5bea")

    f_esq = ttk.Frame(raiz, padding=14)
    f_dir = ttk.Frame(raiz, padding=14)
    f_esq.grid(row=0, column=0, sticky="n")
    f_dir.grid(row=0, column=1, sticky="nsew")
    raiz.columnconfigure(1, weight=1)
    raiz.rowconfigure(0, weight=1)

    ttk.Label(f_esq, text="f = (n / 2L) · √(T / μ)", font=("TkDefaultFont", 13, "bold")
              ).grid(row=0, column=0, columnspan=4, sticky="w")
    ttk.Label(f_esq, text="T = m_obj·g   μ = m_fio / l_fio", foreground="#666"
              ).grid(row=1, column=0, columnspan=4, sticky="w", pady=(0, 10))
    for c, t in enumerate(("Grandeza", "Valor", "Incerteza σ", "Unid.")):
        ttk.Label(f_esq, text=t, foreground="#666").grid(row=2, column=c, sticky="w", padx=3)

    # chave, rótulo, valor padrão, σ padrão, tabela de unidades, unidade padrão
    campos = [
        ("n", "n (nº de ventres)", "3", None, None, None),
        ("L", "L (comp. vibrante)", "80", "0,2", UNIDADES_COMP, "cm"),
        ("mo", "m_obj (massa suspensa)", "200", "1", UNIDADES_MASSA, "g"),
        ("mf", "m_fio (massa do fio)", "1,2", "0,05", UNIDADES_MASSA, "g"),
        ("lf", "l_fio (comp. total)", "200", "0,5", UNIDADES_COMP, "cm"),
        ("g", "g local (m/s²)", "9,78", "0", None, None),
    ]
    val, sig, uni, ents = {}, {}, {}, {}
    for i, (k, rot, v0, s0, tab, u0) in enumerate(campos, start=3):
        ttk.Label(f_esq, text=rot).grid(row=i, column=0, sticky="w", padx=3, pady=3)
        val[k] = tk.StringVar(value=v0)
        ents[k] = ttk.Entry(f_esq, textvariable=val[k], width=11)
        ents[k].grid(row=i, column=1, padx=3)
        if s0 is not None:
            sig[k] = tk.StringVar(value=s0)
            ents["s" + k] = ttk.Entry(f_esq, textvariable=sig[k], width=11)
            ents["s" + k].grid(row=i, column=2, padx=3)
        if tab:
            uni[k] = tk.StringVar(value=u0)
            cb = ttk.Combobox(f_esq, textvariable=uni[k], values=list(tab), width=5, state="readonly")
            cb.grid(row=i, column=3, padx=3)
            cb.bind("<<ComboboxSelected>>", lambda _e: atualizar())

    ttk.Separator(f_esq).grid(row=9, column=0, columnspan=4, sticky="ew", pady=10)
    ttk.Label(f_esq, text="Comparar com valor medido (opcional)").grid(row=10, column=0, columnspan=4, sticky="w")
    fe, sfe = tk.StringVar(), tk.StringVar()
    ttk.Label(f_esq, text="f_exp / σ (Hz)").grid(row=11, column=0, sticky="w", padx=3, pady=3)
    ttk.Entry(f_esq, textvariable=fe, width=11).grid(row=11, column=1, padx=3)
    ttk.Entry(f_esq, textvariable=sfe, width=11).grid(row=11, column=2, padx=3)

    ttk.Label(f_dir, text="Resultado", foreground="#666").pack(anchor="w")
    big = ttk.Label(f_dir, text="", font=("TkDefaultFont", 22, "bold"))
    big.pack(anchor="w", pady=(2, 6))
    corpo = ttk.Label(f_dir, text="", font=("TkFixedFont", 10), justify="left")
    corpo.pack(anchor="w")
    ttk.Label(f_dir, text="Contribuição para σ_f (variância)", foreground="#666"
              ).pack(anchor="w", pady=(14, 4))
    barras = {}
    for nome in ("L", "T (m_obj e g)", "mu (m_fio e l_fio)"):
        linha = ttk.Frame(f_dir)
        linha.pack(fill="x", pady=2)
        lab = ttk.Label(linha, text=nome, width=20)
        lab.pack(side="left")
        pb = ttk.Progressbar(linha, maximum=100, length=220)
        pb.pack(side="left", padx=6)
        pct = ttk.Label(linha, text="", width=8)
        pct.pack(side="left")
        barras[nome] = (pb, pct)
    dica = ttk.Label(f_dir, text="", wraplength=380, foreground="#1d4ed8")
    dica.pack(anchor="w", pady=(10, 0))

    estado = {"texto": "", "mc": None}

    def resumo(texto: str) -> str:
        """Remove do relatório a 1ª linha e a lista de contribuições (já mostradas à parte)."""
        l = texto.split("\n")
        return "\n".join(l[1:3] + l[8:]).strip()

    def ler(k: str) -> float:
        return float(ents[k].get().strip().replace(",", ".")) if k in ents else 0.0

    def atualizar(*_):
        estado["mc"] = None
        try:
            uc = UNIDADES_COMP[uni["L"].get()]
            um = UNIDADES_MASSA[uni["mo"].get()]
            umf = UNIDADES_MASSA[uni["mf"].get()]
            ul = UNIDADES_COMP[uni["lf"].get()]
            e = Entrada(
                n=float(val["n"].get().replace(",", ".")),
                L=float(val["L"].get().replace(",", ".")) * uc,
                sL=float(sig["L"].get().replace(",", ".")) * uc,
                m_obj=float(val["mo"].get().replace(",", ".")) * um,
                sm_obj=float(sig["mo"].get().replace(",", ".")) * um,
                m_fio=float(val["mf"].get().replace(",", ".")) * umf,
                sm_fio=float(sig["mf"].get().replace(",", ".")) * umf,
                l_fio=float(val["lf"].get().replace(",", ".")) * ul,
                sl_fio=float(sig["lf"].get().replace(",", ".")) * ul,
                g=float(val["g"].get().replace(",", ".")),
                sg=float(sig["g"].get().replace(",", ".")),
            )
            r = calcular(e)
            fexp = float(fe.get().replace(",", ".")) if fe.get().strip() else None
            sfexp = float(sfe.get().replace(",", ".")) if sfe.get().strip() else 0.0
        except (ValueError, KeyError) as erro:
            big.config(text="—")
            corpo.config(text=f"Entrada inválida: {erro}" if str(erro) else "Preencha todos os campos com números.")
            for pb, pct in barras.values():
                pb["value"] = 0
                pct.config(text="")
            dica.config(text="")
            estado["texto"] = ""
            return
        estado.update(e=e, r=r, fexp=fexp, sfexp=sfexp)
        v, s = arredonda(r.f, r.sigma_f)
        big.config(text=f"{v} ± {s} Hz")
        texto = relatorio(e, r, estado["mc"], fexp, sfexp)
        corpo.config(text=resumo(texto))
        estado["texto"] = texto
        for nome, (pb, pct) in barras.items():
            frac = r.contribuicoes[nome]
            pb["value"] = 100 * frac
            pct.config(text=f"{100 * frac:.1f} %")
        top = max(r.contribuicoes, key=r.contribuicoes.get)
        dica.config(text=f"Para reduzir σ_f, melhore primeiro a medida de: {top}.")

    def rodar_mc():
        if "e" in estado and estado["texto"]:
            estado["mc"] = monte_carlo(estado["e"])
            texto = relatorio(estado["e"], estado["r"], estado["mc"], estado["fexp"], estado["sfexp"])
            estado["texto"] = texto
            corpo.config(text=resumo(texto))

    def copiar():
        if estado["texto"]:
            raiz.clipboard_clear()
            raiz.clipboard_append(estado["texto"])

    def exemplo():
        for k, (v0, s0) in {"n": ("3", None), "L": ("80", "0,2"), "mo": ("200", "1"),
                            "mf": ("1,2", "0,05"), "lf": ("200", "0,5"), "g": ("9,78", "0")}.items():
            val[k].set(v0)
            if s0 is not None:
                sig[k].set(s0)
        fe.set("")
        sfe.set("")
        atualizar()

    bt = ttk.Frame(f_esq)
    bt.grid(row=12, column=0, columnspan=4, sticky="w", pady=(14, 0))
    ttk.Button(bt, text="Exemplo", command=exemplo).pack(side="left", padx=3)
    ttk.Button(bt, text="Monte Carlo", command=rodar_mc).pack(side="left", padx=3)
    ttk.Button(bt, text="Copiar resultado", command=copiar).pack(side="left", padx=3)

    for var in list(val.values()) + list(sig.values()) + [fe, sfe]:
        var.trace_add("write", atualizar)
    atualizar()
    raiz.mainloop()


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for nome in ("n", "L", "sL", "m-obj", "sm-obj", "m-fio", "sm-fio", "l-fio", "sl-fio"):
        p.add_argument(f"--{nome}", type=float)
    p.add_argument("--g", type=float, default=9.78)
    p.add_argument("--sg", type=float, default=0.0)
    p.add_argument("--uc", choices=UNIDADES_COMP, default="cm", help="unidade de L e l_fio")
    p.add_argument("--um", choices=UNIDADES_MASSA, default="g", help="unidade das massas")
    p.add_argument("--f-exp", type=float, help="frequência medida (Hz)")
    p.add_argument("--sf-exp", type=float, default=0.0)
    p.add_argument("--mc", action="store_true", help="verifica com Monte Carlo")
    p.add_argument("--cli", action="store_true", help="usa o terminal em vez da janela")
    a = p.parse_args()

    obrig = [a.n, a.L, a.sL, a.m_obj, a.sm_obj, a.m_fio, a.sm_fio, a.l_fio, a.sl_fio]
    if all(x is None for x in obrig) and not a.cli:
        try:
            abrir_gui()
            return
        except Exception as erro:  # sem tkinter/display: cai para o terminal
            print(f"(Interface gráfica indisponível: {erro}. Usando o terminal.)\n")
    if any(x is None for x in obrig):
        e = modo_interativo()
    else:
        c, m = UNIDADES_COMP[a.uc], UNIDADES_MASSA[a.um]
        e = Entrada(a.n, a.L * c, a.sL * c, a.m_obj * m, a.sm_obj * m,
                    a.m_fio * m, a.sm_fio * m, a.l_fio * c, a.sl_fio * c, a.g, a.sg)
    try:
        r = calcular(e)
    except ValueError as erro:
        raise SystemExit(f"Erro: {erro}")
    print()
    print(relatorio(e, r, monte_carlo(e) if a.mc else None, a.f_exp, a.sf_exp))


if __name__ == "__main__":
    main()
