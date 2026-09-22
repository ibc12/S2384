from collections import defaultdict
from typing import Dict, List, Tuple, Optional
import pyphysics as phys
import uncertainties as un
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Patch
import copy
import re
import matplotlib as mpl

mpl.rcParams["font.family"] = "sans-serif"


class QState:
    def __init__(self, J: float, pi: int) -> None:
        self.J: float = J
        self.pi: int = pi

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, QState):
            return NotImplemented
        return self.J == other.J and self.pi == other.pi

    def __hash__(self) -> int:
        return hash((self.J, self.pi))

    def __str__(self) -> str:
        return f"QState({self.J}, {self.pi})"

    def __repr__(self) -> str:
        return f"QState({self.J}, {self.pi})"


TwoDict = Dict[QState, Dict[float, List]]
OneDict = Dict[QState, List[Tuple[phys.QuantumNumbers, phys.ShellModelData]]]

# Cada nivel experimental: (Ex, l, j, C2S, C2S_err, Jpi)
# - Si c2s > 0  -> se dibuja como barra coloreada según l (igual criterio que el modelo),
#   normalizada al número de "huecos" disponibles en ese orbital (l, j) — ver
#   VACANCY_EXCEPTIONS / get_vacancies más abajo.
# - Si c2s <= 0 -> se dibuja la barra/línea negra fina de siempre, indicando el nivel
#   sin SF asociado (en ese caso l y j son irrelevantes, pon lo que quieras, p.ej. -1, -1).
# - C2S_err es la incertidumbre absoluta del C2S (misma escala que C2S), o None si no hay.
#   Se escribe en notación compacta tras el valor, con 2 cifras significativas en el error
#   y el C2S redondeado al mismo decimal:
#     (1.22, 0.20)  -> "C²S=1.22(20)"
#     (0.28, 0.056) -> "C²S=0.280(56)"
# - Jpi es un string libre (p.ej. "3/2⁻", "(1,2)⁺") o None si no quieres escribir nada.
#   Se dibuja a la izquierda de la barra, junto con el orbital (si hay C2S).
ExpLevel = Tuple[float, int, float, float, Optional[float], Optional[str]]

# -------------------------
# VACANCIAS (tamaño total de la barra experimental)
# -------------------------
# Por defecto el tamaño total de un orbital (l, j) es su degeneración 2j+1.
# Pero si en el núcleo de referencia ese orbital ya está parcialmente ocupado,
# el número de huecos (vacancies) disponibles para la transferencia es menor.
# Añade aquí las excepciones concretas de tu núcleo/reacción; el resto de
# orbitales usarán 2j+1 automáticamente.
#
# Ejemplo: en 7Li el orbital p3/2 (l=1, j=3/2) tiene degeneración 2j+1=4,
# pero ya hay 2 neutrones ahí, así que solo quedan 2 huecos.
VACANCY_EXCEPTIONS: Dict[Tuple[int, float], float] = {
    (1, 1.5): 2.0,  # p3/2 en 7Li: 4 - 2 neutrones ya presentes = 2 huecos
}


def get_vacancies(l: int, j: float) -> float:
    """Tamaño total (huecos disponibles) del orbital (l, j).

    Usa VACANCY_EXCEPTIONS si el orbital está registrado ahí; si no,
    cae por defecto en la degeneración estándar 2j+1.
    """
    return VACANCY_EXCEPTIONS.get((l, j), 2 * j + 1)


def parse_files(files: list) -> TwoDict:
    ret: TwoDict = {}
    for file in files:
        with open(file, "r") as f:
            n, l, j = -1, -1, -1
            pi = +1
            for lin in f:
                line = lin.strip()
                if not line:
                    continue

                if "orbit" in line:
                    for c, column in enumerate(line.split()):
                        if c == 2:
                            n = int(column)
                        elif c == 3:
                            l = int(column)
                        elif c == 4:
                            j = int(column)

                if "parity" in line:
                    toks = line.split()
                    pi_value = int(toks[-2])  # el bueno
                    pi = +1 if pi_value > 0 else -1

                if re.match(r"^\d+\(", line):
                    J = float(line[0]) / 2
                    ex = float(line[34:41].strip())
                    c2s = float(line[45:51].strip())

                    qstate = QState(J, pi)
                    qnucleon = phys.QuantumNumbers(n, l, j / 2)
                    sm = phys.ShellModelData(ex, c2s)

                    if c2s > 0.1:
                        if qstate not in ret:
                            ret[qstate] = {}
                        if ex not in ret[qstate]:
                            ret[qstate][ex] = []
                        ret[qstate][ex].append((qnucleon, sm))

    return ret


def extract_max_sf(data: TwoDict) -> OneDict:
    ret: OneDict = {}
    for q, inner in data.items():
        for ex, vals in inner.items():
            largestSF = max((tup for tup in vals), key=lambda x: x[1].SF)
            if q not in ret:
                ret[q] = []
            ret[q].append(largestSF)
    return ret


def shift_ex(data: TwoDict | OneDict) -> None:
    maxEx = max(
        tup[1].Ex  # type: ignore
        for v in data.values()
        for tup in (v if isinstance(v, list) else sum(v.values(), []))
    )

    for v in data.values():
        if isinstance(v, list):
            for q, sm in v:
                sm.Ex = maxEx - sm.Ex
        else:
            for ex, lst in v.items():
                for q, sm in lst:
                    sm.Ex = maxEx - sm.Ex


# -------------------------
# COLORES POR L
# -------------------------
l_colors = {
    0: "red",  # s
    1: "dodgerblue",  # p
    2: "green",  # d
    3: "orange",  # f
    4: "purple",  # g
}

l_chars = ["s", "p", "d", "f", "g"]

# -------------------------
# CONVERSIÓN j → FRACCIÓN
# -------------------------
def format_j_fraction(j: float) -> str:
    """Convierte j=0.5→1/2, 1.0→1, 1.5→3/2, etc."""
    twice = round(2 * j)
    if twice % 2 == 0:
        return f"{twice//2}"
    else:
        return f"{twice}/2"


def l_to_color(l: int) -> str:
    return l_colors.get(l, "gray")


def _dodge(y: float, used: List[float], step: float) -> float:
    """Sube el texto en pasos de `step` hasta que no pise a otro ya colocado."""
    offset = 0.0
    while any(abs(y + offset - pos) < step for pos in used):
        offset += step
    used.append(y + offset)
    return y + offset


# -------------------------
# PLOT COMPLETO
# -------------------------
def plot_bars(
    models: List[OneDict],
    ax,
    exp_data: Optional[List[ExpLevel]] = None,
    exclude_l: List[int] = None,  # type: ignore
    **kwargs,
) -> None:
    if exclude_l is None:
        exclude_l = []

    height = 0.1
    pad = 0.05
    fontsize_text = 16
    step = 0.16  # separación vertical mínima entre textos

    # Geometría horizontal
    exp_left = 0.1  # barras experimentales: [0.1, 0.6]
    th_left = 1.5  # barras teóricas:        [1.5, 2.0]
    bar_width = 0.5

    # Posiciones y ya ocupadas por texto (una lista por lado para poder repartirlos)
    used_theory: List[float] = []
    used_exp: List[float] = []

    # -------------------------
    # Barras teóricas (derecha)
    #   - C²S en el hueco central (a la izquierda de la barra)
    #   - nlj y Jπ a la derecha de la barra
    # -------------------------
    for data in models:
        for qstate, vals in data.items():
            for qnucleon, smdata in vals:
                if qnucleon.l in exclude_l:
                    continue

                ex = un.nominal_value(smdata.Ex)
                c2s_theo = un.nominal_value(smdata.SF)

                # Máximo número de huecos/vacantes (respeta VACANCY_EXCEPTIONS)
                max_c2s = get_vacancies(qnucleon.l, qnucleon.j)
                width = (c2s_theo / max_c2s) * bar_width

                # Color base por ℓ
                color = l_colors.get(qnucleon.l, "gray")

                # Oscurecer p1/2
                if qnucleon.l == 1 and abs(qnucleon.j - 0.5) < 1e-6:
                    r, g, b = mcolors.to_rgb(color)
                    factor = 0.6
                    color = (r * factor, g * factor, b * factor)

                # Background barra (representa el 100% del limite max_c2s)
                ax.barh(ex, left=th_left, width=bar_width, height=height, color=color, alpha=0.35, edgecolor="none")

                # Foreground barra (proporcional a C2S / max_c2s)
                ax.barh(ex, left=th_left, width=width, height=height, color=color, alpha=0.75, edgecolor="none")

                # Textos
                nlj = f"{qnucleon.n}{l_chars[qnucleon.l]}{format_j_fraction(qnucleon.j)}"
                pi_char = "⁺" if qstate.pi > 0 else "⁻"
                jpi_txt = f"{format_j_fraction(qstate.J)}{pi_char}"

                y_txt = _dodge(ex, used_theory, step)

                # C²S: centro, pegado a la barra teórica
                ax.annotate(
                    #f"C²S={c2s_theo:.2f}",
                    f"{c2s_theo:.2f}",
                    xy=(th_left - pad, y_txt), ha="right", va="center", fontsize=fontsize_text,
                )
                # nlj + Jπ: a la derecha de la barra
                ax.annotate(
                    f"{nlj}   {jpi_txt}",
                    xy=(th_left + bar_width + pad, y_txt), ha="left", va="center", fontsize=fontsize_text,
                )

    # -------------------------
    # Barras experimentales (izquierda)
    #   - orbital (l j) y Jπ a la izquierda de la barra
    #   - C²S en el hueco central (a la derecha de la barra)
    # -------------------------
    neg_height = height * 0.55  # altura de la línea negra "sin SF"

    if exp_data is not None:
        for ex_val, l, j, c2s, c2s_err, jpi in exp_data:
            has_sf = c2s is not None and c2s > 0

            if has_sf:
                color = l_to_color(l)
                frac = c2s / get_vacancies(l, j)

                # Fondo (referencia del 100% = huecos disponibles en ese orbital)
                ax.barh(
                    ex_val, left=exp_left, width=bar_width, height=height,
                    color=color, alpha=0.35, edgecolor="none",
                )
                # Barra proporcional al C2S respecto a las vacancies del orbital
                ax.barh(
                    ex_val, left=exp_left, width=frac * bar_width, height=height,
                    color=color, alpha=0.85, edgecolor="none",
                )
            else:
                # Sin C2S asignado (o negativo): línea negra fina
                ax.barh(
                    ex_val, left=exp_left, width=bar_width, height=neg_height,
                    color="black", edgecolor="black",
                )

            # Texto de la izquierda: orbital (solo si hay C2S) + Jπ (si se ha dado)
            left_parts = []
            if has_sf:
                l_char = l_chars[l] if 0 <= l < len(l_chars) else str(l)
                left_parts.append(f"{l_char}{format_j_fraction(j)}")
            if jpi:
                left_parts.append(jpi)
            left_txt = "   ".join(left_parts)

            right_txt = ""
            if has_sf:
                if c2s_err is not None and c2s_err > 0:
                    # Notación compacta: el error se da con 2 cifras significativas y el
                    # C2S se redondea al mismo decimal -> 0.28 ± 0.056 -> "0.280(56)"
                    #right_txt = f"C²S={un.ufloat(c2s, c2s_err):.2uS}"
                    right_txt = f"{un.ufloat(c2s, c2s_err):.2uS}"
                else:
                    #right_txt = f"C²S={c2s:.2f}"
                    right_txt = f"{c2s:.2f}"

            if left_txt or right_txt:
                y_txt = _dodge(ex_val, used_exp, step)
                if left_txt:
                    ax.annotate(
                        left_txt,
                        xy=(exp_left - pad, y_txt), ha="right", va="center", fontsize=fontsize_text,
                    )
                if right_txt:
                    ax.annotate(
                        right_txt,
                        xy=(exp_left + bar_width + pad, y_txt), ha="left", va="center", fontsize=fontsize_text,
                    )

    # Quitar eje X
    ax.set_xticks([])

    # Limites
    exp_ex_vals = [lvl[0] for lvl in exp_data] if exp_data else [0]
    max_ex = max(
        [
            un.nominal_value(tup[1].Ex)
            for data in models
            for vals in data.values()
            for tup in vals
            if tup[0].l not in exclude_l
        ]
        + exp_ex_vals
    )
    # x negativo: hace falta sitio para los textos a la izquierda de las barras experimentales
    ax.set_xlim(-0.4, 2.5)
    ax.set_ylim(-0.3, max_ex + 0.5)

    # Etiquetas debajo de cada grupo (centradas en sus barras)
    ax.text(exp_left + bar_width / 2, -0.5, "Experiment", ha="center", va="top", fontsize=16, fontweight="bold")
    ax.text(th_left + bar_width / 2, -0.5, "SFO-tls model", ha="center", va="top", fontsize=16, fontweight="bold")


# -------------------------
# USO
# -------------------------
files = [
    "../../Tese/Li7-Li8/log_Li8_sfop6-16a_Li7_sfop6-16a_tr_j4p_j3n.txt",
    "../../Tese/Li7-Li8/log_Li8_sfop6-16a_Li7_sfop6-16a_tr_j2p_j3n.txt",
    "../../Tese/Li7-Li8/log_Li8_sfop6-16a_Li7_sfop6-16a_tr_j6p_j3n.txt",
    "../../Tese/Li7-Li8/log_Li8_sfop6-16a_Li7_sfop6-16a_tr_j4n_j3n.txt",
    "../../Tese/Li7-Li8/log_Li8_sfop6-16a_Li7_sfop6-16a_tr_j8n_j3n.txt",
    "../../Tese/Li7-Li8/log_Li8_sfop6-16a_Li7_sfop6-16a_tr_j2n_j3n.txt",
    "../../Tese/Li7-Li8/log_Li8_sfop6-16a_Li7_sfop6-16a_tr_j0p_j3n.txt",
]

ret = parse_files(files)
extracted = extract_max_sf(ret)
shift_ex(extracted)

plt.rcParams["ytick.labelsize"] = 18
fig, ax = plt.subplots(1, 1, figsize=(10, 12))  # un poco más ancho por los textos laterales

# exp_data: (Ex, l, j, C2S, C2S_err, Jpi)
# - Si C2S > 0 -> barra coloreada según l, normalizada a get_vacancies(l, j)
# - Si C2S <= 0 -> línea negra fina (l, j no importan en ese caso)
# - C2S_err: incertidumbre absoluta del C2S, o None (no se dibuja paréntesis)
# - Jpi: string o None
# OJO: los Jπ de abajo son de ejemplo (2⁺, 1⁺, 3⁺ típicos de 8Li) — revísalos/complétalos.
exp_data: List[ExpLevel] = [
    (0.0,   1, 1.5,  1.22, 0.24,  "2⁺"),   # -> C²S=1.22(20)
    (0.936, 1, 1.5,  0.65, 0.13,  "1⁺"),   # TODO: rellenar errores
    (2.19,  1, 1.5,  0.28, 0.056, "3⁺"),   # -> C²S=0.280(56)
    (3.29, -1, -1,  -1,    None, None),
    (3.81, -1, -1,  -1,    None, None),
    (5.21,  1, -1,  -1,    None, None),
    (5.765, -1, -1, -1,    None, None),
    (6.24, -1, -1,  -1,    None, None),
]

plot_bars([extracted], ax, exp_data=exp_data, exclude_l=[0,2])  # Ejemplo excluyendo s

ax.set_ylabel(r"$E_{x}$ [MeV]", fontsize=18)
ax.set_title("SF por estado, coloreado por ℓ")

plt.tight_layout()  # antes del savefig, si no el guardado ignora el ajuste
plt.savefig("./SF_bars_Li8.png", bbox_inches="tight")
plt.show()