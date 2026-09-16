from collections import defaultdict
from re import S
from typing import Dict, List, Tuple, Optional
import pyphysics as phys
import uncertainties as un
import pandas as pd
import matplotlib.pyplot as plt
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

# Cada nivel experimental: (Ex, l, j, C2S)
# - Si c2s > 0  -> se dibuja como barra coloreada según l (igual criterio que el modelo),
#   normalizada al número de "huecos" disponibles en ese orbital (l, j) — ver
#   VACANCY_EXCEPTIONS / get_vacancies más abajo.
# - Si c2s <= 0 -> se dibuja la barra/línea negra fina de siempre, indicando el nivel
#   sin SF asociado (en ese caso l y j son irrelevantes, pon lo que quieras, p.ej. -1, -1).
ExpLevel = Tuple[float, int, float, float]

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
        tup[1].Ex #type: ignore
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
    left_padding = 0.05
    right_padding = 0.05
    annotated_positions = []

    for i, data in enumerate(models):
        for qstate, vals in data.items():
            for j, tup in enumerate(vals):
                qnucleon, smdata = tup

                if qnucleon.l in exclude_l:
                    continue

                ex = un.nominal_value(smdata.Ex)
                c2s_theo = un.nominal_value(smdata.SF)

                # Máximo número de huecos/vacantes (respeta VACANCY_EXCEPTIONS)
                max_c2s = get_vacancies(qnucleon.l, qnucleon.j)

                # Barra teoría horizontal de 1.5 a 2.0 (ancho total 0.5)
                left = 1.5
                full_width = 0.5
                width = (c2s_theo / max_c2s) * full_width

                # Color base por ℓ
                color = l_colors.get(qnucleon.l, "gray")

                # Oscurecer p1/2
                if qnucleon.l == 1 and abs(qnucleon.j - 0.5) < 1e-6:
                    import matplotlib.colors as mcolors
                    r, g, b = mcolors.to_rgb(color)
                    factor = 0.6
                    color = (r * factor, g * factor, b * factor)

                # Background barra (representa el 100% del limite max_c2s)
                ax.barh(ex, left=left, width=full_width, height=height, color=color, alpha=0.35, edgecolor="none")

                # Foreground barra (proporcional a C2S / max_c2s)
                ax.barh(ex, left=left, width=width, height=height, color=color, alpha=0.75, edgecolor="none")

                # Texto nlj, SF, Jπ
                n = qnucleon.n
                l = qnucleon.l
                j_orb = qnucleon.j
                j_txt = format_j_fraction(j_orb)
                nlj = f"{n}{l_chars[l]}{j_txt}"
                sf_txt = f"C²S={c2s_theo:.2f}"
                pi_char = "⁺" if qstate.pi > 0 else "⁻"
                jpi_txt = f"{format_j_fraction(qstate.J)}{pi_char}"

                text = f"{nlj}   {sf_txt}   {jpi_txt}"

                fontsize_text = 12
                step = 0.16
                offset = 0.0
                while any(abs(ex + offset - pos) < step for pos in annotated_positions):
                    offset += step
                annotated_positions.append(ex + offset)

                ax.annotate(text, xy=(left - right_padding, ex + offset), ha="right", va="center", fontsize=fontsize_text)

    # -------------------------
    # Barras experimentales a la izquierda
    # -------------------------
    # Altura de la barra/línea "sin SF" (negra), un poco más fina que antes
    neg_height = height * 0.55  # antes era height (0.1); ahora más fina

    exp_left = 0.1
    exp_full_width = 0.5

    if exp_data is not None:
        for ex_val, l, j, c2s in exp_data:
            if c2s is not None and c2s > 0:
                color = l_to_color(l)
                max_c2s = get_vacancies(l, j)
                frac = c2s / max_c2s

                # Fondo (referencia del 100% = huecos disponibles en ese orbital)
                ax.barh(
                    ex_val, left=exp_left, width=exp_full_width, height=height,
                    color=color, alpha=0.35, edgecolor="none",
                )
                # Barra proporcional al C2S respecto a las vacancies del orbital
                ax.barh(
                    ex_val, left=exp_left, width=frac * exp_full_width, height=height,
                    color=color, alpha=0.85, edgecolor="none",
                )

                l_char = l_chars[l] if 0 <= l < len(l_chars) else str(l)
                orbital = f"{l_char}{format_j_fraction(j)}"
                # Experimental C2S label
                label = f"C²S={c2s:.2f}" 
                ax.annotate( label, xy=(exp_left + exp_full_width + right_padding, ex_val), ha="left", va="center", fontsize=12, )
            else:
                # Sin C2S asignado (o negativo): barra/línea negra fina, como antes
                ax.barh(
                    ex_val, left=exp_left, width=exp_full_width, height=neg_height,
                    color="black", edgecolor="black",
                )

    # Quitar eje X
    ax.set_xticks([])

    # Limites
    exp_ex_vals = [ex_val for ex_val, _, _, _ in exp_data] if exp_data else [0]
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
    ax.set_xlim(0, 2.3)
    ax.set_ylim(-0.3, max_ex + 0.5)

    # Etiquetas debajo de cada grupo
    ax.text(0.35, -0.5, "Experiment", ha="center", va="top", fontsize=14, fontweight="bold")
    ax.text(1.75, -0.5, "SFO-tls model", ha="center", va="top", fontsize=14, fontweight="bold")


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

fig, ax = plt.subplots(1, 1, figsize=(9, 12))

# Ejemplo de exp_data: (Ex, l, j, C2S)
# - Si C2S > 0 -> barra coloreada según l, normalizada a get_vacancies(l, j)
# - Si C2S <= 0 -> línea negra fina (l, j no importan en ese caso)
exp_data: List[ExpLevel] = [
    (0.0,   1, 1.5,  1.22),      # sin dato -> línea negra fina
    (0.936,  1, 1.5, 0.65),    # p3/2 (l=1, j=3/2): usa la excepción -> máx=2
    (2.19,   1, 1.5, 0.28),    # f7/2 (l=3, j=7/2): sin excepción -> máx=2j+1=8
    (3.29,  -1, -1,  -1),
    (3.81,  -1, -1,  -1),
    (5.21,   1, -1, -1),    # p1/2 (l=1, j=1/2): sin excepción -> máx=2j+1=2
    (5.765, -1, -1,  -1),
    (6.24,  -1, -1,  -1),
]

plot_bars([extracted], ax, exp_data=exp_data, exclude_l=[0])  # Ejemplo excluyendo s

ax.set_ylabel(r"$E_{x}$ [MeV]")
ax.set_title("SF por estado, coloreado por ℓ")

plt.savefig("./SF_bars_Li8.png")
plt.tight_layout()
plt.show()