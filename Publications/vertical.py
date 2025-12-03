from collections import defaultdict
from re import S
from typing import Dict, List, Tuple
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
        tup[1].Ex
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


# -------------------------
# PLOT COMPLETO
# -------------------------
def plot_bars(models: List[OneDict], ax, exp_ex: List[float] = None, **kwargs) -> None:
    height = 0.1  # grosor vertical
    left_padding = 0.05
    right_padding = 0.05

    annotated_positions = []

    # Barras de teoría a la derecha
    for i, data in enumerate(models):
        for qstate, vals in data.items():
            for j, tup in enumerate(vals):
                qnucleon, smdata = tup
                ex = un.nominal_value(smdata.Ex)
                sf = un.nominal_value(smdata.SF) * qnucleon.degeneracy()
                max_sf = qnucleon.degeneracy()

                # Barra teoría horizontal de 1.5 a 2.0
                left = 1.5
                width = (sf / max_sf) * 0.5
                color = l_colors.get(qnucleon.l, "gray")

                # Background barra
                ax.barh(
                    ex,
                    left=left,
                    width=0.5,
                    height=height,
                    color=color,
                    alpha=0.35,
                    edgecolor="none",
                )
                # Foreground barra proporcional a SF
                ax.barh(
                    ex,
                    left=left,
                    width=width,
                    height=height,
                    color=color,
                    alpha=0.75,
                    edgecolor="none",
                )

                # Texto nlj, SF, Jπ
                n = qnucleon.n
                l = qnucleon.l
                j_orb = qnucleon.j
                l_chars = ["s", "p", "d", "f", "g"]
                nlj = f"{n}{l_chars[l]}{j_orb:.1f}".replace(".0", "")
                sf_txt = f"SF={sf / qnucleon.degeneracy():.2f}"
                pi = "+" if qstate.pi > 0 else "-"
                jpi_txt = rf"$J^\pi={qstate.J}^{ {pi} }$"
                text = f"{nlj}   {sf_txt}   {jpi_txt}"

                fontsize_text = 12
                step = 0.16  # separación mínima adaptada al tamaño de fuente

                offset = 0.0
                while any(abs(ex + offset - pos) < step for pos in annotated_positions):
                    offset += step
                annotated_positions.append(ex + offset)

                ax.annotate(
                    text,
                    xy=(left - right_padding, ex + offset),
                    ha="right",
                    va="center",
                    fontsize=fontsize_text
                )

    # Barras experimentales a la izquierda
    if exp_ex is not None:
        for ex_val in exp_ex:
            ax.barh(
                ex_val,
                left=0.1,
                width=0.5,
                height=0.1,
                color="black",
                edgecolor="black",
            )

    # Quitar eje X
    ax.set_xticks([])

    # Limites
    max_ex = max(
        [
            un.nominal_value(tup[1].Ex)
            for data in models
            for vals in data.values()
            for tup in vals
        ]
        + (exp_ex if exp_ex else [0])
    )
    ax.set_xlim(0, 2.3)
    ax.set_ylim(-0.3, max_ex + 0.5)

    # Etiquetas debajo de cada grupo (en coordenadas de datos)
    ax.text(0.35, -0.5, "Experiment", ha="center", va="top", fontsize=14, fontweight="bold")
    ax.text(1.75, -0.5, "Theory", ha="center", va="top", fontsize=14, fontweight="bold")



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
exp_ex = [0.0, 0.98, 2.1, 3.01654, 5.21178, 5.90033, 6.34973]  # posiciones experimentales
plot_bars([extracted], ax, exp_ex=exp_ex)

ax.set_ylabel(r"$E_{x}$ [MeV]")
ax.set_title("SF por estado, coloreado por ℓ")

plt.tight_layout()
plt.show()
