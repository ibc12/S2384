"""
efficiency_clopper_pearson.py
-------------------------------
Lee un CSV con columnas (run, entry, type, status) y calcula la eficiencia
de reconstruccion con intervalos de Clopper-Pearson a 2 sigma (95.45% CL),
igual que en la app de etiquetado.

Se calculan tres eficiencias, TODAS con el mismo denominador
(Binary + Broken + Holes), tal y como hace App._stats():

    eff_binary          = Binary(ok)                    / (Binary+Broken+Holes)
    eff_binary_holes    = (Binary+Holes)(ok)             / (Binary+Broken+Holes)
    eff_binary_broken_holes = (Binary+Broken+Holes)(ok)  / (Binary+Broken+Holes)

Uso:
    python efficiency_clopper_pearson.py mi_archivo.csv
    python efficiency_clopper_pearson.py mi_archivo.csv --cl 0.6827   # 1 sigma
    python efficiency_clopper_pearson.py mi_archivo.csv --by-run
"""

import argparse
import sys
import pandas as pd
from scipy.stats import beta


# ---------------------------------------------------------------------------
# Clopper-Pearson
# ---------------------------------------------------------------------------
def clopper_pearson(k: int, n: int, cl: float = 0.9545):
    """
    Calcula el intervalo de Clopper-Pearson exacto.

    Parametros
    ----------
    k  : numero de exitos (bien reconstruidos)
    n  : numero total de eventos (denominador)
    cl : nivel de confianza (default 0.9545 = 2 sigma)

    Retorna
    -------
    eff        : eficiencia central k/n
    lo, hi     : limites inferior y superior del intervalo
    delta_lo   : incertidumbre inferior  (eff - lo)
    delta_hi   : incertidumbre superior  (hi  - eff)
    """
    if n == 0:
        # No hay denominador -> no se puede calcular nada
        return float("nan"), float("nan"), float("nan"), float("nan"), float("nan")
    if k < 0 or k > n:
        raise ValueError(f"k={k} fuera del rango [0, n={n}].")

    alpha = 1.0 - cl

    lo = 0.0 if k == 0 else beta.ppf(alpha / 2.0, k, n - k + 1)
    hi = 1.0 if k == n else beta.ppf(1.0 - alpha / 2.0, k + 1, n - k)

    eff = k / n
    return eff, lo, hi, eff - lo, hi - eff


# ---------------------------------------------------------------------------
# Lectura y filtrado del CSV
# ---------------------------------------------------------------------------
def load_events(filepath: str) -> pd.DataFrame:
    """
    Lee el CSV completo y normaliza columnas/tipos.
    Devuelve el dataframe con 'type' limpio (strip) y 'status' como booleano.
    """
    df = pd.read_csv(
        filepath,
        dtype={"run": str, "entry": str, "type": str, "status": str},
        skipinitialspace=True,
    )

    df.columns = df.columns.str.strip().str.lower()

    required = {"run", "entry", "type", "status"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Columnas faltantes en el CSV: {missing}")

    df["type"] = df["type"].str.strip()

    df["status"] = df["status"].astype(str).str.strip().str.lower().map(
        {"true": True, "false": False, "1": True, "0": False}
    )

    n_invalid = df["status"].isna().sum()
    if n_invalid > 0:
        print(
            f"[Aviso] {n_invalid} filas con valor de status no reconocido "
            f"seran ignoradas."
        )
        df = df.dropna(subset=["status"])

    return df


def get_subsets(df: pd.DataFrame):
    """
    Construye los subconjuntos usados en las eficiencias, igual que en la app:

        binary            -> type == 'Binary'
        broken            -> type == 'Broken'
        holes             -> type == 'Holes'
        bin_holes         -> Binary + Holes
        bin_broken_holes  -> Binary + Broken + Holes  (denominador comun)
    """
    type_lower = df["type"].str.lower()

    binary = df[type_lower == "binary"]
    broken = df[type_lower == "broken"]
    holes = df[type_lower == "holes"]

    bin_holes = pd.concat([binary, holes])
    bin_broken_holes = pd.concat([binary, broken, holes])

    return binary, broken, holes, bin_holes, bin_broken_holes


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Calcula eficiencias (denominador comun Binary+Broken+Holes) "
        "con intervalos Clopper-Pearson."
    )
    parser.add_argument("csv", help="Ruta al archivo CSV")
    parser.add_argument(
        "--cl",
        type=float,
        default=0.9545,
        help="Nivel de confianza (default: 0.9545 = 2 sigma)",
    )
    parser.add_argument(
        "--by-run",
        action="store_true",
        help="Calcula las eficiencias tambien por run",
    )
    parser.add_argument(
        "--own-denom",
        action="store_true",
        help=(
            "Usa el denominador PROPIO de cada eficiencia "
            "(Binary/Binary, Bin+Holes/Bin+Holes, Bin+Broken+Holes/Bin+Broken+Holes) "
            "en vez del denominador comun (Binary+Broken+Holes) que se usa por defecto."
        ),
    )
    args = parser.parse_args()

    # Booleano que decide que funcion de reporte se usa
    USE_OWN_DENOMINATOR = args.own_denom

    print(f"\n{'#'*60}")
    if USE_OWN_DENOMINATOR:
        print("  MODO: denominador PROPIO por eficiencia")
        print("        (Binary/Binary, Bin+Holes/Bin+Holes, Bin+Brk+Hol/Bin+Brk+Hol)")
    else:
        print("  MODO: denominador COMUN (Binary+Broken+Holes) [por defecto]")
    print(f"{'#'*60}")

    print(f"\nLeyendo: {args.csv}")
    try:
        df = load_events(args.csv)
    except Exception as e:
        print(f"[Error] {e}")
        sys.exit(1)

    cl = args.cl
    sigma_label = {0.6827: "1σ", 0.9545: "2σ", 0.9973: "3σ"}.get(cl, f"{cl*100:.2f}%")

    def report(df_subset: pd.DataFrame, label: str):
        binary, broken, holes, bin_holes, bin_broken_holes = get_subsets(df_subset)

        n_total = len(bin_broken_holes)  # denominador comun
        n_events = len(df_subset)  # total de eventos (todos los tipos)

        ok_binary = int((binary["status"] == True).sum())
        ok_bin_holes = int((bin_holes["status"] == True).sum())
        ok_bin_broken_holes = int((bin_broken_holes["status"] == True).sum())

        e1 = clopper_pearson(ok_binary, n_total, cl)
        e2 = clopper_pearson(ok_bin_holes, n_total, cl)
        e3 = clopper_pearson(ok_bin_broken_holes, n_total, cl)

        print(f"\n{'='*60}")
        print(f"  {label}  —  Clopper-Pearson {sigma_label}")
        print(f"{'='*60}")
        print(f"  Total de eventos : {n_events}")
        print(f"{'-'*60}")

        # --- Desglose por tipo: buenos / malos / % sobre Binary+Broken+Holes ---
        # (se omite el tipo "Other")
        print(f"  {'Tipo':<10} {'Buenos':>8} {'Malos':>8} {'Total':>8} {'% (Bin+Brk+Hol)':>16}")
        print(f"  {'-'*54}")
        type_order = ["Binary", "Multi", "Broken", "Holes"]
        present_types = list(df_subset["type"].dropna().unique())
        # Tipos conocidos primero, en orden; luego cualquier otro no previsto (salvo "Other")
        ordered_types = [t for t in type_order if t in present_types] + [
            t for t in present_types if t not in type_order and t != "Other"
        ]
        for t in ordered_types:
            grp_t = df_subset[df_subset["type"] == t]
            n_t = len(grp_t)
            good_t = int((grp_t["status"] == True).sum())
            bad_t = n_t - good_t
            pct_t = (n_t / n_total * 100) if n_total > 0 else 0.0
            print(f"  {t:<10} {good_t:>8} {bad_t:>8} {n_t:>8} {pct_t:>15.2f}%")
        print(f"  {'-'*54}")
        print(f"  Denominador eficiencias (Binary+Broken+Holes) : {n_total}")
        print(f"{'-'*60}")

        for name, ok, (eff, lo, hi, dm, dp) in [
            ("Binary                / (Bin+Brk+Hol)", ok_binary, e1),
            ("(Binary+Holes)        / (Bin+Brk+Hol)", ok_bin_holes, e2),
            ("(Binary+Broken+Holes) / (Bin+Brk+Hol)", ok_bin_broken_holes, e3),
        ]:
            if n_total == 0:
                print(f"  {name} : sin eventos en el denominador")
                continue
            print(
                f"  {name} : k={ok:>5}  n={n_total:>5}  "
                f"ε = {eff*100:6.2f} %  +{dp*100:5.2f}  -{dm*100:5.2f}  ({sigma_label})"
            )
        print(f"{'='*60}\n")

    def report_own_denominator(df_subset: pd.DataFrame, label: str):
        """
        Version en la que cada eficiencia usa SU PROPIO denominador, en vez del
        denominador comun (Binary+Broken+Holes). Es el equivalente al primer
        bloque de _stats/_efficiency de la app (antes de pasar al denominador
        comun):

            eff_binary               = Binary(ok)               / Binary
            eff_binary_holes         = (Binary+Holes)(ok)        / (Binary+Holes)
            eff_binary_broken_holes  = (Binary+Broken+Holes)(ok) / (Binary+Broken+Holes)
        """
        binary, broken, holes, bin_holes, bin_broken_holes = get_subsets(df_subset)
        n_events = len(df_subset)

        ok_binary = int((binary["status"] == True).sum())
        ok_bin_holes = int((bin_holes["status"] == True).sum())
        ok_bin_broken_holes = int((bin_broken_holes["status"] == True).sum())

        e1 = clopper_pearson(ok_binary, len(binary), cl)
        e2 = clopper_pearson(ok_bin_holes, len(bin_holes), cl)
        e3 = clopper_pearson(ok_bin_broken_holes, len(bin_broken_holes), cl)

        print(f"\n{'='*60}")
        print(f"  {label} (denominador propio)  —  Clopper-Pearson {sigma_label}")
        print(f"{'='*60}")
        print(f"  Total de eventos : {n_events}")
        print(f"{'-'*60}")

        # --- Desglose por tipo: buenos / malos / % sobre Binary+Broken+Holes ---
        # (se omite el tipo "Other")
        n_total_bbh = len(bin_broken_holes)
        print(f"  {'Tipo':<10} {'Buenos':>8} {'Malos':>8} {'Total':>8} {'% (Bin+Brk+Hol)':>16}")
        print(f"  {'-'*54}")
        type_order = ["Binary", "Multi", "Broken", "Holes"]
        present_types = list(df_subset["type"].dropna().unique())
        ordered_types = [t for t in type_order if t in present_types] + [
            t for t in present_types if t not in type_order and t != "Other"
        ]
        for t in ordered_types:
            grp_t = df_subset[df_subset["type"] == t]
            n_t = len(grp_t)
            good_t = int((grp_t["status"] == True).sum())
            bad_t = n_t - good_t
            pct_t = (n_t / n_total_bbh * 100) if n_total_bbh > 0 else 0.0
            print(f"  {t:<10} {good_t:>8} {bad_t:>8} {n_t:>8} {pct_t:>15.2f}%")
        print(f"  {'-'*54}")

        for name, ok, n, (eff, lo, hi, dm, dp) in [
            ("Binary                / Binary", ok_binary, len(binary), e1),
            ("(Binary+Holes)        / (Binary+Holes)", ok_bin_holes, len(bin_holes), e2),
            ("(Binary+Broken+Holes) / (Binary+Broken+Holes)", ok_bin_broken_holes, len(bin_broken_holes), e3),
        ]:
            if n == 0:
                print(f"  {name} : sin eventos")
                continue
            print(
                f"  {name} : k={ok:>5}  n={n:>5}  "
                f"ε = {eff*100:6.2f} %  +{dp*100:5.2f}  -{dm*100:5.2f}  ({sigma_label})"
            )
        print(f"{'='*60}\n")

    # Elige la funcion de reporte segun el modo pedido por linea de comandos
    report_fn = report_own_denominator if USE_OWN_DENOMINATOR else report

    # --- Global ---
    report_fn(df, "EFICIENCIA DE RECONSTRUCCION (global)")

    # --- Por run (opcional) ---
    if args.by_run:
        for run_id, grp in df.groupby("run"):
            report_fn(grp, f"RUN {run_id}")


if __name__ == "__main__":
    main()