"""
efficiency_clopper_pearson.py
-------------------------------
Lee un CSV con columnas (run, event, Type, status),
filtra eventos de tipo "Binary" y calcula la eficiencia
de reconstruccion con intervalos de Clopper-Pearson a 2 sigma (95.45% CL).

Uso:
    python efficiency_clopper_pearson.py mi_archivo.csv
    python efficiency_clopper_pearson.py mi_archivo.csv --cl 0.6827   # 1 sigma
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
    n  : numero total de eventos
    cl : nivel de confianza (default 0.9545 = 2 sigma)

    Retorna
    -------
    eff        : eficiencia central k/n
    lo, hi     : limites inferior y superior del intervalo
    delta_lo   : incertidumbre inferior  (eff - lo)
    delta_hi   : incertidumbre superior  (hi  - eff)
    """
    if n == 0:
        raise ValueError("El numero total de eventos es 0.")
    if k < 0 or k > n:
        raise ValueError(f"k={k} fuera del rango [0, n={n}].")

    alpha = 1.0 - cl

    # Casos extremos
    lo = 0.0 if k == 0 else beta.ppf(alpha / 2.0, k, n - k + 1)
    hi = 1.0 if k == n else beta.ppf(1.0 - alpha / 2.0, k + 1, n - k)

    eff = k / n
    return eff, lo, hi, eff - lo, hi - eff


# ---------------------------------------------------------------------------
# Lectura y filtrado del CSV
# ---------------------------------------------------------------------------
def load_binary_events(filepath: str) -> pd.DataFrame:
    """
    Lee el CSV y devuelve solo los eventos con Type == 'Binary'.
    La columna Type puede estar vacia; esas filas se descartan.
    La columna status se convierte a booleano:
        True / 'True' / 1 -> bien reconstruido
        False / 'False' / 0 -> mal reconstruido
    """
    df = pd.read_csv(
        filepath,
        dtype={"run": str, "entry": str, "type": str, "status": str},
        skipinitialspace=True,
    )

    # Normalizar nombres de columna (strip espacios, minusculas)
    df.columns = df.columns.str.strip().str.lower()

    # Verificar columnas requeridas
    required = {"run", "entry", "type", "status"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Columnas faltantes en el CSV: {missing}")

    # Filtrar solo Binary (ignorar NaN y cadenas vacias)
    df["type"] = df["type"].str.strip()
    df_binary = df[df["type"].str.lower() == "binary"].copy()

    if df_binary.empty:
        raise ValueError("No se encontraron eventos con Type == 'Binary' en el archivo.")

    # Convertir status a booleano
    df_binary["status"] = df_binary["status"].str.strip().str.lower().map(
        {"true": True, "false": False, "1": True, "0": False}
    )

    n_invalid = df_binary["status"].isna().sum()
    if n_invalid > 0:
        print(
            f"[Aviso] {n_invalid} filas con valor de status no reconocido "
            f"seran ignoradas."
        )
        df_binary = df_binary.dropna(subset=["status"])

    return df_binary


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Calcula eficiencia con intervalos Clopper-Pearson para eventos Binary."
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
        help="Calcula la eficiencia tambien por run",
    )
    args = parser.parse_args()

    # --- Cargar datos ---
    print(f"\nLeyendo: {args.csv}")
    try:
        df = load_binary_events(args.csv)
    except Exception as e:
        print(f"[Error] {e}")
        sys.exit(1)

    n_total_csv  = pd.read_csv(args.csv, skipinitialspace=True).shape[0]
    n_binary     = len(df)
    n_good       = int(df["status"].sum())
    n_bad        = n_binary - n_good
    cl           = args.cl
    sigma_label  = {0.6827: "1σ", 0.9545: "2σ", 0.9973: "3σ"}.get(cl, f"{cl*100:.2f}%")

    print(f"\n{'='*55}")
    print(f"  EFICIENCIA DE RECONSTRUCCION  —  Clopper-Pearson {sigma_label}")
    print(f"{'='*55}")
    print(f"  Filas totales en CSV          : {n_total_csv}")
    print(f"  Eventos Binary                : {n_binary}")
    print(f"  Bien reconstruidos (True)     : {n_good}")
    print(f"  Mal reconstruidos  (False)    : {n_bad}")
    print(f"  Nivel de confianza            : {cl} ({sigma_label})")
    print(f"{'-'*55}")

    # --- Eficiencia global ---
    eff, lo, hi, dm, dp = clopper_pearson(n_good, n_binary, cl)

    print(f"\n  Eficiencia central  : {eff*100:.4f} %")
    print(f"  Intervalo inferior  : {lo*100:.4f} %")
    print(f"  Intervalo superior  : {hi*100:.4f} %")
    print(f"\n  ε = {eff*100:.2f}  +{dp*100:.2f}  -{dm*100:.2f}  %  ({sigma_label})")
    print(f"\n  (Comparacion gaussiana: ±{100*(eff*(1-eff)/n_binary)**0.5:.4f} %)")
    print(f"{'='*55}\n")

    # --- Eficiencia por run (opcional) ---
    if args.by_run:
        print(f"  EFICIENCIA POR RUN  ({sigma_label})\n")
        print(f"  {'Run':<12} {'k':>6} {'n':>6}  {'ε (%)':>8}  {'−Δ':>8}  {'+Δ':>8}")
        print(f"  {'-'*58}")
        for run_id, grp in df.groupby("run"):
            k_r = int(grp["status"].sum())
            n_r = len(grp)
            e_r, lo_r, hi_r, dm_r, dp_r = clopper_pearson(k_r, n_r, cl)
            print(
                f"  {str(run_id):<12} {k_r:>6} {n_r:>6}  "
                f"{e_r*100:>8.3f}  -{dm_r*100:>7.3f}  +{dp_r*100:>7.3f}"
            )
        print()


if __name__ == "__main__":
    main()
