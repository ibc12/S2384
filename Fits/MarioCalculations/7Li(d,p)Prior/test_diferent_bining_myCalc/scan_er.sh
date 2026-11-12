#!/usr/bin/env bash
# ---------------------------------------------------------------
# scan_er.sh
# Crea una carpeta por cada valor de 'er' (primer &Overlap del input),
# genera el input a partir de la plantilla y ejecuta FRESCO.
#
# Uso:   ./scan_er.sh            -> crea carpetas, inputs y ejecuta fresco
#        RUN=0 ./scan_er.sh      -> solo crea carpetas e inputs (sin ejecutar)
#        NJOBS=4 ./scan_er.sh    -> ejecuta 4 cálculos en paralelo
#        LINK=1 ./scan_er.sh     -> enlaza fort.4 en vez de copiarlo
# ---------------------------------------------------------------

# ---- CONFIGURACIÓN ----------------------------------------------
BASE="7Lidp8Li_prior"          # nombre base (plantilla: ${BASE}.in)
TEMPLATE="${BASE}.in"
# Se dan en valor absoluto (MeV): el script escribe siempre er NEGATIVO
EMIN=0.1                       # |er| mínimo (MeV)  -> er = -0.1
EMAX=3.0                       # |er| máximo (MeV)  -> er = -3.0
ESTEP=0.1                      # paso en |er| (MeV)
FRESCO_CMD="fresco"            # comando de FRESCO
EXTRA_FILES="fort.4"           # ficheros a copiar en cada carpeta (potencial externo)
                               # (varios: EXTRA_FILES="fort.4 otro.dat")
LINK="${LINK:-0}"              # 1 = enlace simbólico en vez de copia (ahorra disco)
RUN="${RUN:-1}"                # 1 = ejecutar fresco, 0 = solo preparar
NJOBS="${NJOBS:-1}"            # nº de cálculos simultáneos

# Si prefieres una lista explícita de valores (valor absoluto), descomenta y edita:
# ER_LIST="0.1 0.2 0.3 0.5 0.75 1.0 1.5 2.0 3.0"
# -----------------------------------------------------------------

export LC_NUMERIC=C            # punto decimal siempre

if [[ ! -f "$TEMPLATE" ]]; then
    echo "ERROR: no encuentro la plantilla '$TEMPLATE'" >&2
    exit 1
fi

for F in $EXTRA_FILES; do
    if [[ ! -f "$F" ]]; then
        echo "ERROR: no encuentro '$F' en $(pwd)" >&2
        exit 1
    fi
done

if [[ -z "$ER_LIST" ]]; then
    ER_LIST=$(seq -f "%.2f" "$EMIN" "$ESTEP" "$EMAX")
fi

for ER in $ER_LIST; do
    # fuerza er negativo (acepta valores dados en positivo o en negativo)
    ERVAL=$(awk -v e="$ER" 'BEGIN{ if (e > 0) e = -e; printf "%.5f", e }')
    ERTAG=$(printf "%.2f" "$ERVAL")          # p.ej. -0.10
    NAME="${BASE}_er${ERTAG}"
    DIR="er_${ERTAG}MeV"

    mkdir -p "$DIR"

    # Sustituye SOLO la primera aparición de 'er=' (primer &Overlap)
    ER_NEW="$ERVAL" perl -0pe \
        's/\b(er=\s*)[-+]?[0-9]*\.?[0-9]+(?:[eEdD][-+]?[0-9]+)?/$1.$ENV{ER_NEW}/ie' \
        "$TEMPLATE" > "$DIR/${NAME}.in"

    # copia (o enlaza) los ficheros auxiliares, p.ej. fort.4
    for F in $EXTRA_FILES; do
        if [[ "$LINK" == "1" ]]; then
            ln -sf "$(pwd)/$F" "$DIR/$F"
        else
            cp -f "$F" "$DIR/$F"
        fi
    done

    echo "[$DIR] $(grep -io 'er= *[-+0-9.]*' "$DIR/${NAME}.in" | head -1)"

    if [[ "$RUN" == "1" ]]; then
        # limitar nº de trabajos simultáneos
        while (( $(jobs -rp | wc -l) >= NJOBS )); do sleep 0.5; done
        ( cd "$DIR" && "$FRESCO_CMD" < "${NAME}.in" > "${NAME}.out" 2> "${NAME}.err"
          echo "   terminado: $DIR" ) &
        # Si NJOBS=1 espera a que acabe antes de lanzar el siguiente
        (( NJOBS == 1 )) && wait
    fi
done

wait
echo "Hecho."