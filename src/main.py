from pathlib import Path
from conteo import crear_contador, resumir_contador, guardar_conteo
from gauss import resolver_gauss, guardar_traza
from verificacion import (
    medir_residuo,
    refinar,
    resolver_base,
    medir_errores,
    medir_condicion,
    guardar_solucion,
    guardar_residuo
)


def cargar_datos(ruta_matriz, ruta_vector):
    matriz_texto = []
    matriz_numerica = []
    with open(ruta_matriz, encoding="utf-8") as archivo_matriz:
        for linea in archivo_matriz:
            if linea.strip():
                fila_texto = [campo.strip() for campo in linea.split(",")]
                matriz_texto.append(fila_texto)
                matriz_numerica.append([float(valor) for valor in fila_texto])

    vector_texto = []
    vector_numerico = []
    with open(ruta_vector, encoding="utf-8") as archivo_vector:
        for linea in archivo_vector:
            valor = linea.strip()
            if valor:
                vector_texto.append(valor)
                vector_numerico.append(float(valor))

    return matriz_texto, matriz_numerica, vector_texto, vector_numerico


def main():
    matriz_texto, matriz_numerica, vector_texto, vector_numerico = cargar_datos("datos/matriz_a.csv", "datos/vector_b.csv")
    n = len(vector_numerico)

    conteo_sin = crear_contador(n)
    sin_pivoteo = resolver_gauss(matriz_numerica, vector_numerico, "ninguno", conteo_sin)

    conteo_parcial = crear_contador(n)
    parcial = resolver_gauss(matriz_numerica, vector_numerico, "parcial", conteo_parcial)

    conteo_escalado = crear_contador(n)
    escalado = resolver_gauss(matriz_numerica, vector_numerico, "escalado", conteo_escalado)

    x_referencia = resolver_base(matriz_texto, vector_texto)
    refinamiento = refinar(matriz_texto, vector_texto, escalado)
    residuo_sin, norma_inf_sin, norma_2_sin = medir_residuo(matriz_texto, vector_texto, sin_pivoteo["solucion"])
    condicionamiento = medir_condicion(matriz_numerica, vector_numerico)
    residuo_numpy, norma_inf_numpy, norma_2_numpy = medir_residuo(matriz_texto, vector_texto, condicionamiento["solucion_numpy"])

    errores_sin = medir_errores(sin_pivoteo["solucion"], x_referencia)
    errores_parcial = medir_errores(parcial["solucion"], x_referencia)
    errores_escalado = medir_errores(escalado["solucion"], x_referencia)
    errores_refinados = medir_errores(refinamiento["x_refinada"], x_referencia)
    errores_numpy = medir_errores(condicionamiento["solucion_numpy"], x_referencia)

    Path("resultados").mkdir(exist_ok=True)
    guardar_solucion(escalado["solucion"], errores_escalado, refinamiento["x_refinada"], errores_refinados, "resultados/solucion.csv")
    guardar_traza(escalado, sin_pivoteo, "resultados/traza_eliminacion.txt")
    guardar_conteo(conteo_sin, conteo_escalado, conteo_parcial, "resultados/conteo_flops.md")
    guardar_residuo(
        refinamiento, residuo_sin, norma_inf_sin, norma_2_sin, norma_inf_numpy, norma_2_numpy,
        errores_escalado, errores_sin, errores_refinados, errores_numpy, condicionamiento, "resultados/residuo.txt"
    )

    resumen = resumir_contador(conteo_escalado)
    print(f"La eliminación y la sustitución usaron {resumen['total_flops']} operaciones:")
    print(f"{resumen['total_mult_div']} multiplicaciones o divisiones y {resumen['total_suma_resta']} sumas o restas.")
    print(f"Para elegir los pivotes se hicieron {resumen['comp_pivote']} comparaciones.")

    for nombre, resultado, errores in (
        ("Sin pivoteo", sin_pivoteo, errores_sin),
        ("Con pivoteo parcial", parcial, errores_parcial),
        ("Con pivoteo escalado", escalado, errores_escalado),
    ):
        print(f"\n{nombre}:")
        print(f"  El mayor multiplicador en valor absoluto fue {resultado['max_mult']:.5g}.")
        if nombre == "Con pivoteo escalado":
            print(f"  Al tener en cuenta la escala de cada fila, el mayor fue {resultado['max_mult_esc']:.5g}.")
        print(f"  El factor de crecimiento de la matriz fue {resultado['crecimiento']:.5g}.")
        print(f"  El error relativo fue {errores['error_rel_2']:.4e}.")
        print(f"  Se obtuvieron entre {errores['minimo_cifras']} y {errores['maximo_cifras']} cifras correctas.")

    print(f"\nNúmero de condición de la matriz (norma 2): {condicionamiento['condicion_2']:.6e}")
    print("Antes de refinar:")
    print(f"  Mayor residuo en valor absoluto: {float(refinamiento['norma_inf_antes']):.12e}")
    print(f"  Error relativo: {errores_escalado['error_rel_2']:.6e}")
    print("Después de una pasada de refinamiento:")
    print(f"  Mayor residuo en valor absoluto: {float(refinamiento['norma_inf_despues']):.12e}")
    print(f"  Error relativo: {errores_refinados['error_rel_2']:.6e}")


if __name__ == "__main__":
    main()
