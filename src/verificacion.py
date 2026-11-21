import math
from decimal import Decimal, getcontext
import numpy as np
from gauss import resolver_lu

getcontext().prec = 50


def dar_formato(valor, digitos):
    decimales = digitos - 1
    if decimales < 0:
        decimales = 0
    return f"{valor:.{decimales}e}"


def medir_residuo(matriz_texto, vector_texto, solucion):
    n = len(vector_texto)
    residuo = []
    norma_infinito = Decimal("0")
    suma_cuadrados = Decimal("0")

    for indice_fila in range(n):
        suma_productos = Decimal("0")
        for indice_columna in range(n):
            suma_productos = suma_productos + (Decimal(matriz_texto[indice_fila][indice_columna]) * Decimal(solucion[indice_columna]))
        residuo_fila = Decimal(vector_texto[indice_fila]) - suma_productos
        residuo.append(residuo_fila)
        if abs(residuo_fila) > norma_infinito:
            norma_infinito = abs(residuo_fila)
        suma_cuadrados = suma_cuadrados + (residuo_fila * residuo_fila)

    return residuo, norma_infinito, suma_cuadrados.sqrt()


def refinar(matriz_texto, vector_texto, resultado_gauss):
    x_inicial = resultado_gauss["solucion"]
    n = len(x_inicial)
    residuo_inicial, norma_inf_antes, norma_2_antes = medir_residuo(matriz_texto, vector_texto, x_inicial)

    residuo_float = []
    for indice_fila in range(n):
        residuo_float.append(float(residuo_inicial[indice_fila]))

    correccion = resolver_lu(resultado_gauss["matriz_inferior"], resultado_gauss["matriz_superior"], resultado_gauss["orden_filas"], residuo_float)

    x_refinada = []
    for indice_fila in range(n):
        x_refinada.append(x_inicial[indice_fila] + correccion[indice_fila])

    residuo_refinado, norma_inf_despues, norma_2_despues = medir_residuo(matriz_texto, vector_texto, x_refinada)

    return {
        "x_inicial": x_inicial,
        "residuo_inicial": residuo_inicial,
        "norma_inf_antes": norma_inf_antes,
        "norma_2_antes": norma_2_antes,
        "correccion": correccion,
        "x_refinada": x_refinada,
        "residuo_refinado": residuo_refinado,
        "norma_inf_despues": norma_inf_despues,
        "norma_2_despues": norma_2_despues
    }


def resolver_base(matriz_texto, vector_texto):
    n = len(vector_texto)
    matriz = []
    vector_b = []
    for indice_fila in range(n):
        fila = []
        for indice_columna in range(n):
            fila.append(Decimal(matriz_texto[indice_fila][indice_columna]))
        matriz.append(fila)
        vector_b.append(Decimal(vector_texto[indice_fila]))

    for columna_pivote in range(n - 1):
        fila_pivote = columna_pivote
        for indice_fila in range(columna_pivote + 1, n):
            if abs(matriz[indice_fila][columna_pivote]) > abs(matriz[fila_pivote][columna_pivote]):
                fila_pivote = indice_fila
        if fila_pivote != columna_pivote:
            matriz[columna_pivote], matriz[fila_pivote] = matriz[fila_pivote], matriz[columna_pivote]
            vector_b[columna_pivote], vector_b[fila_pivote] = vector_b[fila_pivote], vector_b[columna_pivote]

        for indice_fila in range(columna_pivote + 1, n):
            multiplicador = matriz[indice_fila][columna_pivote] / matriz[columna_pivote][columna_pivote]
            for indice_columna in range(columna_pivote, n):
                matriz[indice_fila][indice_columna] = matriz[indice_fila][indice_columna] - (multiplicador * matriz[columna_pivote][indice_columna])
            vector_b[indice_fila] = vector_b[indice_fila] - (multiplicador * vector_b[columna_pivote])

    solucion = [Decimal("0")] * n
    for indice_fila in range(n - 1, -1, -1):
        termino = vector_b[indice_fila]
        for indice_columna in range(indice_fila + 1, n):
            termino = termino - (matriz[indice_fila][indice_columna] * solucion[indice_columna])
        solucion[indice_fila] = termino / matriz[indice_fila][indice_fila]
    return solucion


def medir_errores(solucion, x_referencia):
    n = len(solucion)
    epsilon_maquina = 2.0 ** (-52)
    errores_por_x = []
    cifras_por_x = []
    suma_error_2 = 0.0
    suma_ref_2 = 0.0

    for indice_fila in range(n):
        diferencia_decimal = abs(Decimal(solucion[indice_fila]) - x_referencia[indice_fila])
        diferencia = float(diferencia_decimal)
        valor_referencia = float(x_referencia[indice_fila])
        suma_error_2 = suma_error_2 + (diferencia * diferencia)
        suma_ref_2 = suma_ref_2 + (valor_referencia * valor_referencia)

        error_relativo = float(diferencia_decimal / abs(x_referencia[indice_fila]))
        errores_por_x.append(error_relativo)

        if error_relativo > 0.0:
            cifras_correctas = int(math.floor(-math.log10(2.0 * error_relativo)))
            if cifras_correctas < 0:
                cifras_correctas = 0
            if cifras_correctas > 16:
                cifras_correctas = 16
        else:
            cifras_correctas = 16
        cifras_por_x.append(cifras_correctas)

    error_rel_2 = math.sqrt(suma_error_2) / math.sqrt(suma_ref_2)
    if error_rel_2 > epsilon_maquina:
        bits_perdidos = math.log2(error_rel_2 / epsilon_maquina)
    else:
        bits_perdidos = 0.0

    return {
        "error_rel_2": error_rel_2,
        "bits_perdidos": bits_perdidos,
        "errores_por_x": errores_por_x,
        "cifras_por_x": cifras_por_x,
        "minimo_cifras": min(cifras_por_x),
        "maximo_cifras": max(cifras_por_x)
    }


def medir_condicion(matriz_numerica, vector_numerico):
    matriz_numpy = np.array(matriz_numerica, dtype=np.float64)
    vector_numpy = np.array(vector_numerico, dtype=np.float64)

    condicion_2 = float(np.linalg.cond(matriz_numpy, 2))
    matriz_escalada = []
    for indice_fila in range(len(matriz_numerica)):
        escala_fila = abs(matriz_numerica[indice_fila][0])
        for indice_columna in range(1, len(matriz_numerica[indice_fila])):
            if abs(matriz_numerica[indice_fila][indice_columna]) > escala_fila:
                escala_fila = abs(matriz_numerica[indice_fila][indice_columna])
        fila_escalada = []
        for indice_columna in range(len(matriz_numerica[indice_fila])):
            fila_escalada.append(matriz_numerica[indice_fila][indice_columna] / escala_fila)
        matriz_escalada.append(fila_escalada)
    condicion_esc = float(np.linalg.cond(np.array(matriz_escalada, dtype=np.float64), 2))
    epsilon_maquina = 2.0 ** (-52)
    fila_tres = matriz_numerica[2]
    fila_cinco = matriz_numerica[4]
    fila_siete = matriz_numerica[6]

    determinante_base = (fila_tres[1] * fila_cinco[2]) - (fila_tres[2] * fila_cinco[1])
    alfa = ((fila_siete[1] * fila_cinco[2]) - (fila_siete[2] * fila_cinco[1])) / determinante_base
    beta = ((fila_tres[1] * fila_siete[2]) - (fila_tres[2] * fila_siete[1])) / determinante_base

    mayor_desviacion = 0.0
    for indice_columna in range(len(fila_siete)):
        diferencia_columna = abs(fila_siete[indice_columna] - ((alfa * fila_tres[indice_columna]) + (beta * fila_cinco[indice_columna])))
        if diferencia_columna > mayor_desviacion:
            mayor_desviacion = diferencia_columna

    return {
        "determinante": float(np.linalg.det(matriz_numpy)),
        "condicion_2": condicion_2,
        "condicion_esc": condicion_esc,
        "epsilon_maquina": epsilon_maquina,
        "cota": condicion_2 * epsilon_maquina,
        "solucion_numpy": np.linalg.solve(matriz_numpy, vector_numpy).tolist(),
        "alfa": alfa,
        "beta": beta,
        "mayor_desviacion": mayor_desviacion
    }


def guardar_solucion(x_escalada, errores_escalado, x_refinada, errores_refinados, ruta):
    archivo = open(ruta, "w", encoding="utf-8")
    archivo.write("incognita,valor_defendido_refinado,digitos_afirmados_refinado,error_relativo_refinado,valor_eliminacion_directa,digitos_afirmados_directa,error_relativo_directa\n")
    for indice_fila in range(len(x_escalada)):
        cifras_refinadas = errores_refinados["cifras_por_x"][indice_fila]
        cifras_escaladas = errores_escalado["cifras_por_x"][indice_fila]
        archivo.write(
            f"x{indice_fila + 1},{dar_formato(x_refinada[indice_fila], cifras_refinadas)},{cifras_refinadas},{errores_refinados['errores_por_x'][indice_fila]:.4e},"
            f"{dar_formato(x_escalada[indice_fila], cifras_escaladas)},{cifras_escaladas},{errores_escalado['errores_por_x'][indice_fila]:.4e}\n"
        )
    archivo.close()


def cambio(antes, despues):
    if despues < antes:
        return "bajó"
    if despues > antes:
        return "subió"
    return "no cambió"


def guardar_residuo(refinamiento, residuo_sin, norma_inf_sin, norma_2_sin, norma_inf_numpy, norma_2_numpy, errores_escalado, errores_sin, errores_refinados, errores_numpy, condicionamiento, ruta):
    with open(ruta, "w", encoding="utf-8") as archivo:
        archivo.write("Verificación de la solución\n\n")
        archivo.write(f"Se calculó el residuo r = b - Ax con Decimal y {getcontext().prec} cifras de precisión.\n\n")
        for nombre, norma_infinito, norma_2, errores in (
            ("Con pivoteo escalado, antes de refinar", refinamiento["norma_inf_antes"], refinamiento["norma_2_antes"], errores_escalado),
            ("Después de refinar una vez con los mismos factores L y U", refinamiento["norma_inf_despues"], refinamiento["norma_2_despues"], errores_refinados),
            ("Sin pivoteo", norma_inf_sin, norma_2_sin, errores_sin),
            ("Control con numpy.linalg.solve", norma_inf_numpy, norma_2_numpy, errores_numpy),
        ):
            archivo.write(f"{nombre}:\n")
            archivo.write(f"  Mayor residuo en valor absoluto (norma infinito): {float(norma_infinito):.12e}\n")
            archivo.write(f"  Norma 2 del residuo: {float(norma_2):.12e}\n")
            archivo.write(f"  Error relativo de la solución: {errores['error_rel_2']:.6e}\n")
            archivo.write(f"  Bits perdidos estimados: {errores['bits_perdidos']:.2f}\n")
            archivo.write(f"  Cifras correctas: entre {errores['minimo_cifras']} y {errores['maximo_cifras']}\n\n")

        archivo.write(f"Al refinar, el error relativo {cambio(errores_escalado['error_rel_2'], errores_refinados['error_rel_2'])} ")
        archivo.write(f"y la norma infinito del residuo {cambio(refinamiento['norma_inf_antes'], refinamiento['norma_inf_despues'])}.\n\n")
        archivo.write("Residuo por fila:\n")
        archivo.write("Fila   Antes de refinar            Después de refinar          Sin pivoteo\n")
        for fila, (antes, despues, valor_sin) in enumerate(zip(refinamiento["residuo_inicial"], refinamiento["residuo_refinado"], residuo_sin), start=1):
            archivo.write(f"{fila:>4}   {float(antes):>24.18e}   {float(despues):>24.18e}   {float(valor_sin):>20.14e}\n")

        archivo.write("\nCorrección de cada incógnita:\n")
        for indice, valor in enumerate(refinamiento["correccion"]):
            cifras_antes = errores_escalado["cifras_por_x"][indice]
            cifras_despues = errores_refinados["cifras_por_x"][indice]
            antes = dar_formato(refinamiento["x_inicial"][indice], cifras_antes)
            despues = dar_formato(refinamiento["x_refinada"][indice], cifras_despues)
            archivo.write(f"x{indice + 1}: se sumó {valor:.6e}\n")
            archivo.write(f"  Antes:   {antes} ({cifras_antes} cifras correctas)\n")
            archivo.write(f"  Después: {despues} ({cifras_despues} cifras correctas)\n")

        archivo.write("\nCondicionamiento de la matriz:\n")
        archivo.write(f"  Determinante de A: {condicionamiento['determinante']:.6e}\n")
        archivo.write(f"  Número de condición (norma 2): {condicionamiento['condicion_2']:.6e}\n")
        archivo.write(f"  Número de condición al escalar las filas: {condicionamiento['condicion_esc']:.6e}\n")
        archivo.write(f"  Épsilon de máquina: {condicionamiento['epsilon_maquina']:.6e}\n")
        archivo.write(f"  Cota para el error relativo (condición por épsilon): {condicionamiento['cota']:.6e}\n")
        archivo.write(f"La fila 7 se aproxima a {condicionamiento['alfa']:.5f} veces la fila 3 más {condicionamiento['beta']:.5f} veces la fila 5. ")
        archivo.write(f"La mayor diferencia en valor absoluto es {condicionamiento['mayor_desviacion']:.6e}.\n")
