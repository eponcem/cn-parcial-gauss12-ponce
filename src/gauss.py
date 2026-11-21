def copiar_matriz(matriz):
    copia = []
    for indice_fila in range(len(matriz)):
        fila = []
        for indice_columna in range(len(matriz[indice_fila])):
            fila.append(float(matriz[indice_fila][indice_columna]))
        copia.append(fila)
    return copia


def maximo_abs(matriz):
    mayor_valor = 0.0
    for indice_fila in range(len(matriz)):
        for indice_columna in range(len(matriz[indice_fila])):
            valor = abs(matriz[indice_fila][indice_columna])
            if valor > mayor_valor:
                mayor_valor = valor
    return mayor_valor


def resolver_gauss(matriz_original, vector_original, modo, contador):
    n = len(matriz_original)
    matriz = copiar_matriz(matriz_original)
    vector_b = []
    orden_filas = []
    matriz_inferior = []

    for indice_fila in range(n):
        vector_b.append(float(vector_original[indice_fila]))
        orden_filas.append(indice_fila)
        fila_L = [0.0] * n
        fila_L[indice_fila] = 1.0
        matriz_inferior.append(fila_L)

    maximo_inicial = maximo_abs(matriz)
    maximo_alcanzado = maximo_inicial
    max_mult = 0.0
    max_mult_esc = 0.0

    escalas = []
    for indice_fila in range(n):
        escala_fila = abs(matriz[indice_fila][0])
        for indice_columna in range(1, n):
            if modo == "escalado":
                contador["comp_escalas"] = contador["comp_escalas"] + 1
            if abs(matriz[indice_fila][indice_columna]) > escala_fila:
                escala_fila = abs(matriz[indice_fila][indice_columna])
        escalas.append(escala_fila)

    traza = []

    for columna_pivote in range(n - 1):
        fila_pivote = columna_pivote
        cambio = False

        if modo == "parcial":
            mejor = abs(matriz[columna_pivote][columna_pivote])
            for indice_fila in range(columna_pivote + 1, n):
                contador["comp_pivote"][columna_pivote] = contador["comp_pivote"][columna_pivote] + 1
                if abs(matriz[indice_fila][columna_pivote]) > mejor:
                    mejor = abs(matriz[indice_fila][columna_pivote])
                    fila_pivote = indice_fila

        elif modo == "escalado":
            contador["div_escalado"][columna_pivote] = contador["div_escalado"][columna_pivote] + 1
            mejor = abs(matriz[columna_pivote][columna_pivote]) / escalas[columna_pivote]
            for indice_fila in range(columna_pivote + 1, n):
                contador["div_escalado"][columna_pivote] = contador["div_escalado"][columna_pivote] + 1
                contador["comp_pivote"][columna_pivote] = contador["comp_pivote"][columna_pivote] + 1
                razon = abs(matriz[indice_fila][columna_pivote]) / escalas[indice_fila]
                if razon > mejor:
                    mejor = razon
                    fila_pivote = indice_fila

        if fila_pivote != columna_pivote:
            cambio = True
            contador["intercambios"][columna_pivote] = contador["intercambios"][columna_pivote] + 1
            matriz[columna_pivote], matriz[fila_pivote] = matriz[fila_pivote], matriz[columna_pivote]
            vector_b[columna_pivote], vector_b[fila_pivote] = vector_b[fila_pivote], vector_b[columna_pivote]
            escalas[columna_pivote], escalas[fila_pivote] = escalas[fila_pivote], escalas[columna_pivote]
            orden_filas[columna_pivote], orden_filas[fila_pivote] = orden_filas[fila_pivote], orden_filas[columna_pivote]
            for columna in range(0, columna_pivote):
                matriz_inferior[columna_pivote][columna], matriz_inferior[fila_pivote][columna] = matriz_inferior[fila_pivote][columna], matriz_inferior[columna_pivote][columna]

        valor_pivote = matriz[columna_pivote][columna_pivote]
        factores_paso = []

        for indice_fila in range(columna_pivote + 1, n):
            contador["div_factores"][columna_pivote] = contador["div_factores"][columna_pivote] + 1
            multiplicador = matriz[indice_fila][columna_pivote] / valor_pivote
            matriz_inferior[indice_fila][columna_pivote] = multiplicador
            factores_paso.append((indice_fila + 1, multiplicador))

            if abs(multiplicador) > max_mult:
                max_mult = abs(multiplicador)

            factor_escalado = abs(multiplicador) * (escalas[columna_pivote] / escalas[indice_fila])
            if factor_escalado > max_mult_esc:
                max_mult_esc = factor_escalado

            for indice_columna in range(columna_pivote, n):
                contador["mult_matriz"][columna_pivote] = contador["mult_matriz"][columna_pivote] + 1
                producto = multiplicador * matriz[columna_pivote][indice_columna]
                contador["restas_matriz"][columna_pivote] = contador["restas_matriz"][columna_pivote] + 1
                matriz[indice_fila][indice_columna] = matriz[indice_fila][indice_columna] - producto
                if indice_columna == columna_pivote:
                    matriz[indice_fila][columna_pivote] = 0.0

            contador["mult_vector"][columna_pivote] = contador["mult_vector"][columna_pivote] + 1
            producto_vector = multiplicador * vector_b[columna_pivote]
            contador["restas_vector"][columna_pivote] = contador["restas_vector"][columna_pivote] + 1
            vector_b[indice_fila] = vector_b[indice_fila] - producto_vector

        maximo_paso = maximo_abs(matriz)
        if maximo_paso > maximo_alcanzado:
            maximo_alcanzado = maximo_paso

        traza.append({
            "paso": columna_pivote + 1,
            "fila_pivote": fila_pivote + 1,
            "cambio": cambio,
            "valor_pivote": valor_pivote,
            "multiplicadores": factores_paso,
            "maximo_paso": maximo_paso,
            "matriz": copiar_matriz(matriz),
            "vector": list(vector_b)
        })

    matriz_superior = copiar_matriz(matriz)
    solucion = [0.0] * n

    for indice_fila in range(n - 1, -1, -1):
        termino = vector_b[indice_fila]
        for indice_columna in range(indice_fila + 1, n):
            contador["mult_sust"] = contador["mult_sust"] + 1
            producto_solucion = matriz[indice_fila][indice_columna] * solucion[indice_columna]
            contador["restas_sust"] = contador["restas_sust"] + 1
            termino = termino - producto_solucion
        contador["div_sust"] = contador["div_sust"] + 1
        solucion[indice_fila] = termino / matriz[indice_fila][indice_fila]

    no_nulos = 0
    for indice_fila in range(1, n):
        for indice_columna in range(0, indice_fila):
            if abs(matriz_superior[indice_fila][indice_columna]) > 0.0:
                no_nulos = no_nulos + 1

    return {
        "solucion": solucion,
        "matriz_inferior": matriz_inferior,
        "matriz_superior": matriz_superior,
        "orden_filas": orden_filas,
        "max_mult": max_mult,
        "max_mult_esc": max_mult_esc,
        "maximo_inicial": maximo_inicial,
        "maximo_alcanzado": maximo_alcanzado,
        "crecimiento": maximo_alcanzado / maximo_inicial,
        "no_nulos_debajo": no_nulos,
        "traza": traza
    }


def resolver_lu(matriz_inferior, matriz_superior, orden_filas, residuo):
    n = len(residuo)
    y_intermedia = [0.0] * n
    for indice_fila in range(n):
        suma = 0.0
        for indice_columna in range(0, indice_fila):
            suma = suma + matriz_inferior[indice_fila][indice_columna] * y_intermedia[indice_columna]
        y_intermedia[indice_fila] = residuo[orden_filas[indice_fila]] - suma

    correccion = [0.0] * n
    for indice_fila in range(n - 1, -1, -1):
        suma = 0.0
        for indice_columna in range(indice_fila + 1, n):
            suma = suma + matriz_superior[indice_fila][indice_columna] * correccion[indice_columna]
        correccion[indice_fila] = (y_intermedia[indice_fila] - suma) / matriz_superior[indice_fila][indice_fila]
    return correccion


def escribir_pasos(archivo, titulo, resultado):
    archivo.write(f"{titulo}\n\n")
    for paso in resultado["traza"]:
        numero = paso["paso"]
        archivo.write(f"Paso {numero} de {len(resultado['traza'])}\n")
        if paso["cambio"]:
            archivo.write(f"Se intercambiaron las filas {numero} y {paso['fila_pivote']}.\n")
        else:
            archivo.write("No se intercambiaron filas.\n")
        archivo.write(f"Pivote usado: {paso['valor_pivote']:.8e}\n")
        archivo.write(f"Mayor valor absoluto después del paso: {paso['maximo_paso']:.8e}\n")
        archivo.write("Multiplicadores de este paso:\n")
        for fila, valor in paso["multiplicadores"]:
            archivo.write(f"  Fila {fila}: {valor:.4e}\n")

        archivo.write("Matriz y vector b después del paso:\n")
        for fila, (valores, termino) in enumerate(zip(paso["matriz"], paso["vector"]), start=1):
            texto_fila = " ".join(f"{valor:13.4e}" for valor in valores)
            archivo.write(f"F{fila:02d}: [{texto_fila}] | b={termino:13.5e}\n")
        archivo.write("\n")

    archivo.write(f"Al terminar ({titulo.lower()}):\n")
    archivo.write(f"  Mayor valor absoluto de A al inicio: {resultado['maximo_inicial']:.8e}\n")
    archivo.write(f"  Mayor valor absoluto durante la eliminación: {resultado['maximo_alcanzado']:.8e}\n")
    archivo.write(f"  Factor de crecimiento: {resultado['crecimiento']:.6e}\n")
    archivo.write(f"  Mayor multiplicador en valor absoluto: {resultado['max_mult']:.6e}\n")
    archivo.write(f"  Mayor multiplicador al ajustar por la escala: {resultado['max_mult_esc']:.6e}\n")
    n = len(resultado["solucion"])
    archivo.write(f"  Elementos distintos de cero bajo la diagonal: {resultado['no_nulos_debajo']} de {n * (n - 1) // 2}.\n")


def guardar_traza(escalado, sin_pivoteo, ruta):
    with open(ruta, "w", encoding="utf-8") as archivo:
        archivo.write("Pasos de la eliminación gaussiana\n\n")
        escribir_pasos(archivo, "Con pivoteo escalado", escalado)
        archivo.write("\n")
        escribir_pasos(archivo, "Sin pivoteo", sin_pivoteo)
