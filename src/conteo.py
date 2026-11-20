def crear_contador(n):
    pasos = n - 1
    contador = {
        "n": n,
        "pasos": pasos,
        "div_factores": [0] * pasos,
        "mult_matriz": [0] * pasos,
        "restas_matriz": [0] * pasos,
        "mult_vector": [0] * pasos,
        "restas_vector": [0] * pasos,
        "comp_pivote": [0] * pasos,
        "div_escalado": [0] * pasos,
        "intercambios": [0] * pasos,
        "mult_sust": 0,
        "restas_sust": 0,
        "div_sust": 0,
        "comp_escalas": 0
    }
    return contador


def resumir_contador(contador):
    div_factores = sum(contador["div_factores"])
    mult_matriz = sum(contador["mult_matriz"])
    restas_matriz = sum(contador["restas_matriz"])
    mult_vector = sum(contador["mult_vector"])
    restas_vector = sum(contador["restas_vector"])

    mult_div_elim = div_factores + mult_matriz + mult_vector
    suma_resta_elim = restas_matriz + restas_vector
    flops_elim = mult_div_elim + suma_resta_elim

    mult_div_sust = contador["mult_sust"] + contador["div_sust"]
    suma_resta_sust = contador["restas_sust"]
    flops_sust = mult_div_sust + suma_resta_sust

    return {
        "div_factores": div_factores,
        "mult_matriz": mult_matriz,
        "restas_matriz": restas_matriz,
        "mult_vector": mult_vector,
        "restas_vector": restas_vector,
        "mult_div_elim": mult_div_elim,
        "suma_resta_elim": suma_resta_elim,
        "flops_elim": flops_elim,
        "mult_div_sust": mult_div_sust,
        "suma_resta_sust": suma_resta_sust,
        "flops_sust": flops_sust,
        "total_mult_div": mult_div_elim + mult_div_sust,
        "total_suma_resta": suma_resta_elim + suma_resta_sust,
        "total_flops": flops_elim + flops_sust,
        "comp_pivote": sum(contador["comp_pivote"]),
        "div_escalado": sum(contador["div_escalado"]),
        "intercambios": sum(contador["intercambios"]),
        "comp_escalas": contador["comp_escalas"]
    }


def guardar_conteo(conteo_sin, conteo_escalado, conteo_parcial, ruta):
    total = resumir_contador(conteo_escalado)
    n = conteo_escalado["n"]
    referencia = (4 * n**3 + 15 * n**2 - 13 * n) // 6
    diferencia = 100 * (total["total_flops"] - referencia) / referencia

    with open(ruta, "w", encoding="utf-8") as archivo:
        archivo.write("# Conteo de operaciones\n\n")
        archivo.write("Cada suma, resta, multiplicación o división cuenta como 1 flop. ")
        archivo.write("En cada paso se calcula un multiplicador por fila y se usa para actualizar esa fila y su valor de b. ")
        archivo.write("Las columnas anteriores no se tocan.\n\n")
        archivo.write("| Paso (pivoteo escalado) | Multiplicaciones y divisiones | Sumas y restas | Total |\n")
        archivo.write("| --- | ---: | ---: | ---: |\n")
        for paso in range(conteo_escalado["pasos"]):
            productos = conteo_escalado["div_factores"][paso] + conteo_escalado["mult_matriz"][paso] + conteo_escalado["mult_vector"][paso]
            sumas = conteo_escalado["restas_matriz"][paso] + conteo_escalado["restas_vector"][paso]
            archivo.write(f"| {paso + 1} | {productos} | {sumas} | {productos + sumas} |\n")
        archivo.write(f"| Eliminación | {total['mult_div_elim']} | {total['suma_resta_elim']} | {total['flops_elim']} |\n")
        archivo.write(f"| Sustitución regresiva | {total['mult_div_sust']} | {total['suma_resta_sust']} | {total['flops_sust']} |\n")
        archivo.write(f"| Total | {total['total_mult_div']} | {total['total_suma_resta']} | {total['total_flops']} |\n\n")

        archivo.write(f"El contador dio {total['total_flops']} flops para la eliminación y la sustitución. ")
        archivo.write(f"La referencia del enunciado es {referencia} (diferencia: {diferencia:.2f} %). ")
        archivo.write(f"El término principal, 2n³/3, por sí solo da {2 * n**3 / 3:g}.\n\n")
        archivo.write("Si se escribiera el cero directamente en la columna del pivote, se ahorrarían ")
        archivo.write(f"{2 * total['div_factores']} operaciones. Ese ahorro no se descontó del total medido.\n\n")
        archivo.write("Costo de elegir los pivotes (se cuenta aparte):\n\n")
        
        for nombre, contador in (("Sin pivoteo", conteo_sin), ("Pivoteo parcial", conteo_parcial), ("Pivoteo escalado", conteo_escalado)):
            resumen = resumir_contador(contador)
            archivo.write(f"- {nombre}: {resumen['total_flops']} flops de eliminación y sustitución; ")
            archivo.write(f"{resumen['comp_pivote']} comparaciones para elegir pivote y {resumen['intercambios']} intercambios de filas.")
            if contador is conteo_escalado:
                archivo.write(f" Además, {resumen['comp_escalas']} comparaciones para calcular las escalas ")
                archivo.write(f"y {resumen['div_escalado']} divisiones para comparar los candidatos.")
            archivo.write("\n")
