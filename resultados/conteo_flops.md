# Conteo de operaciones

Cada suma, resta, multiplicación o división cuenta como 1 flop. En cada paso se calcula un multiplicador por fila y se usa para actualizar esa fila y su valor de b. Las columnas anteriores no se tocan.

| Paso (pivoteo escalado) | Multiplicaciones y divisiones | Sumas y restas | Total |
| --- | ---: | ---: | ---: |
| 1 | 154 | 143 | 297 |
| 2 | 130 | 120 | 250 |
| 3 | 108 | 99 | 207 |
| 4 | 88 | 80 | 168 |
| 5 | 70 | 63 | 133 |
| 6 | 54 | 48 | 102 |
| 7 | 40 | 35 | 75 |
| 8 | 28 | 24 | 52 |
| 9 | 18 | 15 | 33 |
| 10 | 10 | 8 | 18 |
| 11 | 4 | 3 | 7 |
| Eliminación | 704 | 638 | 1342 |
| Sustitución regresiva | 78 | 66 | 144 |
| Total | 782 | 704 | 1486 |

El contador dio 1486 flops para la eliminación y la sustitución. La referencia del enunciado es 1486 (diferencia: 0.00 %). El término principal, 2n³/3, por sí solo da 1152.

Si se escribiera el cero directamente en la columna del pivote, se ahorrarían 132 operaciones. Ese ahorro no se descontó del total medido.

Costo de elegir los pivotes (se cuenta aparte):

- Sin pivoteo: 1486 flops de eliminación y sustitución; 0 comparaciones para elegir pivote y 0 intercambios de filas.
- Pivoteo parcial: 1486 flops de eliminación y sustitución; 66 comparaciones para elegir pivote y 9 intercambios de filas.
- Pivoteo escalado: 1486 flops de eliminación y sustitución; 66 comparaciones para elegir pivote y 8 intercambios de filas. Además, 132 comparaciones para calcular las escalas y 77 divisiones para comparar los candidatos.
