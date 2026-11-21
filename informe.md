# Eliminación gaussiana con mínimo costo operacional y máxima exactitud de mantisa

## 1. Identificación

- Nombre: Ernesto Maximiliano Ponce Muñoz
- Asignatura: Computación Numérica (INF-321)
- Carrera: Ingeniería Civil Informática
- Facultad de Ciencias de la Ingeniería, Universidad Católica del Maule
- Docentes: Juan Manuel Cabezas Contardo y Pedro Pesante
- Fecha: Talca, noviembre de 2026
- [Repositorio del trabajo](https://github.com/eponcem/cn-parcial-gauss12-ponce)

## 2. Planteamiento

Se resuelve $Ax=b$ con 12 incógnitas. Al ser una matriz completamente densa (144 de 144 entradas no nulas), no hay ceros previos que permitan saltar operaciones: el costo mínimo queda fijado en orden $\frac{2}{3}n^3$ flops y cualquier ahorro debe provenir únicamente de eliminar trabajo redundante. Los valores de $b$ se leen con todas las cifras del enunciado y $\det(A)\approx 1.364551\times10^9\neq 0$ garantiza solución única. En doble precisión, tres propiedades dificultan calcularla con exactitud.

- **Escalas de fila.** Las escalas $s_i=\max_j\lvert a_{ij}\rvert$ dividen la matriz en dos regímenes: las filas impares ($F_1, F_3, F_5, F_7, F_9, F_{11}$) se ubican entre $8.050\times10^{-5}$ ($F_5$) y $6.334\times10^{-2}$ ($F_7$), mientras que las filas pares ($F_2, F_4, F_6, F_8, F_{10}, F_{12}$) van de $8.383\times10^3$ ($F_4$) a $9.669\times10^5$ ($F_6$). Esta brecha de hasta diez órdenes de magnitud hace que comparar pivotes en valor absoluto sin escalar elija mal y degrade la mantisa.
- **Pivote inicial.** Dividir $a_{6,1}=3.69040121\times10^5$ por $a_{1,1}=1.3407\times10^{-8}$ produce un multiplicador de orden $10^{13}$ en el primer paso.
- **Casi dependencia lineal.** A partir de las columnas 2 y 3 se obtiene $F_7\approx 4.173F_3+582.7F_5$. Las diferencias más grandes aparecen en las columnas 1, 6 y 12 y valen aproximadamente $+2.0\times10^{-4}$, $-6.0\times10^{-5}$ y $+3.0\times10^{-5}$. Esta cercanía a una dependencia lineal contribuye al mal condicionamiento.

## 3. Estrategia de eliminación

Los 11 pasos anulan las 66 entradas bajo la diagonal sin exceder el mínimo teórico de la eliminación densa porque se elimina todo trabajo redundante:

- **Reutilización del multiplicador:** para cada fila $i>k$, src/gauss.py calcula $m_{ik}=a_{ik}/a_{kk}$ una sola vez, lo guarda en $L$ y lo reutiliza en toda la fila activa y en $b_i$; recalcularlo por columna añadiría 506 divisiones redundantes en $A$ (572 en vez de 66).
- **Barrido restringido a la submatriz activa:** el ciclo interno parte en la columna $k$ y no vuelve a operar las columnas $1,\dots,k-1$ ya anuladas, evitando 440 flops redundantes.
- **Una sola factorización $PA=LU$ y nunca $A^{-1}$:** nunca se calcula la matriz inversa (lo que elevaría el costo a $\ge 2n^3=3456$ flops, o $\approx 7000+$ sin optimizar, y deterioraría el error). La corrida principal obtiene $PA=LU$ una sola vez y reutiliza esos factores durante el refinamiento iterativo; las corridas sin pivoteo y con pivoteo parcial simple se ejecutan aparte únicamente para comparar el efecto del pivoteo.

En el pivoteo parcial escalado se calcula una escala inicial $s_i=\max_j\lvert a_{ij}\rvert$ por fila y en cada paso se elige la fila con mayor razón $\lvert a_{ik}\rvert/s_i$, permutando también $b$, las escalas, la permutación y los multiplicadores previos de $L$. El orden de filas originales elegidas como pivote a lo largo de los 11 pasos es $F_8, F_2, F_1, F_6, F_9, F_{11}, F_{12}, F_5, F_4, F_{10}, F_3$ (quedando la fila casi dependiente $F_7$ al final, en la posición 12). La traza completa de los 11 pasos se guarda en resultados/traza_eliminacion.txt.

## 4. Conteo de operaciones

Cada suma, resta, multiplicación o división cuenta como un flop. El contador registra las operaciones sobre $A$ y $b$ en cada paso, más la sustitución regresiva. La búsqueda e intercambio de pivotes se contabilizan aparte.

| Paso | Mult. y div. | Sumas y restas | Total |
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

El término dominante $2n^3/3$ vale 1152 para $n=12$. Con la convención del enunciado, el contador da **1486 flops** (1342 de eliminación y 144 de sustitución), exactamente la referencia. Asignar el cero directamente en la columna del pivote ahorraría dos operaciones por cada una de las 66 entradas, lo que reduciría el total a 1354.

El pivoteo parcial simple añade 66 comparaciones y 9 intercambios. El escalado requiere 66 comparaciones de candidatos, 8 intercambios, 132 comparaciones iniciales para obtener las escalas y 77 divisiones para calcular las razones. El detalle generado por el programa se guarda en resultados/conteo_flops.md.

## 5. Análisis de mantisa y pivoteo

Las tres corridas usan el mismo código y parten del mismo máximo inicial, $9.668707\times10^5$. El factor de crecimiento $\rho$ divide el máximo observado durante la eliminación por ese valor. Los bits perdidos se estiman a partir del error relativo $e$ mediante $\max(0,\log_2(e/\varepsilon))$, con $\varepsilon=2^{-52}$.

| Indicador numérico | Sin pivoteo | Parcial simple | Parcial escalado |
| :--- | :---: | :---: | :---: |
| Primer pivote elegido | $a_{1,1} = 1.3407 \times 10^{-8}$ ($F_1$) | $a_{6,1} = 3.6904 \times 10^{5}$ ($F_6$) | $a_{8,1} = -8.295805 \times 10^{4}$ ($F_8$) |
| Multiplicador bruto máximo $\max \lvert m_{ik}\rvert$ | $2.752593 \times 10^{13}$ | $0.619882$ | $2.465558 \times 10^{9}$ |
| Multiplicador escalado máximo $\max \lvert m_{ik}\rvert\frac{s_p}{s_i}$ | $6.917560 \times 10^{3}$ | $8.225432$ | $\mathbf{0.962129}$ |
| Máximo histórico en eliminación $\max_{k,i,j} \lvert a^{(k)}_{ij}\rvert$ | $2.553722 \times 10^{9}$ | $9.668707 \times 10^{5}$ | $1.724126 \times 10^{6}$ |
| Factor de crecimiento $\rho = \frac{\max \lvert a^{(k)}_{ij}\rvert}{\max \lvert a_{ij}\rvert}$ | $\mathbf{2.641224 \times 10^{3}}$ | $1.000000$ | $\mathbf{1.783202}$ |
| Error relativo final $\lVert x - x^*\rVert_2 / \lVert x^*\rVert_2$ | $2.380178 \times 10^{-10}$ | $6.331029 \times 10^{-13}$ | $\mathbf{2.701713 \times 10^{-13}}$ |
| Error relativo máximo por componente | $7.5062 \times 10^{-8}$ ($x_6$) | $2.0007 \times 10^{-10}$ ($x_6$) | $\mathbf{8.5359 \times 10^{-11}}$ ($x_6$) |
| Bits de mantisa perdidos (norma 2, de 53 bits) | $\mathbf{20.03}$ **bits** | $11.48$ bits | $\mathbf{10.25}$ **bits** |
| Bits de mantisa perdidos en peor componente ($x_6$) | $\mathbf{28.33}$ **bits** | $19.78$ bits | $\mathbf{18.55}$ **bits** |
| Dígitos significativos correctos por incógnita | **6 a 9 dígitos** | 9 a 12 dígitos | **9 a 12 dígitos** |

**Por qué se descarta el primer pivote.** En la corrida sin pivoteo, el multiplicador de $F_6$ alcanza $2.752593\times10^{13}$. En la columna 6, el producto que se resta supera al valor original por un factor de 4832, aproximadamente $2^{12.24}$. Al alinear los exponentes se pierden bits bajos. Las filas quedan dominadas por múltiplos de $F_1$ y las restas posteriores entre términos cercanos agravan la cancelación. En $x_6$ se estiman 28.33 bits perdidos.

**Efecto del escalado.** El pivoteo parcial simple empieza por $F_6$, que tiene la mayor entrada absoluta en la primera columna. El escalado elige $F_8$, cuya razón vale 1. Aunque los multiplicadores brutos pueden ser grandes por la diferencia de escalas, el mayor multiplicador ajustado por las escalas es 0.962129. Frente a la corrida sin pivoteo, el error relativo baja por un factor de 881 y se recuperan unos 9.8 bits.

## 6. Condicionamiento

La estimación por valores singulares da $\kappa_2(A)\approx3.084744\times10^{13}$. Con $\varepsilon=2^{-52}\approx2.22\times10^{-16}$, la cota clásica de peor caso es $\kappa_2(A)\varepsilon\approx6.8495\times10^{-3}$. Bajo el criterio de cifras de la sección 7, esta estimación solo respalda **un dígito** a priori.

El pivoteo escalado alcanza entre 9 y 12 dígitos antes de refinar. Parte de la diferencia se explica por las unidades de las filas: al dividir cada fila y su término independiente por $s_i$ ($D=\operatorname{diag}(s_1,\dots,s_{12})$), la solución se conserva y el número de condición efectivo cae a $\kappa_2(D^{-1}A)\approx 8.132348\times10^3$. Con ese valor, la cota escalada es $\kappa_2(D^{-1}A)\varepsilon\approx 1.8057\times10^{-12}$, que bajo $\lfloor-\log_{10}(2e)\rfloor$ respalda alrededor de **11 dígitos significativos**, coherente con las 9 a 12 cifras obtenidas de forma directa. Además, la cota supone que el error actúa en la dirección de mayor amplificación, mientras que un residuo concreto puede proyectarse muy poco sobre esa dirección.

## 7. Resultados

Se defiende la solución refinada $x^{(1)}$. Para cada componente se compara el valor calculado con la referencia $x^*$ obtenida con Decimal a 50 cifras. Si $e_i=\lvert x_i-x_i^*\rvert/\lvert x_i^*\rvert$, se declaran $d_i=\lfloor-\log_{10}(2e_i)\rfloor$ cifras, con un máximo de 16 en doble precisión.

| Incógnita | Resultado directo | Cifras | Resultado refinado | Cifras defendidas |
| :--- | ---: | ---: | ---: | ---: |
| $x_{1}$ | 7.38491562700e-01 | 12 | 7.38491562700000e-01 | 15 |
| $x_{2}$ | -1.2059374168e+00 | 11 | -1.20593741680000e+00 | 15 |
| $x_{3}$ | 4.38719520e-03 | 9 | 4.387195200000000e-03 | 16 |
| $x_{4}$ | -9.1726453800e-02 | 11 | -9.172645380000000e-02 | 16 |
| $x_{5}$ | 2.84610937250e+00 | 12 | 2.846109372500000e+00 | 16 |
| $x_{6}$ | -6.28154090e-03 | 9 | -6.281540900000000e-03 | 16 |
| $x_{7}$ | 1.47382059160e+00 | 12 | 1.473820591600000e+00 | 16 |
| $x_{8}$ | -3.8492716400e-02 | 11 | -3.849271640000000e-02 | 16 |
| $x_{9}$ | 9.15264803700e-01 | 12 | 9.152648037000000e-01 | 16 |
| $x_{10}$ | -2.60713948520e+00 | 12 | -2.60713948520000e+00 | 15 |
| $x_{11}$ | 5.1873629400e-02 | 11 | 5.18736294000000e-02 | 15 |
| $x_{12}$ | -1.08493716250e+00 | 12 | -1.08493716250000e+00 | 15 |

En la solución directa, $x_3$ y $x_6$ son las más sensibles en términos relativos. Sus magnitudes son del orden de $10^{-3}$ y ambas conservan 9 cifras. Después de refinar, las 12 componentes coinciden con el redondeo de la referencia a float64. Los errores relativos por incógnita se guardan en resultados/solucion.csv.

## 8. Verificación

Se calcula $r=b-Ax$ con decimal.Decimal a 50 cifras a partir de los datos originales. Esto permite medir residuos de orden $10^{-10}$ en filas cuyos productos alcanzan $10^6$.

| Método evaluado | Norma infinito $\lVert r\rVert_\infty$ | Norma euclidiana $\lVert r\rVert_2$ | Error relativo $\lVert x - x^*\rVert_2 / \lVert x^*\rVert_2$ |
| :--- | :---: | :---: | :---: |
| Eliminación sin pivoteo | $1.771752003123 \times 10^{-7}$ | $2.172195951037 \times 10^{-7}$ | $2.380178 \times 10^{-10}$ |
| Parcial escalado ($x^{(0)}$) | $\mathbf{2.290876438498 \times 10^{-10}}$ | $2.369583447014 \times 10^{-10}$ | $\mathbf{2.701713 \times 10^{-13}}$ |
| Tras refinar ($x^{(1)}$) | $2.603187369605 \times 10^{-10}$ | $2.605075943953 \times 10^{-10}$ | $\mathbf{4.855364 \times 10^{-17}}$ |
| Control con NumPy | $7.756086660965 \times 10^{-11}$ | $1.095217213724 \times 10^{-10}$ | $2.926926 \times 10^{-13}$ |

Para refinar se resuelve $LU\,\delta x=P r^{(0)}$ con los factores de la corrida escalada y se actualiza $x^{(1)}=x^{(0)}+\delta x$. La corrección cuesta 300 flops en las dos sustituciones y 12 sumas para actualizar $x$, además del cálculo del residuo.

El error relativo de la solución disminuye de $2.70\times10^{-13}$ a $4.86\times10^{-17}$, mientras que la norma infinito del residuo sube de $2.29\times10^{-10}$ a $2.60\times10^{-10}$. La fila 6 determina esa norma. Sus coeficientes de orden $10^6$ amplifican pequeñas diferencias del redondeo de $x$ a doble precisión. En la fila 7, en cambio, el residuo absoluto baja de $1.39\times10^{-17}$ a $9.81\times10^{-19}$. Como control independiente, numpy.linalg.solve da un error relativo de $2.93\times10^{-13}$, del mismo orden que la resolución directa con pivoteo escalado. El residuo por fila y las correcciones se guardan en resultados/residuo.txt.

## 9. Conclusiones

Para este sistema, elegir el pivote según la escala de su fila reduce el crecimiento de 2641.22 a 1.78 con el mismo costo aritmético de eliminación y sustitución. El sobrecosto de seleccionar pivotes permite pasar de 6–9 a 9–12 cifras correctas. Reutilizar los factores en una pasada de refinamiento lleva el resultado a 15–16 cifras.

Con $n=1000$, el término dominante ya ronda $6.67\times10^8$ flops. Usaría una LU por bloques (BLAS de nivel 3) para aprovechar la caché en las actualizaciones de la matriz. Equilibraría filas y columnas con potencias de dos y reutilizaría los factores con el estimador de Hager–Higham para aproximar $\kappa_1(A)$ en $O(n^2)$ después de factorizar.

## 10. Reproducibilidad

Entorno probado en Windows de 64 bits con Python 3.14.4 y NumPy 2.5.2. El archivo requirements.txt fija NumPy 2.2.6 para Python 3.10–3.11 y 2.5.2 para Python 3.12 o superior. Desde la raíz del repositorio se ejecutan, en orden:

```bash
pip install -r requirements.txt
python src/main.py
```

**Uso de IA.** Se utilizó un asistente de IA para transcribir los valores de la matriz $A$ y el vector $b$ desde el enunciado hacia los archivos CSV y para revisar errores ortográficos en el informe.

**Referencias.** (1) Higham, N. J. (2002), *Accuracy and Stability of Numerical Algorithms* (2.ª ed.), SIAM, caps. 9, 12 y 15; (2) Golub, G. H. & Van Loan, C. F. (2013), *Matrix Computations* (4.ª ed.), Johns Hopkins University Press, cap. 3; (3) Chapra, S. C. & Canale, R. P. (2015), *Métodos Numéricos para Ingenieros* (7.ª ed.), McGraw-Hill.
