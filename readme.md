# El pivote · Eliminación gaussiana densa

Solución de un sistema lineal denso de 12 incógnitas con mínimo costo operacional (1486 flops), pivoteo parcial escalado y una pasada de refinamiento iterativo.

## Autoría y asignatura

- **Estudiante:** Ernesto Ponce
- **Asignatura:** Computación Numérica (INF-321) — Ingeniería Civil Informática, Universidad Católica del Maule
- **Docentes:** Juan Manuel Cabezas Contardo · Pedro Pesante
- **Fecha:** Talca, noviembre de 2026

## Requisitos

- **Python:** 3.10 o superior (probado en Python 3.14.4).
- **Dependencias:** requirements.txt fija numpy==2.2.6 para Python 3.10–3.11 y numpy==2.5.2 para Python 3.12 o superior.

```bash
pip install -r requirements.txt
```

## Ejecución

Desde la raíz del repositorio, ejecuta:

```bash
python src/main.py
```

*(En Windows con el lanzador de Python instalado también puede ejecutarse con py src/main.py.)*

El comando genera automáticamente cuatro archivos en resultados/: solucion.csv (incógnitas con sus dígitos defendidos), traza_eliminacion.txt (matriz y vector tras cada uno de los 11 pasos), conteo_flops.md (operaciones paso a paso) y residuo.txt (residuo en precisión extendida y refinamiento).

## Estructura de carpetas

```text
cn-parcial-gauss12-ponce/
├── datos/        Matriz A y vector b
├── src/          Eliminación, conteo y verificación
└── resultados/   Archivos generados por main.py
```

## Resultado resumido

Solución refinada a 10 dígitos significativos.

| Incógnita | Valor | Incógnita | Valor |
| :---: | :---: | :---: | :---: |
| $x_1$ | 0.7384915627 | $x_7$ | 1.473820592 |
| $x_2$ | -1.205937417 | $x_8$ | -0.03849271640 |
| $x_3$ | 0.004387195200 | $x_9$ | 0.9152648037 |
| $x_4$ | -0.09172645380 | $x_{10}$ | -2.607139485 |
| $x_5$ | 2.846109372 | $x_{11}$ | 0.05187362940 |
| $x_6$ | -0.006281540900 | $x_{12}$ | -1.084937162 |

**Costo de la resolución directa:** 1486 flops. **Norma infinito del residuo refinado:** 2.603187369605e-10.

## Informe técnico

Para el análisis completo de costo operacional, pérdida de mantisa, condicionamiento y verificación, ver [informe.md](informe.md).
