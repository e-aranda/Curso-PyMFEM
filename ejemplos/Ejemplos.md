## Listado de Ejemplos

* Ejemplo 1: Laplace con condiciones Dirichlet homogéneas

---

* Laplace-Dirichlet: Laplace con condiciones Dirichlet no homogéneas, dadas por dos funciones distintas.

  En paralelo: medición de tiempo
---

* Conductor: Conducción del calor con coeficiente de difusión dependiente del dominio y temperatura fija en la frontera (dos valores). El coeficiente dependiente del dominio se puede obtener mediante `PyCoefficient`, `RestrictedCoefficient` o `PWCoefficient` (los dos últimos más eficientes)

  En paralelo, se puede ver la comparativa de tiempo con las diferentes definiciones del coeficiente de difusión.
---

* Elasticidad2D: Elasticidad 2D con fuerza aplicada en todo el dominio, condición Dirichlet nula en una parte de la frontera. Refinamiento de mallas. Visualización de mallas deformadas, extracción de componentes y interpolación de `GridFunction` entre espacios.

---

* Elasticidad3D: Elasticidad 3D sin fuerza en el dominio, con condición Dirichlet nula en un parte de la frontera y con fuerza aplicada en otra parte de la frontera. Cálculo de la tensión de Von Misses a partir de los desplazamientos (`Coefficient` calculado mediante `GridFunction`). Obtención de la solución en puntos de la malla (localización con `FindPoints` y extracción de valores). Evaluación de integrales. Dependencia de los coeficientes calculado con respecto a la `GridFunction` asociada.

En el paralelo, atención al uso de `FindPoints`. Uso de MPIAllreduce para el cálculo de integrales.

---

* Termico: Problema parabólico: conducción térmica. Ecuación del calor con coeficiente depediente del dominio, condición Dirichlet no nula en una parte de la frontera y condición Robin en otra. Resuelto mediante discretización en tiempo, y bucle implícito. Se modifican los atributos de la malla. Visualización de todos los instantes temporales

---

* Ejemplo 7: Problema de Darcy: formulación variacional mixta en espacio Raviart-Thomas. Sin codición Dirichlet sobre `u`, la condición frontera Dirichlet de `p` resulta en condición natural. Se usa `MINRES` como *solver*. Se compara con la solución exacta (`ComputeL2Error` y `ComputeLpNorm`).

---

* Termico-C: Mismo problema parabólico de conducción térmica resuelto usando un `TimeDependentOperator`.

La paralelización exige cambios significativos

---

* Ejemplo 9: Problema de Stokes: formulación variacional mixta en (H2,H1) con condiciones Dirichlet no nulas en una parte de la frontera. Variante del ej 7.

La convergencia en este ejemplo es más complicada. Por alguna razón, es preciso multiplicar por (-1) la ecuación de divergencia. Solo `MINRES` logra convergencia.

---

* Ejemplo 10: Problema no lineal: difusión con término adicional no lineal.

---

* Ejemplo 11: Problema no lineal: difusión con coeficiente no lineal. Idéntico al anterior, salvo que el integrador no lineal es más complicado

---

* Ejemplo 12: Ecuación de advección (transporte): ecuación de primer orden sin condiciones frontera, sobre una malla periódica. La formulación variacional usa DG (Discontinuous Galerking) Se usa un `TimeDependentOperator`. Atención a la convergencia con métodos explícitos (necesario paso de tiempo pequeño)
  
---

* Ejemplo 13: Ecuación de Laplace con condiciones de frontera periódicas. La periodicidad está exclusivamente en la malla considerada.

---

* Ejemplo 14: Ecuación de Laplace mediante formulación variacional mixta con elementos RT (similar al ej 7)

---

* Ejemplo 15: Bilaplaciano con condiciones frontera sobre $u$ y $\Delta u$. Resuelto de forma desacoplada, lo que introduce una `GridFunctionCoefficient`

---

* Ejemplo 16: Igual al ej.16 pero resuelto con matriz formada por bloques.

---

* Ejemplo 17: Laplace con condiciones Neumann en toda la frontera. El problema está mal planteado y requiere una restricción adicional para garantizar existencia $\int_\Omega u\,dx = 0$. Se usa `BlockOperator`. La extensión al paralelo no funciona pues las submatrices resultantes son singulares. El ej17p solo da resultado correcto con 1 procesador. No obstante, es interesante ver cómo se ha realizado el paralelo.

--- 

* Ejemplo 18: Mismo problema que ej17, pero corrigiendo el probleam usando un `OrthoSolver`.