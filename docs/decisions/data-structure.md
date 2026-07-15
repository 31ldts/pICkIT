# Decisión: estructura de datos interna de la matriz de interacciones

**Fecha:** 2026-07-14 (durante la extracción de `models.py` y `filter_mixin.py`, Fase 2)
**Estado:** Decidido — se mantiene lista de listas (`list[list[str]]`)

## Contexto

El Anexo B del plan de acción planteaba evaluar `list[list[str]]` (estructura
actual) contra `pandas.DataFrame` o un `dict`/JSON anidado, con el criterio
explícito de no cambiar la estructura salvo que las operaciones reales lo
justificaran una vez el código estuviera aislado en `filter_mixin.py`.

## Qué se observó al aislar `filter_mixin.py`

- `transpose_matrix` es efectivamente trivial con listas anidadas (una
  comprensión de listas), tal como anticipaba el Anexo B.
- Las operaciones de filtrado (`filter_by_interaction`, `filter_by_subunit`,
  `filter_by_residue`) no indexan por nombre de residuo/ligando de forma
  repetida de un modo que un `DataFrame.loc` simplificara de verdad: lo que
  hacen es parsear el *contenido* de cada celda (una cadena codificada tipo
  `"13 |CA-O1(A)|; 8 |CB-C2(A)|"` con múltiples interacciones, átomos y
  subunidades empaquetados). Ese formato de celda es el verdadero punto de
  complejidad, no el acceso a filas/columnas.
- No se encontró duplicación significativa de "parsear el string de la celda
  una y otra vez" que por sí sola justificara introducir una clase `Cell`
  dedicada — cada método (`filter_by_interaction`, `filter_by_subunit`,
  `filter_by_residue`, `_stack_reactives`/`_get_interactions`,
  `heatmap`/`process_matrix`) parsea la celda con una necesidad ligeramente
  distinta (por tipo de interacción, por subunidad, por átomo/cadena
  principal-lateral, por conteo), y forzar una abstracción común en este
  punto habría sido un rediseño no justificado por el propio criterio del
  Anexo B.
- `pandas` sí se usa, pero solo en la frontera de exportación
  (`get_dataframe`, `save_interaction_data` vía `export_mixin.py`) y en
  `plot_mixin.py` para construir el `DataFrame` que consume `seaborn`. La
  matriz interna nunca necesita operaciones vectorizadas de pandas durante
  el filtrado/ordenado — solo al final, para presentación.

## Decisión

Se mantiene `InteractionData.matrix` como `list[list[str]]`, sin rediseño.
Revisar esta decisión si en el futuro se detecta que una parte
significativa de los bugs o de la complejidad de `filter_mixin.py` viene
específicamente de re-parsear el string de celda de forma repetida (el
criterio práctico que marcaba el propio Anexo B) — no se ha observado eso
en esta extracción.
