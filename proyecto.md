

## Rúbrica de Evaluación del Proyecto de Análisis Funcional y Propagación en Redes

Cada proyecto será evaluado sobre **100 puntos**, distribuidos en **6 criterios principales**. Se valorará tanto la **calidad técnica** como la **claridad del diseño y documentación**.

| Criterio | Descripción | Puntos |
|---------|-------------|--------|
| **1. Funcionalidad del script (core)** | El script debe ejecutar correctamente el análisis funcional y/o la propagación en redes sobre una lista de genes diferencialmente expresados. Se valorará que el código sea robusto, modular y que maneje errores. | 25 |
| **2. Elección y justificación de técnicas** | Se evaluará la elección de métodos (e.g., ORA, GSEA, RWR, DIAMOnD, etc.) y su adecuación al problema. Deben justificar por qué han elegido esas técnicas y no otras. | 15 |
| **3. Automatización y diseño del flujo de trabajo** | Se valorará que el script permita una ejecución sencilla (por CLI o interfaz), que automatice tareas como la conversión de IDs, descarga de datos, preprocesamiento, etc. | 15 |
| **4. Documentación y reproducibilidad** | El repositorio debe incluir un `README.md` claro con instrucciones de uso, dependencias, ejemplos de ejecución y explicación del flujo de trabajo. Se valorará el uso de `requirements.txt` o `environment.yml`. | 15 |
| **5. Calidad del código y buenas prácticas** | Código limpio, modular, con funciones bien definidas, comentarios útiles, uso adecuado de estructuras de datos y convenciones de estilo (PEP8). | 10 |
| **6. Análisis y visualización de resultados** | Se valorará que el script genere salidas interpretables (tablas, gráficos, informes) y que se incluyan ejemplos de resultados con interpretación biológica. | 10 |
| **Bonus: originalidad y extensión** | Uso de técnicas avanzadas, integración de múltiples fuentes de datos, visualizaciones interactivas, uso de APIs externas, etc. | +10 |

---

### Estructura mínima esperada del repositorio

```
/nombre-del-proyecto/
├── README.md
├── script_principal.py
├── utils/
│   └── conversion.py
├── data/
│   └── genes_input.txt
├── results/
│   └── output.tsv
├── requirements.txt
└── LICENSE (opcional)
```
