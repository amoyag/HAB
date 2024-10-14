
## Rúbrica para evaluación de tareas de análisis funcional

1.	Funcionalidad del Script (40%)	(30% ponderado)
	-	El script cumple con los requisitos solicitados (análisis funcional de una lista de genes, exportación del archivo con los resultados, posibilidad de ordenar por FDR, configuración de categorías funcionales...)
	
2.	Eficiencia y Organización del Código (20%)	(15% ponderado)
	- EL código es eficiente en cuanto a tiempos de ejecución, uso de bibliotecas adecuadas, y claridad en la lógica implementada.
	-	El código está bien organizado y sigue buenas prácticas de programación (uso de funciones, modularización, etc.).
3.	Flexibilidad y Configurabilidad (15%)	(11% ponderado)
	-	El script es flexible y permite la configuración adecuada de parámetros, como el umbral de FDR o la categoría funcional a elegir.
4.	Manejo de Errores y Validación (15%)	(11% ponderado)
	- 	El script maneja errores comunes (como entradas incorrectas o falta de conexión a bases de datos) de manera adecuada, con mensajes de error útiles para el usuario.
	- 	El código incluye validaciones para las entradas del usuario (ver explicación).
5.	Documentación y Comentarios (10%)	(8% ponderado)
	-	El código está bien documentado, incluyendo comentarios que expliquen las secciones clave del código.
	-	Se ha proporcionado una descripción clara de cómo usar el script, ya sea dentro del código o como un archivo separado.

### Ejercicio práctico
El script deberá usarse pararesponder a una pregunta. La respuesta correcta es el 25% de la calificación.
### Validaciones para la entrada de usuario
“Validaciones para la entrada del usuario”, se refiere a verificar si el script maneja adecuadamente las entradas que recibe del usuario para evitar errores o problemas de ejecución:
- **Validación de tipo de datos:** el script verifica que el tipo de datos introducidos por el usuario es el correcto. Por ejemplo, si se espera un número (como el valor de FDR o un identificador de especie), el script debería asegurarse de que el usuario ingresa un número y no un texto.
- **Validación de formato:** Algunos parámetros pueden tener un formato específico. Por ejemplo, si el usuario ingresa una lista de genes, esta debería tener el formato adecuado (por ejemplo, separados por comas). El script debería validar esto y mostrar un mensaje de error claro si el formato no es correcto.
- **Verificación de rangos y valores permitidos:** Si el usuario introduce un valor numérico, como el umbral de FDR, el script debería asegurarse de que ese valor está dentro de un rango razonable (por ejemplo, entre 0 y 1). Si el valor está fuera de este rango, el script debería advertir al usuario.
- **Manejo de entradas faltantes:**
Si un parámetro es obligatorio (como el nombre de la especie o los genes a analizar), el script debería comprobar que el usuario lo ha proporcionado y mostrar un error si no lo ha hecho.
- **Mensajes de error claros:**
Cuando el usuario introduce una entrada incorrecta, el script debería generar mensajes de error claros que expliquen el problema y, si es posible, indiquen cómo corregirlo. Esto mejora la experiencia del usuario y ayuda a evitar confusiones.

Ejemplo de validación:

Validar que el valor de FDR es un número entre 0 y 1

```

fdr = input("Introduce el umbral de FDR (entre 0 y 1): ")
try:
    fdr = float(fdr)
    if fdr < 0 or fdr > 1:
        print("Error: El FDR debe estar entre 0 y 1.")
except ValueError:
    print("Error: Debes introducir un número.")
``` 

