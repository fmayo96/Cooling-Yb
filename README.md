# Cooling-Yb

Cambié todas las unidades para que no haya que hacer tantas conversiones. Ahora todo está en cm,s,J,K. 

Todas las magnitudes que definen ellos en el paper están en las unidades que corresponden. Únicamente modifiqué la pump intensity para que esté en las unidades correctas. 

Hay dos programas distintos para calcular las cosas. 

## Código en C.

Uno es el código en C que está en la carpeta ```C_module```. Este código se compila con el comando ```make``` y se ejecuta con ``` ./main.e```.

El programa resuelve la ecuación maestra que usan ellos en el paper integrando numéricamente con un método Runge-Kutta 4, y se obtienen como resultados el estado inicial y final del sistema, que es todo lo necesario para calcular las diferentes contribuciones a la potencia neta. Cada ejecución del programa genera un archivo txt con el estado final del sistema, que es lo necesario para calcular la potencia de enfriamiento, esa parte la hacemos en python. ***Por ahora estamos teniendo algún error que hace que den mal las potencias. Ese error tiene que estar en el archivo ```time_evol.c``` donde definimos la ecuación diferencial***. Por ahora está chequeado que el estado del sistema es hermítico, definido positivo y mantiene la traza igual a uno (dentro de los errores numéricos).

En el archivo aux definimos algunas funciones auxiliares, como el cálculo de la frecuencia de Rabi, el estado térmico y la distribución de Bose-Einstein. 

En time_evol.c está definida la ecuación maestra y el método de integración. 

Todos los parámetros del problema se definen y se modifican en en archivo main.c. Una vez que estemos seguros de que esto esté funcionando bien podemos empezar a armar los distintos barridos. 

## Código en Python.

En la carpeta ```Ivanov``` está el código que resuelve calculando el estado estacionario, que no es estrictamente lo que hacen ellos (según dicen), pero es bastante parecido porque ellos integran hasta un tiempo muy largo. Este programa es rápido y permite calcular las curvas en un tiempo razonable. Por ahora está implementado con un único láser, y en processo de debuggeo. 
