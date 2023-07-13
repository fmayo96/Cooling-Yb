# Cooling-Yb

Cambié todas las unidades para que no haya que hacer tantas conversiones. Ahora todo está en cm,s,J,K. 

Todas las magnitudes que definen ellos en el paper están en las unidades que corresponden. Únicamente modifiqué la pump intensity para que esté en las unidades correctas. 

El código que hay que ejecutar es el de la carpeta C_module. Se compila con el comando ```make``` y se ejecuta con ``` ./main.e```.

El programa resuelve la ecuación maestra que usan ellos en el paper integrando numéricamente con un método Runge-Kutta 4, y se obtienen como resultados el estado inicial y final del sistema, que es todo lo necesario para calcular las diferentes contribuciones a la potencia neta. 

En el archivo aux definimos algunas funciones auxiliares, como el cálculo de la frecuencia de Rabi, el estado térmico y la distribución de Bose-Einstein. 

En time_evol.c está definida la ecuación maestra y el método de integración. 

Todos los parámetros del problema se definen y se modifican en en archivo main.c. Una vez que estemos seguros de que esto esté funcionando bien podemos empezar a armar los distintos barridos. 
