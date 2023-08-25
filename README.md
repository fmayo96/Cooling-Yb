# Cooling-Yb

Cambié todas las unidades para que no haya que hacer tantas conversiones. Ahora todo está en cm,s,J,K. 

Todas las magnitudes que definen ellos en el paper están en las unidades que corresponden. Únicamente modifiqué la pump intensity para que esté en las unidades correctas. 

Hay dos programas distintos para calcular las cosas. 

# Cómo ejecutar.
El código que funciona está metido en la carpeta Ivanov. El programa que hay que ejecutar es ```steady_state.py```. Es programa resuelve un ciclo de cooling para un set de parámetros (las intensidades de los láseres y las temperaturas) y genera un archivo txt con la diagonal del estado final, que es lo que se necesita para calcular la potencia. Después en la misma carpeta está el archivo ```first_results.py``` que agarra los datos de los txt y hace las curvas de potencia neta. 

Para ejecutar en segundo plano y poder dejar varias instancias corriendo lo que hay que poner en la terminal es ``` nohup python steady_state.py &> salida &```. Las librerías necesarias para ejecutar son numpy y numba. 
