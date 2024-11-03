from bacteria import bacteria
from chemiotaxis import chemiotaxis
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Configuración inicial
poblacion = []
path = "C:\\Proyecto Bacteria\\multiFasta.fasta"
numeroDeBacterias = 5
numRandomBacteria = 1
iteraciones = 50
tumbo = 1  # Número de gaps a insertar
nado = 3
chemio = chemiotaxis()
veryBest = bacteria(path)  # Mejor bacteria
tempBacteria = bacteria(path)  # Bacteria temporal para validaciones
original = bacteria(path)  # Bacteria original sin gaps
globalNFE = 0  # Número de evaluaciones de la función objetivo

dAttr = 0.1  # 0.1
wAttr = 0.2  # 0.2
hRep = dAttr
wRep = 10  # 10

# Inicialización de listas para recopilar datos
fitness_vals = []
nfe_vals = []
population_sizes = []
very_best_fitness = []

# Funciones auxiliares
def clonaBest(veryBest, best):
    veryBest.matrix.seqs = np.array(best.matrix.seqs)
    veryBest.blosumScore = best.blosumScore
    veryBest.fitness = best.fitness
    veryBest.interaction = best.interaction


def validaSecuencias(path, veryBest):
    # Clona a veryBest en tempBacteria
    tempBacteria.matrix.seqs = np.array(veryBest.matrix.seqs)
    # Descartar los gaps de cada secuencia
    for i in range(len(tempBacteria.matrix.seqs)):
        tempBacteria.matrix.seqs[i] = tempBacteria.matrix.seqs[i].replace("-", "")
    # Valida que las secuencias originales sean iguales a las secuencias de tempBacteria
    for i in range(len(tempBacteria.matrix.seqs)):
        if tempBacteria.matrix.seqs[i] != original.matrix.seqs[i]:
            print("*****************Secuencias no coinciden********************")
            return


# Inicialización de la población
for i in range(numeroDeBacterias):
    poblacion.append(bacteria(path))

# Ciclo de optimización
for iteracion in range(iteraciones):
    for bact in poblacion:
        bact.tumboNado(tumbo)
        bact.autoEvalua()

    chemio.doChemioTaxis(poblacion, dAttr, wAttr, hRep, wRep)
    globalNFE += chemio.parcialNFE

    # Guardar datos de cada iteración
    best = max(poblacion, key=lambda x: x.fitness)
    fitness_vals.append(best.fitness)
    nfe_vals.append(globalNFE)
    population_sizes.append(len(poblacion))
    very_best_fitness.append(veryBest.fitness)

    if (veryBest is None) or (best.fitness > veryBest.fitness):
        clonaBest(veryBest, best)
    print("Interacción: ", veryBest.interaction, "Fitness: ", veryBest.fitness, "NFE:", globalNFE)

    chemio.eliminarClonar(path, poblacion)
    chemio.insertRamdomBacterias(path, numRandomBacteria, poblacion)
    print("Población actual:", len(poblacion))

    # Llama al nuevo método mutacionDirigida
    for bacteria in poblacion:
        bacteria.mutacionDirigida(num_mutaciones=3, probabilidad_mutacion=0.5)

    # Aplica migración cada 10 iteraciones
    if iteracion % 10 == 0 and iteracion > 0:
        chemio.migracion(poblacion, porcentaje_migracion=0.2)
        print("Migración aplicada para aumentar la diversidad genética.")

# Finalización: mostrar el mejor genoma y validar secuencias
veryBest.showGenome()
validaSecuencias(path, veryBest)

# Crear tabla de datos recopilados
tabla_datos = pd.DataFrame({
    'Iteración': range(1, iteraciones + 1),
    'Fitness Mejor Bacteria': fitness_vals,
    'NFE': nfe_vals,
    'Tamaño Población': population_sizes,
    'Fitness VeryBest': very_best_fitness
})

# Muestra la tabla
print("\nTabla de datos recopilados:")
print(tabla_datos)

# Gráfica del Fitness de la Mejor Bacteria en Cada Iteración
plt.figure(figsize=(10, 6))
plt.plot(tabla_datos['Iteración'], tabla_datos['Fitness Mejor Bacteria'], label='Fitness Mejor Bacteria')
plt.plot(tabla_datos['Iteración'], tabla_datos['Fitness VeryBest'], label='Fitness VeryBest', linestyle='--')
plt.xlabel('Iteración')
plt.ylabel('Fitness')
plt.title('Progreso del Fitness de la Mejor Bacteria en Cada Iteración')
plt.legend()
plt.grid()
plt.show()

# Gráfica de NFE a lo Largo de las Iteraciones
plt.figure(figsize=(10, 6))
plt.plot(tabla_datos['Iteración'], tabla_datos['NFE'], color='purple')
plt.xlabel('Iteración')
plt.ylabel('Número de Evaluaciones de la Función Objetivo (NFE)')
plt.title('Progreso de NFE a lo Largo de las Iteraciones')
plt.grid()
plt.show()

# Gráfica del Tamaño de la Población
plt.figure(figsize=(10, 6))
plt.plot(tabla_datos['Iteración'], tabla_datos['Tamaño Población'], color='orange')
plt.xlabel('Iteración')
plt.ylabel('Tamaño de Población')
plt.title('Tamaño de la Población en Cada Iteración')
plt.grid()
plt.show()


