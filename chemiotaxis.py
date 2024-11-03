import math
import random
import numpy
import copy

from bacteria import bacteria


class chemiotaxis():
    def __init__(self):
       parcialNFE = 0  
    
    def compute_cell_interaction(self, bacteria, poblacion, d, w):
      total = 0.0
      for other in poblacion:
        diff = 0.0
        diff += (bacteria.blosumScore - other.blosumScore) ** 2.0
        total += d * math.exp(w * diff)
        #print("diff: ", diff, "total: ", total)
      return total

    def attract_repel(self, bacteria, poblacion, d_attr, w_attr, h_rep, w_rep):
      attract = self.compute_cell_interaction(bacteria, poblacion,  -d_attr, -w_attr)
      repel = self.compute_cell_interaction(bacteria, poblacion,  h_rep, -w_rep)
      #print("attract: ", attract, "repel: ", repel)
      return attract + repel   #interaction
    
    
    def chemio(self, bacteria, poblacion, d_attr, w_attr, h_rep, w_rep):
      bacteria.interaction = self.attract_repel(bacteria, poblacion, d_attr, w_attr, h_rep, w_rep)
      bacteria.fitness = bacteria.blosumScore + bacteria.interaction
      

    
    def doChemioTaxis(self, poblacion, d_attr, w_attr, h_rep, w_rep):
      self.parcialNFE = 0
      for bacteria in poblacion:
        self.chemio(bacteria, poblacion, d_attr, w_attr, h_rep, w_rep)
        self.parcialNFE += bacteria.NFE
        bacteria.NFE = 0
        

    def eliminarClonar(self, path, poblacion):
      """elimina el 50% de las bacterias con menos fitness """
      poblacion.sort(key=lambda x: x.fitness)
      for i in range(int(len(poblacion)/2)):
          del poblacion[0]
          
      clones = self.clonacion(path, poblacion)
      for clone in clones:
            poblacion.append(clone)
    
    def clonacion(self, path, poblacion):
       poblacionClones = []
       best = max(poblacion, key=lambda x: x.fitness)
       for bacteria in poblacion:
         newBacteria = bacteria.clonar(path)
         mutacion = int((best.fitness - bacteria.fitness)/10)    #mutacion en funcion de la diferencia de fitness
         newBacteria.tumboNado(mutacion)
         newBacteria.autoEvalua()
         poblacionClones.append(newBacteria)
       return poblacionClones

    
         
    def randomBacteria(self, path):
       bact = bacteria(path)
       bact.tumboNado(random.randint(1, 10))
       return bact 
   
    def insertRamdomBacterias(self, path, num, poblacion):
      for i in range(num):
         poblacion.append(self.randomBacteria(path))
         #eliminar la bacteria con menos fitness
         poblacion.sort(key=lambda x: x.fitness)
         del poblacion[0]
    
    def mutacionAdaptativa(self, path, poblacion, max_mutacion=5):
        
        # Encuentra el mejor y el peor fitness en la población
        best_fitness = max(bacteria.fitness for bacteria in poblacion)
        worst_fitness = min(bacteria.fitness for bacteria in poblacion)

        for bact in poblacion:
            # Calcula la intensidad de mutación adaptativa basada en el fitness
            if best_fitness != worst_fitness:  # Evita división por cero
                intensidad_mutacion = max_mutacion * (1 - (bact.fitness - worst_fitness) / (best_fitness - worst_fitness))
            else:
                intensidad_mutacion = max_mutacion / 2  # Si todos tienen el mismo fitness, usar una mutación media
            
            intensidad_mutacion = max(1, int(intensidad_mutacion))  # Asegura un mínimo de mutación

            # Aplica la mutación adaptativa
            bact.tumboNado(intensidad_mutacion)
            bact.autoEvalua()  # Reevaluar el fitness después de la mutación
    
    def balanceaDiversidad(self, poblacion, tolerancia=0.8):
      # Ordena bacterias por similitud en función de su fitness
        poblacion.sort(key=lambda x: x.fitness)
    
      # Identifica bacterias con alta similitud
        for i in range(len(poblacion)-1):
          if abs(poblacion[i].fitness - poblacion[i+1].fitness) < tolerancia:
            # Reemplaza la bacteria por una nueva y diversa
            nueva_bacteria = self.randomBacteria(poblacion[i].matrix.path)
            poblacion[i] = nueva_bacteria
    
    def mutaBacteriasBajas(self, poblacion, factor_intensidad=5):
        promedio_fitness = sum(b.fitness for b in poblacion) / len(poblacion)
    
        for bacteria in poblacion:
          if bacteria.fitness < promedio_fitness:
            # Mutación intensiva en bacterias con bajo fitness
            bacteria.tumboNado(factor_intensidad)
            bacteria.autoEvalua()
            
    def migracion(self, poblacion, porcentaje_migracion=0.2):
        
        # Ordena la población por fitness de menor a mayor
        poblacion.sort(key=lambda x: x.fitness)
        
        # Determina el número de bacterias a migrar
        num_migracion = int(len(poblacion) * porcentaje_migracion)
        
        # Selecciona las bacterias con menor y mayor fitness
        bacterias_bajas = poblacion[:num_migracion]
        bacterias_altas = poblacion[-num_migracion:]
        
        # Realiza la migración, reemplazando las bacterias con bajo fitness por clones de las de mayor fitness
        for baja, alta in zip(bacterias_bajas, bacterias_altas):
            baja.matrix.seqs = numpy.array(copy.deepcopy(alta.matrix.seqs))
            baja.autoEvalua()  # Reevaluar el fitness después de la migración        
    

# Path: BFOA_MSAv2/evaluadorBlosum.py