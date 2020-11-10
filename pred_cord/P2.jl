using CSVFiles, JuMP, Pkg, Ipopt, DataFrames   #Importamos los paquetes a ocupar
const MOI = JuMP.MathOptInterface
#"coord_X","coord_Y" , types=[Int8,Float64,Float64,Int64]
DATA = DataFrame(load("pred_cord.csv", delim=';'))    #Lectura del CSV
DATA = convert(Matrix, DATA)
COOR = DATA[:,2:3]              #Guardamos las columnas con las coordenadas como una matriz
PROD = DATA[:,4]/10000         #Guardamos la producción de todas las plantas como un vector con valores de 1 a 3
PROD = convert.(Int64,PROD)     #Se definen los elementos del vector como números enteros
costos = [250,300,350]          #Se define el vector de costos de transporte
N = size(COOR)[1]               #Guardamos el número de plantas

#PREGUNTA 2.1:
function primer_centro()            #Función para encontrar la ubicación del primer centro de acopio
    modelo = Model(with_optimizer(Ipopt.Optimizer, print_level=0))      #Creamos el modelo y le asignamos el solver Ipopt

    @variable(modelo,0<=x)          #Creamos dos variables no negativas
    @variable(modelo,0<=y)

    @NLobjective(modelo,Min,sum(sqrt((COOR[i,1]-x)^2+(COOR[i,2]-y)^2)*1.1*costos[PROD[i]]*2 for i=1:N))     #Definimos la función objetivo no lineal
    optimize!(modelo)               #Optimizamos el modelo, busca minimizar el costo de ir a todas las plantas desde las coordenadas x e y

    coordx = round(value(x),digits = 4)
    coordy = round(value(y),digits = 4)                 #Guardamos la solución y el valor objetivo en nuevas variables
    valor = round(objective_value(modelo),digits = 0)

    println("-Las coordenadas que minimizan el costo, son X e Y iguales a $coordx y $coordy respectivamente para el primer centro de acopio.\n")
    println("-El valor de la funcion objetivo es US $valor\n")      #Mostramos los resultados en el terminal
    return coordx, coordy, valor       #La función devuelve los valores obtenidos
end

#Pregunta 2.2:
function segundo_centro(coordx, coordy, valor)      #Definimos la función para buscar la ubicación del segundo centro de acopio

    modelo2 = Model(with_optimizer(Ipopt.Optimizer, print_level=0))     #Definimos otro modelo con el mismo solver

    costoTotal(a,b) = sum(menor(i,a,b)*costos[PROD[i]]*2 for i=1:N)       #Definimos otra función para calcular el costo total de visitar todas las plantas teniendo los dos centros disponibles

    function menor(d,e,f)                       #Definimos otra función para determinar la distancia entre la planta d y el centro de acopio más cercano
        if(sqrt((COOR[d,1]-e)^2+(COOR[d,2]-f)^2)*1.1 < sqrt((COOR[d,1]-coordx)^2+(COOR[d,2]-coordy)^2)*1.1)
            return sqrt((COOR[d,1]-e)^2+(COOR[d,2]-f)^2)*1.1        #Si el nuevo centro está más cerca de esta planta, devuelve esta distancia
        else
            return sqrt((COOR[d,1]-coordx)^2+(COOR[d,2]-coordy)^2)*1.1  #Sino, devuelve la distancia al primer centro
        end
    end
    register(modelo2, :costoTotal, 2, costoTotal, autodiff=true)    #Registramos las dos funciones reién definidas para poder ocuparlas en el modelo
    register(modelo2, :menor, 3, menor, autodiff=true)

    @variable(modelo2,0<=x2)      #Definimos las variables x2 e y2 que determinan las coordenadas de la nueva planta
    @variable(modelo2,0<=y2)

    @NLobjective(modelo2,Min,costoTotal(x2,y2))       #Se define la función objetivo no lineal para el modeolo
    optimize!(modelo2)                      #Se optimiza el modelo minimizando el costo total de visitar todas las plantas teniendo los dos centros

    coordx2 = round(value(x2), digits=4)
    coordy2 = round(value(y2), digits=4)                #Guardamos la solución y el valor objetivo del modelo en nuevas variables
    valor2 = round(objective_value(modelo2), digits=0)
    reduccion = valor - valor2
    println("-Las coordenadas que minimizan el costo, son X e Y iguales a $coordx2 y $coordy2 respectivamente para el segundo centro.\n")
    println("-El valor de la funcion objetivo es US $valor2\n")              #Se muestran los resultados en el terminal
    println(reduccion > 400000 ? "-Vale la pena invertir en una segundo centro\n" : "-No vale la pena invertir en un segundo centro\n")    #Se recomienda construir solo si la reduccion es mayor a 400000
    println("-El costo total se redujo en US $reduccion")
end

coordx,coordy,valor = primer_centro()
segundo_centro(coordx,coordy,valor)
