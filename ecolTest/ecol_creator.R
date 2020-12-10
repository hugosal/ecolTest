s# parte 1. Instlación del paquete para prueba
setwd("C:/Users/user/Desktop/R/THutcheson/Hutcheson-T-test/trunk/ecolTest") # dar la dirección donde  se descargó el fichero del paquete "ecoltest"
library("devtools")
setwd("./ecolTest/")
setwd("..")
install("ecolTest")
library("ecolTest")
# como verificación de que se instaló correctamente el paquete prueba, comprobar que se carque la ayuda: 
help(Hutcheson.T.test)
help(multiple.Hutcheson.T.test)
#=============================================================
# parte 2.  PRUEBA CON DATOS REALES
#   2.1. Prueba para la función Hutcheson-T-test. Esta función calcula el p-val de Hutcheson T para los pares VAij
# cargar el archivo de prueba: prueba1.csv
# correr la función comparando por pares las VA's, ejemplo de sintaxis
# Hutcheson.T.test(data.frame$VAi,data.frame$VAj,alternative = "g, l, default")
#-------------------------------------------------------------
#   parte 2.2. Hutcheson T test múltiple.
# NOTA: el data frame debe de contener solo objestos numéricos
#ejemplo de sintaxis
# multiple.Hutcheson.T.test(nombre de data.frame con solo objetos numéricos)
