#Script TFM Master en Bioinformatica avanzado

# Javier Pineda Pampliega


# 1) Extraccion de taxonomia a partir de los resultados de Kraken2 --------


direct=choose.dir() # Hay que elegir la carpeta con todos los archivos output.

total_tax=NULL
grupo_estudiar="P" #P para filo, C para clase, O para orden.

for(pedazos in 1:length(list.files(direct))) ###Se extraen todos los nombre del grupo en cuestion
{
  files_read=read.delim(header = F,paste(direct,list.files(direct)[pedazos],sep="\\"))
  # paste es para unir vectores y luego convertirlos a caracter. Así, direct es la direccion de la carpeta, y list.files(direct) genera los nombres de todos los archivos.
  # Con este bucle, va abriendo los archivos 1 a 1, y los llama "files_read".
  
  for(i in 1:length(files_read[,1])) # Esta parte permite ignorar los valores que no sean bacteria, ya que los resultados de Kraken tambien incluyen eucariota.
  {
    if(files_read[i,6]=="  Bacteria")
    {
      start_bacteria=i  
    }else if (files_read[i,6]=="  Eukaryota")
    {
      end_bacteria=i
    }
  }
  files_read=files_read[start_bacteria:end_bacteria-1,]
  
  ###Filo
  
  files_filo = files_read[files_read[,4]==grupo_estudiar,] # Del archivo completo, crea otro donde la cuarta columna sea igual a la letra indicada arriba en "grupo_estudiar".
  
  total_tax = c(total_tax,files_filo[,6]) # Al archivo total_tax (vacio) le añades las sexta columna del archivo que acabas de leer.
}

total_tax=unique(gsub(" ", "", total_tax, fixed = TRUE))

Output=data.frame(matrix(nrow = length(total_tax)))

rownames(Output)=total_tax ##Copiamos los nombres a Output

for(pedazos in 1:length(list.files(direct)))
{
  files_read=read.delim(header = F,paste(direct,list.files(direct)[pedazos],sep="\\"))
  
  for(i in 1:length(files_read[,1]))
  {
    if(files_read[i,6]=="  Bacteria")
    {
      start_bacteria=i  
    }else if (files_read[i,6]=="  Eukaryota")
    {
      end_bacteria=i
    }
  }
  files_read=files_read[start_bacteria:end_bacteria-1,]
  
  #Filo
  
  files_filo=files_read[files_read[,4]==grupo_estudiar,] 
  
  files_filo[,6]=gsub(" ", "", files_filo[,6], fixed = TRUE)
  
  agg = aggregate(files_filo$V2,by = list(files_filo$V6),FUN = sum)
  
  Output[match(agg[,1],row.names(Output)),pedazos]=agg[,2]/sum(agg[,2]) 
  
  colnames(Output)[pedazos]=list.files(direct)[pedazos] 
}


colSums(Output,na.rm = T)

write.table(Output, file="Tabla_filo.txt")




# 2) Comparacion resultados amplicones y RNASeq ------------------------------

library(tidyverse) # Libreria para trabajar con datos
library(ggplot2) # Libreria para crear graficas

# Empleando RStudio, se cargan las tablas:Kraken_filo, Kraken_genero, QIIME2_filo_porcentaje y QIIME2_genero_porcentaje

# El siguiente paso es combinar las tablas (filo con filo y clase con clase)

tabla_filo = inner_join(QIIME2_filo_percentage, Kraken_filo,  by = "Filo")

tabla_genero = inner_join(QIIME2_genero_percentage, Kraken_genero,  by = "Genero")

# Ahora se evalua la correlacion muestra a muestra

Correlacion_CRC01_filo = cor.test(tabla_filo$CRC01,tabla_filo$CRC01_kraken)
Correlacion_CRC02_filo = cor.test(tabla_filo$CRC02,tabla_filo$CRC02_kraken)
Correlacion_CRC03_filo = cor.test(tabla_filo$CRC03,tabla_filo$CRC03_kraken)
Correlacion_CRC04_filo = cor.test(tabla_filo$CRC04,tabla_filo$CRC04_kraken)
Correlacion_CRC05_filo = cor.test(tabla_filo$CRC05,tabla_filo$CRC05_kraken)
Correlacion_CRC06_filo = cor.test(tabla_filo$CRC06,tabla_filo$CRC06_kraken)
Correlacion_CRC07_filo = cor.test(tabla_filo$CRC07,tabla_filo$CRC07_kraken)
Correlacion_CRC08_filo = cor.test(tabla_filo$CRC08,tabla_filo$CRC08_kraken)
Correlacion_CRC09_filo = cor.test(tabla_filo$CRC09,tabla_filo$CRC09_kraken)
Correlacion_CRC10_filo = cor.test(tabla_filo$CRC10,tabla_filo$CRC10_kraken)
Correlacion_CRC11_filo = cor.test(tabla_filo$CRC11,tabla_filo$CRC11_kraken)
Correlacion_CRC12_filo = cor.test(tabla_filo$CRC12,tabla_filo$CRC12_kraken)
Correlacion_CRC13_filo = cor.test(tabla_filo$CRC13,tabla_filo$CRC13_kraken)
Correlacion_CRC15_filo = cor.test(tabla_filo$CRC15,tabla_filo$CRC15_kraken)
Correlacion_CRC16_filo = cor.test(tabla_filo$CRC16,tabla_filo$CRC16_kraken)
Correlacion_CRC17_filo = cor.test(tabla_filo$CRC17,tabla_filo$CRC17_kraken)
Correlacion_CRC18_filo = cor.test(tabla_filo$CRC18,tabla_filo$CRC18_kraken)
Correlacion_CRC19_filo = cor.test(tabla_filo$CRC19,tabla_filo$CRC19_kraken)
Correlacion_CRC20_filo = cor.test(tabla_filo$CRC20,tabla_filo$CRC20_kraken)
Correlacion_CRC21_filo = cor.test(tabla_filo$CRC21,tabla_filo$CRC21_kraken)
Correlacion_CRC22_filo = cor.test(tabla_filo$CRC22,tabla_filo$CRC22_kraken)
Correlacion_CRC23_filo = cor.test(tabla_filo$CRC23,tabla_filo$CRC23_kraken)
Correlacion_CRC24_filo = cor.test(tabla_filo$CRC24,tabla_filo$CRC24_kraken)
Correlacion_CRC25_filo = cor.test(tabla_filo$CRC25,tabla_filo$CRC25_kraken)
Correlacion_CRC26_filo = cor.test(tabla_filo$CRC26,tabla_filo$CRC26_kraken)
Correlacion_CRC27_filo = cor.test(tabla_filo$CRC27,tabla_filo$CRC27_kraken)
Correlacion_CRC28_filo = cor.test(tabla_filo$CRC28,tabla_filo$CRC28_kraken)
Correlacion_CRC29_filo = cor.test(tabla_filo$CRC29,tabla_filo$CRC29_kraken)
Correlacion_CRC30_filo = cor.test(tabla_filo$CRC30,tabla_filo$CRC30_kraken)
Correlacion_CRC31_filo = cor.test(tabla_filo$CRC31,tabla_filo$CRC31_kraken)
Correlacion_CRC32_filo = cor.test(tabla_filo$CRC32,tabla_filo$CRC32_kraken)
Correlacion_CRC33_filo = cor.test(tabla_filo$CRC33,tabla_filo$CRC33_kraken)
Correlacion_CRC34_filo = cor.test(tabla_filo$CRC34,tabla_filo$CRC34_kraken)

Correlacion_CRC01_filo
Correlacion_CRC02_filo
Correlacion_CRC03_filo
Correlacion_CRC04_filo
Correlacion_CRC05_filo
Correlacion_CRC06_filo
Correlacion_CRC07_filo
Correlacion_CRC08_filo
Correlacion_CRC09_filo
Correlacion_CRC10_filo
Correlacion_CRC11_filo
Correlacion_CRC12_filo
Correlacion_CRC13_filo
Correlacion_CRC15_filo
Correlacion_CRC16_filo
Correlacion_CRC17_filo
Correlacion_CRC18_filo
Correlacion_CRC19_filo
Correlacion_CRC20_filo
Correlacion_CRC21_filo
Correlacion_CRC22_filo
Correlacion_CRC23_filo
Correlacion_CRC24_filo
Correlacion_CRC25_filo
Correlacion_CRC26_filo
Correlacion_CRC27_filo
Correlacion_CRC28_filo
Correlacion_CRC29_filo
Correlacion_CRC30_filo
Correlacion_CRC31_filo
Correlacion_CRC32_filo
Correlacion_CRC33_filo
Correlacion_CRC34_filo


# Genero ------------------------------------------------------------------


Correlacion_CRC01_genero = cor.test(tabla_genero$CRC01,tabla_genero$CRC01_kraken)
Correlacion_CRC02_genero = cor.test(tabla_genero$CRC02,tabla_genero$CRC02_kraken)
Correlacion_CRC03_genero = cor.test(tabla_genero$CRC03,tabla_genero$CRC03_kraken)
Correlacion_CRC04_genero = cor.test(tabla_genero$CRC04,tabla_genero$CRC04_kraken)
Correlacion_CRC05_genero = cor.test(tabla_genero$CRC05,tabla_genero$CRC05_kraken)
Correlacion_CRC06_genero = cor.test(tabla_genero$CRC06,tabla_genero$CRC06_kraken)
Correlacion_CRC07_genero = cor.test(tabla_genero$CRC07,tabla_genero$CRC07_kraken)
Correlacion_CRC08_genero = cor.test(tabla_genero$CRC08,tabla_genero$CRC08_kraken)
Correlacion_CRC09_genero = cor.test(tabla_genero$CRC09,tabla_genero$CRC09_kraken)
Correlacion_CRC10_genero = cor.test(tabla_genero$CRC10,tabla_genero$CRC10_kraken)
Correlacion_CRC11_genero = cor.test(tabla_genero$CRC11,tabla_genero$CRC11_kraken)
Correlacion_CRC12_genero = cor.test(tabla_genero$CRC12,tabla_genero$CRC12_kraken)
Correlacion_CRC13_genero = cor.test(tabla_genero$CRC13,tabla_genero$CRC13_kraken)
Correlacion_CRC15_genero = cor.test(tabla_genero$CRC15,tabla_genero$CRC15_kraken)
Correlacion_CRC16_genero = cor.test(tabla_genero$CRC16,tabla_genero$CRC16_kraken)
Correlacion_CRC17_genero = cor.test(tabla_genero$CRC17,tabla_genero$CRC17_kraken)
Correlacion_CRC18_genero = cor.test(tabla_genero$CRC18,tabla_genero$CRC18_kraken)
Correlacion_CRC19_genero = cor.test(tabla_genero$CRC19,tabla_genero$CRC19_kraken)
Correlacion_CRC20_genero = cor.test(tabla_genero$CRC20,tabla_genero$CRC20_kraken)
Correlacion_CRC21_genero = cor.test(tabla_genero$CRC21,tabla_genero$CRC21_kraken)
Correlacion_CRC22_genero = cor.test(tabla_genero$CRC22,tabla_genero$CRC22_kraken)
Correlacion_CRC23_genero = cor.test(tabla_genero$CRC23,tabla_genero$CRC23_kraken)
Correlacion_CRC24_genero = cor.test(tabla_genero$CRC24,tabla_genero$CRC24_kraken)
Correlacion_CRC25_genero = cor.test(tabla_genero$CRC25,tabla_genero$CRC25_kraken)
Correlacion_CRC26_genero = cor.test(tabla_genero$CRC26,tabla_genero$CRC26_kraken)
Correlacion_CRC27_genero = cor.test(tabla_genero$CRC27,tabla_genero$CRC27_kraken)
Correlacion_CRC28_genero = cor.test(tabla_genero$CRC28,tabla_genero$CRC28_kraken)
Correlacion_CRC29_genero = cor.test(tabla_genero$CRC29,tabla_genero$CRC29_kraken)
Correlacion_CRC30_genero = cor.test(tabla_genero$CRC30,tabla_genero$CRC30_kraken)
Correlacion_CRC31_genero = cor.test(tabla_genero$CRC31,tabla_genero$CRC31_kraken)
Correlacion_CRC32_genero = cor.test(tabla_genero$CRC32,tabla_genero$CRC32_kraken)
Correlacion_CRC33_genero = cor.test(tabla_genero$CRC33,tabla_genero$CRC33_kraken)
Correlacion_CRC34_genero = cor.test(tabla_genero$CRC34,tabla_genero$CRC34_kraken)

Correlacion_CRC01_genero
Correlacion_CRC02_genero
Correlacion_CRC03_genero
Correlacion_CRC04_genero
Correlacion_CRC05_genero
Correlacion_CRC06_genero
Correlacion_CRC07_genero
Correlacion_CRC08_genero
Correlacion_CRC09_genero
Correlacion_CRC10_genero
Correlacion_CRC11_genero
Correlacion_CRC12_genero
Correlacion_CRC13_genero
Correlacion_CRC15_genero
Correlacion_CRC16_genero
Correlacion_CRC17_genero
Correlacion_CRC18_genero
Correlacion_CRC19_genero
Correlacion_CRC20_genero
Correlacion_CRC21_genero
Correlacion_CRC22_genero
Correlacion_CRC23_genero
Correlacion_CRC24_genero
Correlacion_CRC25_genero
Correlacion_CRC26_genero
Correlacion_CRC27_genero
Correlacion_CRC28_genero
Correlacion_CRC29_genero
Correlacion_CRC30_genero
Correlacion_CRC31_genero
Correlacion_CRC32_genero
Correlacion_CRC33_genero
Correlacion_CRC34_genero

# Grafica que muestra los valores de correlacion

ggplot(data=Resultados, aes(x=Sample, y=Correlation, fill=Significant)) +
  geom_bar(stat="identity")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



# 3) Creacion de las redes neuronales -------------------------------------

library(nnet)
library(NeuralNetTools)
library(ggplot2)

# A, B y C son con filo, clase y orden respectivamente.
# 1 y 2 Son con un datasheet y con dos, respectivamente.
# El que es solo un datasheet, se le ha quitado la muestra única con valor "Stage" = 4.
# Respecto a "Stage" son 6 tipo I, 13 tipo II, 13 tipo III y 19 tipo IV.

# A1) Filo un datasheet ---------------------------------------------------

#Cargar tabla

tabla_A1 <- read.table("Tabla_filo_primer_datasheet.txt", row.names = 1, header=TRUE)
tabla_A1 <- tabla_A1[!(row.names(tabla_A1) %in% "CRC32"), ] # Así eliminamos la única muestra con "Stage = 4".
tabla_A1$Stage = as.factor(tabla_A1$Stage)

# B) Configurar la parte para el entrenamiento

tabla_A1_sort <-tabla_A1[order(tabla_A1$Stage),]
set.seed(1000)
entren_tabla_A1 <-c(sample(1:6,2),sample(7:19,3),sample(20:32,3))

# C) Crear la red
set.seed(1000)
red_tabla_A1 = nnet(Stage ~ Bacteroidetes + Firmicutes +
                      Fusobacteria + Proteobacteria + Tenericutes +
                      Chloroflexi + Marinimicrobia.SAR406clade. +
                      Elusimicrobia + Cyanobacteria + Actinobacteria +
                      Patescibacteria + Epsilonbacteraeota + Verrucomicrobia +
                      Planctomycetes + Acidobacteria + Halanaerobiaeota +
                      Spirochaetes + Atribacteria + Latescibacteria +
                      LCP.89 + Gemmatimonadetes + Kiritimatiellaeota +
                      Omnitrophicaeota + PAUC34f + Deinococcus.Thermus +
                      Dependentiae + Fibrobacteres + Synergistetes +
                      Margulisbacteria + Nitrospinae + BRC1 +
                      Schekmanbacteria + Deferribacteres + Nitrospinae +
                      Aquificae + TA06 + Cloacimonetes +
                      Lentisphaerae + WPS.2 + Chlamydiae +
                      Armatimonadetes + Hydrothermae + FBP +
                      WS4 + Desantisbacteria + Zixibacteria +
                      Dadabacteria + Coprothermobacteraeota +
                      Calditrichaeota + CK.2C2.2 +
                      WS2 + Thermotogae + Aerophobetes +
                      Aegiribacteria + Caldiserica + Fervidibacteria +
                      GN01 + Hydrogenedentes + Rokubacteria +
                      Modulibacteria + Acetothermia + Entotheonellaeota +
                      WS1 + Poribacteria + WOR.1 + AncK6 +
                      FCPU426 + BHI80.139 + MAT.CR.M4.B07,data=tabla_A1_sort,subset = entren_tabla_A1, size=2, decay=1.0e-5, maxit=1000)

# D) Evaluar la red con una matriz de confusion

mc_tabla_A1 <-table(tabla_A1_sort$Stage[-entren_tabla_A1],predict(red_tabla_A1,tabla_A1_sort[-entren_tabla_A1,],type="class"))
predict(red_tabla_A1, Stage = tabla_A1_sort$Stage[-entren_tabla_A1], type="class")
actual_tabla_A1 <- tabla_A1_sort$Stage[-entren_tabla_A1]
preds_tabla_A1 <- predict(red_tabla_A1, tabla_A1_sort[-entren_tabla_A1, ], type="class")
mc_tabla_A1 <- table(actual_tabla_A1, preds_tabla_A1)
print(mc_tabla_A1)

aciertos = mc_tabla_A1[1,1] + mc_tabla_A1[2,2] + mc_tabla_A1[3,3]
precision = (aciertos / 24) * 100
message("Precision de la red neuronal del ", round(precision)," %")

olden(red_tabla_A1) + theme(axis.text.x = element_text(angle = 90)) + ggtitle("Importancia de cada variable: Filo, un datasheet")

neuralweights(red_tabla_A1)


# A2) Filo dos datasheets ---------------------------------------------------

#Cargar tabla

tabla_A2 <- read.table("Tabla_filo_primer_segundo_datasheets.txt", row.names = 1, header=TRUE)
tabla_A2$Stage = as.factor(tabla_A2$Stage)

# B) Configurar la parte para el entrenamiento

tabla_A2_sort <-tabla_A2[order(tabla_A2$Stage),]
set.seed(1000)
entren_tabla_A2 <-c(sample(1:6,2),sample(7:19,3),sample(20:32,3),sample(33:51,4))

# C) Crear la red
set.seed(1000)
red_tabla_A2 = nnet(Stage ~ Bacteroidetes + Firmicutes +
                      Fusobacteria + Proteobacteria + Tenericutes +
                      Chloroflexi + Marinimicrobia.SAR406clade. +
                      Elusimicrobia + Cyanobacteria + Actinobacteria +
                      Patescibacteria + Epsilonbacteraeota + Verrucomicrobia +
                      Planctomycetes + Acidobacteria + Halanaerobiaeota +
                      Spirochaetes + Atribacteria + Latescibacteria +
                      LCP.89 + Gemmatimonadetes + Kiritimatiellaeota +
                      Omnitrophicaeota + PAUC34f + Deinococcus.Thermus +
                      Dependentiae + Fibrobacteres + Synergistetes +
                      Margulisbacteria + Nitrospinae + BRC1 +
                      Schekmanbacteria + Deferribacteres + Nitrospinae +
                      Aquificae + TA06 + Cloacimonetes +
                      Lentisphaerae + WPS.2 + Chlamydiae +
                      Armatimonadetes + Hydrothermae + FBP +
                      WS4 + Desantisbacteria + Zixibacteria +
                      Dadabacteria + Coprothermobacteraeota +
                      Calditrichaeota + CK.2C2.2 +
                      WS2 + Thermotogae + Aerophobetes +
                      Aegiribacteria + Caldiserica + Fervidibacteria +
                      GN01 + Hydrogenedentes + Rokubacteria +
                      Modulibacteria + Acetothermia + Entotheonellaeota +
                      WS1 + Poribacteria + WOR.1 + AncK6 +
                      FCPU426 + BHI80.139 + MAT.CR.M4.B07,data=tabla_A2_sort,subset = entren_tabla_A2, size=2, decay=1.0e-5, maxit=1000)

# D) Evaluar la red con una matriz de confusion

mc_tabla_A2 <-table(tabla_A2_sort$Stage[-entren_tabla_A2],predict(red_tabla_A2,tabla_A2_sort[-entren_tabla_A2,],type="class"))
predict(red_tabla_A2, Stage = tabla_A2_sort$Stage[-entren_tabla_A2], type="class")
actual_tabla_A2 <- tabla_A2_sort$Stage[-entren_tabla_A2]
preds_tabla_A2 <- predict(red_tabla_A2, tabla_A2_sort[-entren_tabla_A2, ], type="class")
mc_tabla_A2 <- table(actual_tabla_A2, preds_tabla_A2)
print(mc_tabla_A2)

aciertos = mc_tabla_A2[1,1] + mc_tabla_A2[2,2] + mc_tabla_A2[3,3] + mc_tabla_A2[4,4]
precision = (aciertos / 39) * 100
message("Precision de la red neuronal del ", round(precision)," %")

olden(red_tabla_A2) + theme(axis.text.x = element_text(angle = 90)) + ggtitle("Importancia de cada variable: Filo, dos datasheets")


# B1) Clase un datasheet ---------------------------------------------------

#Cargar tabla

tabla_B1 <- read.table("tabla_clase_primer_datasheet.txt", row.names = 1, header=TRUE)
tabla_B1 <- tabla_B1[!(row.names(tabla_B1) %in% "CRC32"), ] # Así eliminamos la única muestra con "Stage = 4".
tabla_B1$Stage = as.factor(tabla_B1$Stage)

# B) Configurar la parte para el entrenamiento

tabla_B1_sort <-tabla_B1[order(tabla_B1$Stage),]
set.seed(1000)
entren_tabla_B1 <-c(sample(1:6,2),sample(7:19,3),sample(20:32,3))

# C) Crear la red
set.seed(1000)
red_tabla_B1 = nnet(Stage ~ Bacteroidia + Rhodothermia + Ignavibacteria + Chlorobia + Clostridia + Bacilli +
                      Negativicutes + Erysipelotrichia + Limnochordia + Fusobacteriia + Gammaproteobacteria +
                      Alphaproteobacteria + Deltaproteobacteria + Mollicutes + Dehalococcoidia + Anaerolineae +
                      Chloroflexia + Ktedonobacteria + AD3 + TK10 + SHA.26 + P2.11E + JG30.KF.CM66 +
                      Endomicrobia + Elusimicrobia + LineageIIb + LineageIIc + Oxyphotobacteria + Melainabacteria +
                      Sericytochromatia + Actinobacteria + Coriobacteriia + Acidimicrobiia + Thermoleophilia +
                      Saccharimonadia + Parcubacteria + Microgenomatia + ABY1 + Gracilibacteria + WWE3 + CPR2 +
                      WS6_Dojkabacteria + Campylobacteria + Verrucomicrobiae + Planctomycetacia + OM190 +
                      Phycisphaerae + Pla3lineage + vadinHA49 + Brocadiae + Holophagae + Acidobacteriia +
                      Subgroup6 + Subgroup18 + Blastocatellia_Subgroup4 + Aminicenantia + Subgroup22 +
                      Subgroup15 + Subgroup21 + Subgroup17 + Halanaerobiia + Spirochaetia + Brachyspirae + 
                      JS1 + Latescibacteria + BD2.11terrestrialgroup + Gemmatimonadetes + Longimicrobia + 
                      MD2902.B12 + Kiritimatiellae + Omnitrophia + Deinococci + Babeliae + Fibrobacteria + 
                      Synergistia + Thermodesulfovibrionia + Nitrospira + Deferribacteres + Nitrospinia + 
                      Aquificae + Cloacimonadia + Lentisphaeria + Chlamydiae + Chthonomonadetes + 
                      Fimbriimonadia + Dadabacteriia + Coprothermobacteria + Calditrichia + Thermotogae + 
                      KIST.JJY010 + Rubrobacteria + Berkelbacteria + Kazania + FFCH5909 + BD7.11 + Pla4lineage + 
                      uncultured + Caldisericia + NC10 + Nitriliruptoria + RBG.16.55.12 + KD4.96 + ODP123 + 
                      Leptospirae + Thermoanaerobaculia + Subgroup19 + Subgroup26 + PAUC43fmarinebenthicgroup + 
                      Armatimonadia + LD1.PA32 + Moduliflexia + Acetothermiia + Oligosphaeria + Entotheonellia + 
                      Zetaproteobacteria + X4.29.1 + WCHB1.81 + OC31 + Magnetococcia + MB.A2.108 + Gitt.GS.136 + 
                      Subgroup5 + d142 + Chitinivibrionia + BMS9AB35 + S0134terrestrialgroup + Desulfurobacteriia + 
                      ODP1230B23.02 + LineageIIa + SGST604 + Caldatribacteriia + Subgroup9 + BRH.c20a + N9D0 + 
                      Desulfurellia + X0319.7L14 + TK30 + MVP.15 + NLS2.31 + SPG12.343.353.B69 + Fischerbacteria + 
                      P9X2b3D02 + GBS.L1.B05 + Subgroup11 + DG.56 + Subgroup25 + AT.s3.28 + X028H05.P.BN.P5,data=tabla_B1_sort,subset = entren_tabla_B1, size=2, decay=1.0e-5, maxit=1000)

# D) Evaluar la red con una matriz de confusion

mc_tabla_B1 <-table(tabla_B1_sort$Stage[-entren_tabla_B1],predict(red_tabla_B1,tabla_B1_sort[-entren_tabla_B1,],type="class"))
predict(red_tabla_B1, Stage = tabla_B1_sort$Stage[-entren_tabla_B1], type="class")
actual_tabla_B1 <- tabla_B1_sort$Stage[-entren_tabla_B1]
preds_tabla_B1 <- predict(red_tabla_B1, tabla_B1_sort[-entren_tabla_B1, ], type="class")
mc_tabla_B1 <- table(actual_tabla_B1, preds_tabla_B1)
print(mc_tabla_B1)

aciertos = mc_tabla_B1[1,1] + mc_tabla_B1[2,2] + mc_tabla_B1[3,3]
precision = (aciertos / 24) * 100
message("Precision de la red neuronal del ", round(precision)," %")

olden(red_tabla_B1) + theme(axis.text.x = element_text(angle = 90)) + ggtitle("Importancia de cada variable: Clase, un datasheet")


# B2) Clase dos datasheets ---------------------------------------------------

#Cargar tabla

tabla_B2 <- read.table("tabla_clase_primer_segundo_datasheets.txt", row.names = 1, header=TRUE)
tabla_B2$Stage = as.factor(tabla_B2$Stage)

# B) Configurar la parte para el entrenamiento

tabla_B2_sort <-tabla_B2[order(tabla_B2$Stage),]
set.seed(1000)
entren_tabla_B2 <-c(sample(1:6,2),sample(7:19,3),sample(20:32,3),sample(33:51,4))

# C) Crear la red
set.seed(1000)
red_tabla_B2 = nnet(Stage ~ Bacteroidia + Rhodothermia + Ignavibacteria + Chlorobia + Clostridia + Bacilli +
                      Negativicutes + Erysipelotrichia + Limnochordia + Fusobacteriia + Gammaproteobacteria +
                      Alphaproteobacteria + Deltaproteobacteria + Mollicutes + Dehalococcoidia + Anaerolineae +
                      Chloroflexia + Ktedonobacteria + AD3 + TK10 + SHA.26 + P2.11E + JG30.KF.CM66 +
                      Endomicrobia + Elusimicrobia + LineageIIb + LineageIIc + Oxyphotobacteria + Melainabacteria +
                      Sericytochromatia + Actinobacteria + Coriobacteriia + Acidimicrobiia + Thermoleophilia +
                      Saccharimonadia + Parcubacteria + Microgenomatia + ABY1 + Gracilibacteria + WWE3 + CPR2 +
                      WS6_Dojkabacteria + Campylobacteria + Verrucomicrobiae + Planctomycetacia + OM190 +
                      Phycisphaerae + Pla3lineage + vadinHA49 + Brocadiae + Holophagae + Acidobacteriia +
                      Subgroup6 + Subgroup18 + Blastocatellia_Subgroup4 + Aminicenantia + Subgroup22 +
                      Subgroup15 + Subgroup21 + Subgroup17 + Halanaerobiia + Spirochaetia + Brachyspirae + 
                      JS1 + Latescibacteria + BD2.11terrestrialgroup + Gemmatimonadetes + Longimicrobia + 
                      MD2902.B12 + Kiritimatiellae + Omnitrophia + Deinococci + Babeliae + Fibrobacteria + 
                      Synergistia + Thermodesulfovibrionia + Nitrospira + Deferribacteres + Nitrospinia + 
                      Aquificae + Cloacimonadia + Lentisphaeria + Chlamydiae + Chthonomonadetes + 
                      Fimbriimonadia + Dadabacteriia + Coprothermobacteria + Calditrichia + Thermotogae + 
                      KIST.JJY010 + Rubrobacteria + Berkelbacteria + Kazania + FFCH5909 + BD7.11 + Pla4lineage + 
                      uncultured + Caldisericia + NC10 + Nitriliruptoria + RBG.16.55.12 + KD4.96 + ODP123 + 
                      Leptospirae + Thermoanaerobaculia + Subgroup19 + Subgroup26 + PAUC43fmarinebenthicgroup + 
                      Armatimonadia + LD1.PA32 + Moduliflexia + Acetothermiia + Oligosphaeria + Entotheonellia + 
                      Zetaproteobacteria + X4.29.1 + WCHB1.81 + OC31 + Magnetococcia + MB.A2.108 + Gitt.GS.136 + 
                      Subgroup5 + d142 + Chitinivibrionia + BMS9AB35 + S0134terrestrialgroup + Desulfurobacteriia + 
                      ODP1230B23.02 + LineageIIa + SGST604 + Caldatribacteriia + Subgroup9 + BRH.c20a + N9D0 + 
                      Desulfurellia + X0319.7L14 + TK30 + MVP.15 + NLS2.31 + SPG12.343.353.B69 + Fischerbacteria + 
                      P9X2b3D02 + GBS.L1.B05 + Subgroup11 + DG.56 + Subgroup25 + AT.s3.28 + X028H05.P.BN.P5,data=tabla_B2_sort,subset = entren_tabla_B2, size=2, decay=1.0e-5, maxit=1000)

# D) Evaluar la red con una matriz de confusion

mc_tabla_B2 <-table(tabla_B2_sort$Stage[-entren_tabla_B2],predict(red_tabla_B2,tabla_B2_sort[-entren_tabla_B2,],type="class"))
predict(red_tabla_B2, Stage = tabla_B2_sort$Stage[-entren_tabla_B2], type="class")
actual_tabla_B2 <- tabla_B2_sort$Stage[-entren_tabla_B2]
preds_tabla_B2 <- predict(red_tabla_B2, tabla_B2_sort[-entren_tabla_B2, ], type="class")
mc_tabla_B2 <- table(actual_tabla_B2, preds_tabla_B2)
print(mc_tabla_B2)

aciertos = mc_tabla_B2[1,1] + mc_tabla_B2[2,2] + mc_tabla_B2[3,3] + mc_tabla_B2[4,4]
precision = (aciertos / 39) * 100
message("Precision de la red neuronal del ", round(precision)," %")

olden(red_tabla_B2) + theme(axis.text.x = element_text(angle = 90)) + ggtitle("Importancia de cada variable: Clase, dos datasheets")



# C1) Orden un datasheet ---------------------------------------------------

#Cargar tabla

tabla_C1 <- read.table("tabla_orden_primer_datasheet.txt", row.names = 1, header=TRUE)
tabla_C1 <- tabla_C1[!(row.names(tabla_C1) %in% "CRC32"), ] # Así eliminamos la única muestra con "Stage = 4".
tabla_C1$Stage = as.factor(tabla_C1$Stage)

# B) Configurar la parte para el entrenamiento

tabla_C1_sort <-tabla_C1[order(tabla_C1$Stage),]
set.seed(1000)
entren_tabla_C1 <-c(sample(1:6,2),sample(7:19,3),sample(20:32,3))

# C) Crear la red
set.seed(1000)
red_tabla_C1 = nnet(Stage ~ Absconditabacteriales_SR1 + Acanthopleuribacterales + Acetobacterales + Acholeplasmatales +
                      Acidiferrobacterales + Acidimicrobiales + Acidithiobacillales + Acidobacteriales +
                      Actinomarinales + Actinomycetales + ADurbBin180 + Aeromonadales + AlphaproteobacteriaIncertaeSedis +
                      Alteromonadales + Aminicenantales + Anaerolineales + Anaeroplasmatales + Apal.E12 + Aquificales +
                      Arctic97B.4marinegroup + Ardenticatenales + Arenicellales + Armatimonadales +
                      AT.s2.59 + Azospirillales + B12.WMSP1 + B2M28 + Babeliales + Bacillales + Bacteroidales +
                      BacteroidetesVC2.1Bac22 + Balneolales + Bdellovibrionales + Beggiatoales + Betaproteobacteriales +
                      Bifidobacteriales + Blastocatellales + Brachyspirales + Bradymonadales + Brevinematales +
                      Brocadiales + C0119 + C86 + Caedibacterales + Caenarcaniphilales + Caldatribacteriales +
                      Caldilineales + Caldisericales + Calditrichales + Campylobacterales +
                      CandidatusAbawacabacteria + CandidatusAdlerbacteria + CandidatusAmesbacteria + 
                      CandidatusAzambacteria + CandidatusBrennerbacteria + CandidatusBuchananbacteria + 
                      CandidatusCampbellbacteria + CandidatusCollierbacteria + CandidatusColwellbacteria + 
                      CandidatusDaviesbacteria + CandidatusDoudnabacteria + CandidatusFalkowbacteria + 
                      CandidatusGottesmanbacteria + CandidatusJacksonbacteria + CandidatusJorgensenbacteria + 
                      CandidatusKaiserbacteria + CandidatusKerfeldbacteria + CandidatusKomeilibacteria + 
                      CandidatusKuenenbacteria + CandidatusLevybacteria + CandidatusMagasanikbacteria + 
                      CandidatusMoranbacteria + CandidatusNomurabacteria + CandidatusPacebacteria + 
                      CandidatusPeregrinibacteria + CandidatusPeribacteria + CandidatusPortnoybacteria +
                      CandidatusRoizmanbacteria + CandidatusShapirobacteria + 
                      CandidatusSpechtbacteria + CandidatusUhrbacteria + CandidatusVogelbacteria + CandidatusWoesebacteria +
                      CandidatusWolfebacteria + CandidatusWoykebacteria + CandidatusYanofskybacteria +
                      CandidatusZambryskibacteria + Cardiobacteriales + Catenulisporales + Caulobacterales +
                      CCD24 + CCM11a + CCM19a + Cellvibrionales + CG1.02.40.25 + CG1.02.42.13 + CHAB.XI.27 +
                      ChitinivibrioniaIncertaeSedis + Chitinophagales + Chlamydiales + Chloracidobacteriales +
                      Chlorobiales + Chloroflexales + Chloroplast + Chromatiales + Chthoniobacterales + Chthonomonadales +
                      Cloacimonadales + ClostridiaIncertaeSedis + Clostridiales + Competibacterales + Coprothermobacterales +
                      Coriobacteriales + Corynebacteriales + Coxiellales + Cytophagales + D8A.2 + Dadabacteriales +
                      Deferribacterales + Deinococcales + DeltaproteobacteriaIncertaeSedis + Desulfarculales +
                      Desulfobacterales + Desulfovibrionales + Desulfurellales + Desulfurobacteriales +
                      Desulfuromonadales + DG.20 + Diplorickettsiales + DMI + Dongiales + DS.100 + DscP2 + DTB120 +
                      DTU014 + EC3 + Ectothiorhodospirales + Elsterales + EMP.G18 + Endomicrobiales + Enterobacteriales +
                      Entomoplasmatales + Entotheonellales + EPR3968.O8a.Bc78 + Erysipelotrichales + EUB33.2 + eub62A3 +
                      Eurycoccales + Euzebyales + EV818SWSAP88 + FCPU453 + Fibrobacterales + Fimbriimonadales +
                      Flavobacteriales + Francisellales + Frankiales + FS117.23B.02 + FS118.23B.02 + Fusobacteriales +
                      FW113 + FW22 + Ga0077536 + Gaiellales + GammaproteobacteriaIncertaeSedis + Gastranaerophilales +
                      Gemmatales + Gemmatimonadales + GIF3 + GIF9 + Gloeobacterales + Glycomycetales + GWA2.38.13b +
                      GWB1.42.6 + H3.93 + Halanaerobiales + Haloplasmatales + Halothiobacillales + HOC36 + Holophagales +
                      Holosporales + Hydrogenedentiales + Hydrogenothermales + Ignavibacteriales + IMCC26256 +
                      Immundisolibacterales + Isosphaerales + Izimaplasmatales + JB111 + JG36.TzT.191 + JGI0000069.P22 +
                      JTB23 + Kallotenuales + KI89Aclade + Kineosporiales + Kiritimatiellales + Kordiimonadales +
                      Kosmotogales + Kryptoniales + Ktedonobacterales + Lactobacillales + Latescibacterales + LD1.PB3 +
                      Legionellales + Lentisphaerales + Leptolyngbyales + Leptospirales + Limnochordales + Limnotrichales +
                      LineageIV + Longimicrobiales + Magnetococcales + Mariprofundales + MBA03 + MBAE14 + MBMPE27 +
                      MBNT15 + Mesoaciditogales + Methylacidiphilales + Methylococcales + Methylomirabilales +
                      Micavibrionales + Micrococcales + Micromonosporales + Micropepsales + Microtrichales +
                      Milano.WF1B.44 + ML602M.17 + ML635J.38 + mle1.8 + Moduliflexales + MollicutesIncertaeSedis +
                      MollicutesRF39 + MSB.5B2 + MSB.5E12 + MSBL2 + MSBL5 + MSBL9 + MVP.21 + MVP.88 + Mycoplasmatales + 
                      Myxococcales + Napoli.4B.65 + Natranaerobiales + Nautiliales + NB1.j + Nitriliruptorales + 
                      Nitrococcales + Nitrosococcales + Nitrospinales + Nitrospirales + NKB15 + Nostocales + 
                      Obscuribacterales + Oceanospirillales + Oligoflexales + Oligosphaerales + OM182clade + 
                      Omnitrophales + OPB41 + OPB56 + Opitutales + Orbales + OxyphotobacteriaIncertaeSedis + P.palmC41 + 
                      Paracaedibacterales + Parvibaculales + Pasteurellales + PB19 + Pedosphaerales + PeM15 + 
                      Petrotogales + Phormidesmiales + Phycisphaerales + Pirellulales + Piscirickettsiales + 
                      Pla1lineage + Planctomycetales + PLTA13 + Propionibacteriales + Pseudanabaenales + 
                      Pseudomonadales + Pseudonocardiales + Puniceispirillales + Pyrinomonadales + R7C24 +
                      RBG.13.46.9 + RBG.13.54.9 + RCP2.54 + RD011 + Reyranellales + Rhizobiales +
                      Rhodobacterales + Rhodospirillales + Rhodothermales + Rhodovibrionales + Rickettsiales +
                      Rokubacteriales + Rubrobacterales + Run.SP154 + S085 + S.70 + Saccharimonadales +
                      Salinisphaerales + SAR11clade + SAR202clade + SAR324clade_MarinegroupB + SAR86clade +
                      S.BQ2.57soilgroup + SBR1031 + sediment.surface35 + Selenomonadales + Sh765B.TzT.20 + SJA.15 +
                      SJA.28 + SM1A07 + SM23.32 + Sneathiellales + Solibacterales + Solirubrobacterales +
                      Sphingobacteriales + Sphingomonadales + Spirochaetales + SS1.B.02.17 + SS1.B.04.55 +
                      Steroidobacterales + Streptomycetales + Streptosporangiales + Subgroup13 + Subgroup2 +
                      Subgroup7 + Sva0485 + Synechococcales + Synergistales + Syntrophobacterales + Tenderiales +
                      Tepidisphaerales + Thalassobaculales + Thermales + Thermoanaerobacterales + Thermoanaerobaculales +
                      Thermobaculales + Thermodesulfobacteriales + Thermodesulfovibrionales + Thermoflexales +
                      Thermolithobacterales + Thermomicrobiales + Thermosynechococcales + Thermotogales +
                      Thiohalorhabdales + Thiomicrospirales + Thiotrichales + Tistrellales + UBA10353marinegroup +
                      UBA9983 + uncultured + UnknownOrder + vadinBA26 + Vampirovibrionales + Verrucomicrobiales +
                      Vibrionales + Victivallales + WCHB1.41 + WN.HWB.116 + Xanthomonadales,data=tabla_C1_sort,subset = entren_tabla_C1, size=2, decay=1.0e-5, maxit=1000)


# D) Evaluar la red con una matriz de confusion

mc_tabla_C1 <-table(tabla_C1_sort$Stage[-entren_tabla_C1],predict(red_tabla_C1,tabla_C1_sort[-entren_tabla_C1,],type="class"))
predict(red_tabla_C1, Stage = tabla_C1_sort$Stage[-entren_tabla_C1], type="class")
actual_tabla_C1 <- tabla_C1_sort$Stage[-entren_tabla_C1]
preds_tabla_C1 <- predict(red_tabla_C1, tabla_C1_sort[-entren_tabla_C1, ], type="class")
mc_tabla_C1 <- table(actual_tabla_C1, preds_tabla_C1)
print(mc_tabla_C1)

aciertos = mc_tabla_C1[1,1] + mc_tabla_C1[2,2] + mc_tabla_C1[3,3]
precision = (aciertos / 24) * 100
message("Precision de la red neuronal del ", round(precision)," %")

olden(red_tabla_C1) + theme(axis.text.x = element_text(angle = 90)) + ggtitle("Importancia de cada variable: Orden, un datasheet")


# C2) Orden dos datasheets ---------------------------------------------------

#Cargar tabla

tabla_C2 <- read.table("tabla_orden_primer_segundo_datasheets.txt", row.names = 1, header=TRUE)
tabla_C2$Stage = as.factor(tabla_C2$Stage)

# B) Configurar la parte para el entrenamiento

tabla_C2_sort <-tabla_C2[order(tabla_C2$Stage),]
set.seed(1000)
entren_tabla_C2 <-c(sample(1:6,2),sample(7:19,3),sample(20:32,3),sample(33:51,4))

# C) Crear la red
set.seed(1000)
red_tabla_C2 = nnet(Stage ~ Absconditabacteriales_SR1 + Acanthopleuribacterales + Acetobacterales + Acholeplasmatales +
                      Acidiferrobacterales + Acidimicrobiales + Acidithiobacillales + Acidobacteriales +
                      Actinomarinales + Actinomycetales + ADurbBin180 + Aeromonadales + AlphaproteobacteriaIncertaeSedis +
                      Alteromonadales + Aminicenantales + Anaerolineales + Anaeroplasmatales + Apal.E12 + Aquificales +
                      Arctic97B.4marinegroup + Ardenticatenales + Arenicellales + Armatimonadales +
                      AT.s2.59 + Azospirillales + B12.WMSP1 + B2M28 + Babeliales + Bacillales + Bacteroidales +
                      BacteroidetesVC2.1Bac22 + Balneolales + Bdellovibrionales + Beggiatoales + Betaproteobacteriales +
                      Bifidobacteriales + Blastocatellales + Brachyspirales + Bradymonadales + Brevinematales +
                      Brocadiales + C0119 + C86 + Caedibacterales + Caenarcaniphilales + Caldatribacteriales +
                      Caldilineales + Caldisericales + Calditrichales + Campylobacterales +
                      CandidatusAbawacabacteria + CandidatusAdlerbacteria + CandidatusAmesbacteria + 
                      CandidatusAzambacteria + CandidatusBrennerbacteria + CandidatusBuchananbacteria + 
                      CandidatusCampbellbacteria + CandidatusCollierbacteria + CandidatusColwellbacteria + 
                      CandidatusDaviesbacteria + CandidatusDoudnabacteria + CandidatusFalkowbacteria + 
                      CandidatusGottesmanbacteria + CandidatusJacksonbacteria + CandidatusJorgensenbacteria + 
                      CandidatusKaiserbacteria + CandidatusKerfeldbacteria + CandidatusKomeilibacteria + 
                      CandidatusKuenenbacteria + CandidatusLevybacteria + CandidatusMagasanikbacteria + 
                      CandidatusMoranbacteria + CandidatusNomurabacteria + CandidatusPacebacteria + 
                      CandidatusPeregrinibacteria + CandidatusPeribacteria + CandidatusPortnoybacteria +
                      CandidatusRoizmanbacteria + CandidatusShapirobacteria + 
                      CandidatusSpechtbacteria + CandidatusUhrbacteria + CandidatusVogelbacteria + CandidatusWoesebacteria +
                      CandidatusWolfebacteria + CandidatusWoykebacteria + CandidatusYanofskybacteria +
                      CandidatusZambryskibacteria + Cardiobacteriales + Catenulisporales + Caulobacterales +
                      CCD24 + CCM11a + CCM19a + Cellvibrionales + CG1.02.40.25 + CG1.02.42.13 + CHAB.XI.27 +
                      ChitinivibrioniaIncertaeSedis + Chitinophagales + Chlamydiales + Chloracidobacteriales +
                      Chlorobiales + Chloroflexales + Chloroplast + Chromatiales + Chthoniobacterales + Chthonomonadales +
                      Cloacimonadales + ClostridiaIncertaeSedis + Clostridiales + Competibacterales + Coprothermobacterales +
                      Coriobacteriales + Corynebacteriales + Coxiellales + Cytophagales + D8A.2 + Dadabacteriales +
                      Deferribacterales + Deinococcales + DeltaproteobacteriaIncertaeSedis + Desulfarculales +
                      Desulfobacterales + Desulfovibrionales + Desulfurellales + Desulfurobacteriales +
                      Desulfuromonadales + DG.20 + Diplorickettsiales + DMI + Dongiales + DS.100 + DscP2 + DTB120 +
                      DTU014 + EC3 + Ectothiorhodospirales + Elsterales + EMP.G18 + Endomicrobiales + Enterobacteriales +
                      Entomoplasmatales + Entotheonellales + EPR3968.O8a.Bc78 + Erysipelotrichales + EUB33.2 + eub62A3 +
                      Eurycoccales + Euzebyales + EV818SWSAP88 + FCPU453 + Fibrobacterales + Fimbriimonadales +
                      Flavobacteriales + Francisellales + Frankiales + FS117.23B.02 + FS118.23B.02 + Fusobacteriales +
                      FW113 + FW22 + Ga0077536 + Gaiellales + GammaproteobacteriaIncertaeSedis + Gastranaerophilales +
                      Gemmatales + Gemmatimonadales + GIF3 + GIF9 + Gloeobacterales + Glycomycetales + GWA2.38.13b +
                      GWB1.42.6 + H3.93 + Halanaerobiales + Haloplasmatales + Halothiobacillales + HOC36 + Holophagales +
                      Holosporales + Hydrogenedentiales + Hydrogenothermales + Ignavibacteriales + IMCC26256 +
                      Immundisolibacterales + Isosphaerales + Izimaplasmatales + JB111 + JG36.TzT.191 + JGI0000069.P22 +
                      JTB23 + Kallotenuales + KI89Aclade + Kineosporiales + Kiritimatiellales + Kordiimonadales +
                      Kosmotogales + Kryptoniales + Ktedonobacterales + Lactobacillales + Latescibacterales + LD1.PB3 +
                      Legionellales + Lentisphaerales + Leptolyngbyales + Leptospirales + Limnochordales + Limnotrichales +
                      LineageIV + Longimicrobiales + Magnetococcales + Mariprofundales + MBA03 + MBAE14 + MBMPE27 +
                      MBNT15 + Mesoaciditogales + Methylacidiphilales + Methylococcales + Methylomirabilales +
                      Micavibrionales + Micrococcales + Micromonosporales + Micropepsales + Microtrichales +
                      Milano.WF1B.44 + ML602M.17 + ML635J.38 + mle1.8 + Moduliflexales + MollicutesIncertaeSedis +
                      MollicutesRF39 + MSB.5B2 + MSB.5E12 + MSBL2 + MSBL5 + MSBL9 + MVP.21 + MVP.88 + Mycoplasmatales + 
                      Myxococcales + Napoli.4B.65 + Natranaerobiales + Nautiliales + NB1.j + Nitriliruptorales + 
                      Nitrococcales + Nitrosococcales + Nitrospinales + Nitrospirales + NKB15 + Nostocales + 
                      Obscuribacterales + Oceanospirillales + Oligoflexales + Oligosphaerales + OM182clade + 
                      Omnitrophales + OPB41 + OPB56 + Opitutales + Orbales + OxyphotobacteriaIncertaeSedis + P.palmC41 + 
                      Paracaedibacterales + Parvibaculales + Pasteurellales + PB19 + Pedosphaerales + PeM15 + 
                      Petrotogales + Phormidesmiales + Phycisphaerales + Pirellulales + Piscirickettsiales + 
                      Pla1lineage + Planctomycetales + PLTA13 + Propionibacteriales + Pseudanabaenales + 
                      Pseudomonadales + Pseudonocardiales + Puniceispirillales + Pyrinomonadales + R7C24 +
                      RBG.13.46.9 + RBG.13.54.9 + RCP2.54 + RD011 + Reyranellales + Rhizobiales +
                      Rhodobacterales + Rhodospirillales + Rhodothermales + Rhodovibrionales + Rickettsiales +
                      Rokubacteriales + Rubrobacterales + Run.SP154 + S085 + S.70 + Saccharimonadales +
                      Salinisphaerales + SAR11clade + SAR202clade + SAR324clade_MarinegroupB + SAR86clade +
                      S.BQ2.57soilgroup + SBR1031 + sediment.surface35 + Selenomonadales + Sh765B.TzT.20 + SJA.15 +
                      SJA.28 + SM1A07 + SM23.32 + Sneathiellales + Solibacterales + Solirubrobacterales +
                      Sphingobacteriales + Sphingomonadales + Spirochaetales + SS1.B.02.17 + SS1.B.04.55 +
                      Steroidobacterales + Streptomycetales + Streptosporangiales + Subgroup13 + Subgroup2 +
                      Subgroup7 + Sva0485 + Synechococcales + Synergistales + Syntrophobacterales + Tenderiales +
                      Tepidisphaerales + Thalassobaculales + Thermales + Thermoanaerobacterales + Thermoanaerobaculales +
                      Thermobaculales + Thermodesulfobacteriales + Thermodesulfovibrionales + Thermoflexales +
                      Thermolithobacterales + Thermomicrobiales + Thermosynechococcales + Thermotogales +
                      Thiohalorhabdales + Thiomicrospirales + Thiotrichales + Tistrellales + UBA10353marinegroup +
                      UBA9983 + uncultured + UnknownOrder + vadinBA26 + Vampirovibrionales + Verrucomicrobiales +
                      Vibrionales + Victivallales + WCHB1.41 + WN.HWB.116 + Xanthomonadales,data=tabla_C2_sort,subset = entren_tabla_C2, size=2, decay=1.0e-5, maxit=1000)


# D) Evaluar la red con una matriz de confusion

mc_tabla_C2 <-table(tabla_C2_sort$Stage[-entren_tabla_C2],predict(red_tabla_C2,tabla_C2_sort[-entren_tabla_C2,],type="class"))
predict(red_tabla_C2, Stage = tabla_C2_sort$Stage[-entren_tabla_C2], type="class")
actual_tabla_C2 <- tabla_C2_sort$Stage[-entren_tabla_C2]
preds_tabla_C2 <- predict(red_tabla_C2, tabla_C2_sort[-entren_tabla_C2, ], type="class")
mc_tabla_C2 <- table(actual_tabla_C2, preds_tabla_C2)
print(mc_tabla_C2)

aciertos = mc_tabla_C2[1,1] + mc_tabla_C2[2,2] + mc_tabla_C2[3,3] + mc_tabla_C2[4,4]
precision = (aciertos / 39) * 100
message("Precision de la red neuronal del ", round(precision)," %")

olden(red_tabla_C2) + theme(axis.text.x = element_text(angle = 90)) + ggtitle("Importancia de cada variable: Orden, dos datasheets")

archivo = olden(red_tabla_C2)





