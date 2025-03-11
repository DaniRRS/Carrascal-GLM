#####################################
##  Modelos de efectos aleatorios  ##
##    MODELOS MIXTOS GENERALES     ##
##   lme4::lmer   y   nlme::lme    ##
#####################################

## Luis M. Carrascal
## https://lmcarrascal.eu
## 24/11/2024
## corrido bajo R-4.3.0
##    https://cran.r-project.org/bin/windows/base/old/
## en RStudio 2023.12.1 (Build 402)
##    https://dailies.rstudio.com/rstudio/ocean-storm/electron/windows/2023-12-1-402/


#### REFERENCIAS EN LINEA ####
## https://en.wikipedia.org/wiki/Multilevel_model
## https://www.sciencedirect.com/science/article/pii/S0169534709000196
## https://peerj.com/articles/4794/
## https://peerj.com/articles/9522/
## http://www.john-ros.com/Rcourse/lme.html?utm_source=pocket_mylist
## https://m-clark.github.io/mixed-models-with-R/
## http://www.bristol.ac.uk/cmm/learning/videos/random-slopes.html?utm_source=pocket_mylist
## http://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf (paginas 6 y 30)
## https://rpsychologist.com/r-guide-longitudinal-lme-lmer
## https://fromthebottomoftheheap.net/2021/02/02/random-effects-in-gams/
## https://www.crumplab.com/psyc7709_2019/book/docs/a-tutorial-for-using-the-lme-function-from-the-nlme-package-.html
## https://stat.ethz.ch/pipermail/r-help/2006-October/115572.html
## https://peerj.com/articles/12794/
## https://cran.r-project.org/web/packages/glmmTMB/vignettes/glmmTMB.pdf
## https://backend.orbit.dtu.dk/ws/files/154739064/Publishers_version.pdf  (mirad los ejemplos del Appendix A)
## https://cran.r-project.org/web/packages/glmmTMB/vignettes/model_evaluation.pdf
## https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html
## https://besjournals.onlinelibrary.wiley.com/doi/pdfdirect/10.1111/2041-210X.13434
## https://m-clark.github.io/posts/2019-05-14-shrinkage-in-mixed-models/
## https://stats.stackexchange.com/questions/37647/what-is-the-minimum-recommended-number-of-groups-for-a-random-effects-factor?utm_source=pocket_mylist
## https://errickson.net/stats-notes/vizrandomeffects.html




#### IMPORTACION DE DATOS ####
## para importar los datos fuente
## datos <- read.table("clipboard", header=T, sep="\t", dec=".")    ## importar desde portapapeles
## importar en MAC desde el portapapeles: datos <- read.table(pipe("pbpaste"), header=T, sep="\t", dec=".")
## para exportar 
## write.table(datos, "c:/datos/nombre.txt", sep="\t", row.names=FALSE)                          ## lo guarda en el directorio especificado

## datos de ejemplo obtenidos de una URL
meconecto <- url("http://www.lmcarrascal.eu/cursos/nwayREPEATED.RData")  ## abro una conexion con internet
load(meconecto)                                                          ## cargo el archivo de la conexion
close(meconecto); rm(meconecto)                                          ## cierro la conexion
## convertimos el identificador numerico de los individuos (id; factor aleatorio) en un factor
datos$id <- as.factor(datos$id)


#### CARGAMOS PAQUETES DE ANALISIS ####
library(lme4)        ## modelos generales lmer
library(nlme)        ## modelos generales lineales lme
library(glmmTMB)     ## generalized mixed models (glmm) con TMB
library(DHARMa)      ## para diagnosis de modelos
library(moments)     ## para obtener el sesgo, kurtosis y sus significaciones: kurtosis, anscombe.test, skewness, agostino.test
library(lmerTest)    ## para MS, df, p ... usando type 3/type 1 hypotheses with "Satterthwaite" and "Kenward-Roger"
library(pbkrtest)    ## necesario para lmerTest
library(LMERConvenienceFunctions)
library(merTools)    ## para residuos de DHARMa
library(car)         ## para Anova(modelo, type=3) para equivalente a suma de tipo III; para boxCox(modelo, lambda=seq(-2,2, 1/100))
library(MuMIn)       ## para AICc
library(lmtest)      ## para lrtest
library(lattice)     ## para plots de residuos
library(psych)       ## para estadisticas resumen de tablas
library(performance) ## para pruebas de los modelos
library(sjPlot)      ## para efectos parciales
library(corrplot)    ## para graficos de correlaciones
library(parameters)  ## para salidas de tablas de coeficientes
library(effectsize)  ## para effectos parciales
library(effects)     ## para visualizacion de efectos parciales


#### Definimos los TIPOS DE CONTRASTES del modelo ####
##  https://rstudio-pubs-static.s3.amazonaws.com/84177_4604ecc1bae246c9926865db53b6cc29.html
##  https://rpubs.com/monajhzhu/608609
##  https://www.clayford.net/statistics/tag/sum-contrasts/
##
options(contrasts=c(factor="contr.sum", ordered="contr.poly"))
##
## otra posibilidad es incluir lo siguiente en la funcion:
##   contrasts=list(factor.mio1=contr.sum, factor.mio2=contr.sum, factor.mio3=contr.sum, ...)
## Por defecto R opera con contr.treatment que no es correcto: 
##   "[It] contrasts each level with the baseline level (specified by base); the baseline level is omitted. 
##    Note that this does not produce 'contrasts' as defined in the standard theory for linear models
##    as they are not orthogonal to the intercept"
## contr.treatment: obtiene para los niveles de los factores los EFECTOS MARGINALES respecto al intercepto, que es el valor medio en el nivel de referencia
##                  los coeficientes de los niveles de los factores miden los cambios respecto al valor del nivel de referencia
## contr.sum: obtiene para los niveles de los factores los EFECTOS PRINCIPALES respecto a la gran media (promedio de todos los niveles de los factores)
##            los coeficientes de los niveles de los factores miden los cambios respecto a la gran media.
## contr.sum asegura que todos los contrastes sumen cero, de manera que el "intercepto" sea la gran media. 
##   Los efectos se resumen con coeficientes que representan el numero de niveles de factor (k) menos 1, eliminando el ultimo nivel. 
##   Los coeficientes en contr.sum representan lna diferencia media respecto a la "gran media" para los primeros niveles del factor (del 1 al k-1).
## ved la diferencia para los contrastes en un factor con 5 niveles
##   el baseline es el ultimo nivel, y las columnas suman CERO (ortogonalidad)
contr.sum(5)
mat.desv1 <- matrix(c(1,0,0,0,-1, 0,1,0,0,-1, 0,0,1,0,-1, 0,0,0,1,-1), ncol=4)               ## igual que el anterior
mat.desv2 <- matrix(c(4,-1,-1,-1,-1, -1,4,-1,-1,-1, -1,-1,4,-1,-1, -1,-1,-1,4,-1), ncol=4)   ## equivalente al anterior pero con diferentes coeficiente
##   el baseline es el primer nivel, y las columnas NO suman CERO (ausencia de ortogonalidad)
contr.treatment(5) 
## para ver la secuencia ordenada de los niveles de un factor
levels(datos$especie)
levels(datos$fuente)
## para cambiar la ordenacion de los niveles
datos$fuente2 <- factor(datos$fuente, levels=c("Luz", "Infra"))
levels(datos$fuente2)



#### gestion de problemas de estimas y convergencia ####
##   introducid en lmer(..., control=control.lmer)      lme(..., control=control.lme)
control.lmer <- lmerControl(check.conv.grad=.makeCC(action ="ignore", tol=1e-6, relTol=NULL), optimizer="bobyqa", optCtrl=list(maxfun=100000))
control.lme <- lmeControl(maxIter=200, msMaxIter=200, msVerbose=TRUE, opt="optim")    ## el optimizador  opt="nlminb"  hay veces que da problemas
## control.lmer <- lmerControl(check.conv.grad=.makeCC(action ="ignore", tol=1e-6, relTol=NULL), optCtrl=list(maxfun=100000))
## control.lmer <- lmerControl(check.conv.grad=.makeCC(action ="ignore", tol=1e-6, relTol=NULL), optimizer="Nelder_Mead", optCtrl=list(maxfun=100000))



#### COMENTARIOS SOBRE MODELOS lmer Y lme ####
## Los modelos nlme::lme no gestionan factores aleatorios perfectamente cruzados; solo efectos anidados
##   si en una estructura aleatoria se incluyen dos efectos aleatorios el segundo se anida dentro del primero
##     random=list(especie=~1, id=~1)  implica especie/id   (estructura (1|especie/id) en lme4::lmer)
##   si una estructura random=list(...) incluye el mismo factor aleatorio, esos efectos no se anidan
##
## El metodo ML produce estimas (negativamente) sesgadas de las componentes de la varianza.
##     ML produce menores valores MSE (mean-squared errors) que REML
## El metodo REML es preferido para bajos tamannios muestrales donde el problema de ML se ve incrementado
##     REML se usa para comparar modelos que comparten exactamente la misma parte fija 
##     REML estima mejores valores medios de los coeficientes, pero con mayores variabilidades en sus estimas
##
## En modelos lmer(...) frecuentemente aparece esta alarma: boundary (singular) fit
##   en este caso se obtiene un ajuste "singular", indicativo de que el modelo esta sobreajustado 
##     o que los efectos aleatorios son muy pequennios (<1.0e-3)
##     (i.e., la estructura de efectos aleatorios es demasiado compleja para ser apoyada por los datos)
##   esto sugiere eliminar la parte mas compleja de la estructura de los efectos aleatorios (generalmente "random slopes"),
##     que conduce a un modelo mas parsimonioso, no sobreajustado.



#### DIFERENTES MODELOS ####
##  desarrollados con: lme4::lmer  y nlme::lme

## en el ejemplo, el factor fijo especie no puede cruzarse con el factor aleatorio individuos (id)

## *** PONED AQUI VUESTROS DATOS *** (las variables respuesta y predictoras a utilizar)

## Intraclass correlation: the "unconditional means model" (or "null model"); RANDOM INTERCEPTS WITH FIXED MEAN
##   diferencias entre los niveles del factor aleatorio (individuo)
##   lmer.0 aproximadamente como todos los modelos lme(...)
##   los cuatro modelos lme(...) producen los mismos resultados
lmer.0 <- lmer(a_por_b ~ 1 + (1|id), data=datos, control=control.lmer, REML=TRUE)
##
lme.0 <- lme(a_por_b ~ 1, random=~1|id, data=datos, control=control.lme, method="REML")     ## identico a las dos lineas siguientes 
lme.0 <- lme(a_por_b ~ 1, random=list(id=~1), data=datos, control=control.lme, method="REML")
lme.0 <- lme(a_por_b ~ 1, random=list(id=pdLogChol(~1)), data=datos, control=control.lme, method="REML")
lme.0a <- lme(a_por_b ~ 1, random=list(id=pdSymm(~1)), data=datos, control=control.lme, method="REML")
lme.0b <- lme(a_por_b ~ 1, random=list(id=pdDiag(~1)), data=datos, control=control.lme, method="REML")
lme.0c <- lme(a_por_b ~ 1, random=list(id=pdIdent(~1)), data=datos, control=control.lme, method="REML")

## modelo de RANDOM INTERCEPTS & FIXED SLOPES; lmer.1 es aproximadamente como lme.1
##   lmer.1 aproximadamente como todos los modelos lme(...)
##   los cuatro modelos lme(...) producen los mismos resultados
lmer.1 <- lmer(a_por_b ~ t_sonda+fuente*especie + (1|id), data=datos, control=control.lmer, REML=TRUE)
##
lme.1 <- lme(a_por_b ~ t_sonda+fuente*especie, random=list(id=~1), data=datos, control=control.lme, method="REML")
lme.1a <- lme(a_por_b ~ t_sonda+fuente*especie, random=list(id=pdSymm(~1)), data=datos, control=control.lme, method="REML")
lme.1b <- lme(a_por_b ~ t_sonda+fuente*especie, random=list(id=pdDiag(~1)), data=datos, control=control.lme, method="REML")
lme.1c <- lme(a_por_b ~ t_sonda+fuente*especie, random=list(id=pdIdent(~1)), data=datos, control=control.lme, method="REML")

## modelo de CORRELATED RANDOM INTERCEPTS & RANDOM SLOPES FOR TWO CORRELATED FIXED EFFECTS
##   lmer.2 "Kenward-Roger" aproximadamente como el modelo lme(...) lme.2b pdCompSymm
##   lmer.2 "Satterthwaite" aproximadamente como el modelo lme(...) lme.2a pdSymm
lmer.2 <- lmer(a_por_b ~ t_sonda+fuente*especie + (1+t_sonda+fuente|id), data=datos, control=control.lmer, REML=TRUE)
##
lme.2 <- lme(a_por_b ~ t_sonda+fuente*especie, random=~1+t_sonda+fuente|id, data=datos, control=control.lme, method="REML")        ## identico a la linea siguiente
lme.2 <- lme(a_por_b ~ t_sonda+fuente*especie, random=list(id=~1+t_sonda+fuente), data=datos, control=control.lme, method="REML")  ## identico a la linea siguiente
lme.2 <- lme(a_por_b ~ t_sonda+fuente*especie, random=list(id=pdLogChol(~1+t_sonda+fuente)), data=datos, control=control.lme, method="REML")
lme.2a <- lme(a_por_b ~ t_sonda+fuente*especie, random=list(id=pdSymm(~1+t_sonda+fuente)), data=datos, control=control.lme, method="REML")
lme.2b <- lme(a_por_b ~ t_sonda+fuente*especie, random=list(id=pdCompSymm(~1+t_sonda+fuente)), data=datos, control=control.lme, method="REML")

## modelo de CORRELATED RANDOM INTERCEPTs & RANDOM SLOPES FOR TWO UNCORRELATED FIXED EFFECTS
##   el orden de los terminos aleatorios dentro de lme(..., random=list(...)) no altera los resultados
##   lmer.3 "Satterthwaite" es aproximadamente como lme.3a pdSymm
lmer.3 <- lmer(a_por_b ~ t_sonda+fuente*especie + (1+t_sonda|id)+(1+fuente|id), data=datos, control=control.lmer, REML=TRUE)
##
lme.3 <- lme(a_por_b ~ t_sonda+fuente*especie, random=list(id=~1+t_sonda, id=~1+fuente), data=datos, control=control.lme, method="REML")
lme.3a <- lme(a_por_b ~ t_sonda+fuente*especie, random=list(id=pdSymm(~1+t_sonda), id=pdSymm(~1+fuente)), data=datos, control=control.lme, method="REML")
lme.3b <- lme(a_por_b ~ t_sonda+fuente*especie, random=list(id=pdCompSymm(~1+t_sonda), id=pdCompSymm(~1+fuente)), data=datos, control=control.lme, method="REML")

## modelo de UNCORRELATED RANDOM INTERCEPT & RANDOM SLOPES FOR TWO UNCORRELATED FIXED EFFECTS
##   elimina la covarianza intercepto-pendientes en cada efecto fijo 
##   el orden de los terminos aleatorios dentro de lme(..., random=list(...)) no altera los resultados
##   lmer.4 "Satterthwaite" es aproximadamente como los cuatro modelos lme(...) 
##   los modelos lme pdDiag y pdIdent asumen ausencia de correlacion entre interceptos y pendientes para los efectos fijos
##      ademas, estiman menos parametros al asumir que las covarianzas son 0, y pdIdent asumir que las varianzas son identicas
lmer.4 <- lmer(a_por_b ~ t_sonda+fuente*especie + (1|id)+(t_sonda-1|id)+(fuente-1|id), data=datos, control=control.lmer, REML=TRUE)
##
lme.4 <- lme(a_por_b ~ t_sonda+fuente*especie, random=list(id=~1, id=~t_sonda-1, id=~fuente-1), data=datos, control=control.lme, method="REML")
lme.4a <- lme(a_por_b ~ t_sonda+fuente*especie, random=list(id=pdSymm(~1), id=pdSymm(~t_sonda-1), id=pdSymm(~fuente-1)), data=datos, control=control.lme, method="REML")
lme.4b <- lme(a_por_b ~ t_sonda+fuente*especie, random=list(id=pdDiag(~1+t_sonda+fuente)), data=datos, control=control.lme, method="REML")
lme.4c <- lme(a_por_b ~ t_sonda+fuente*especie, random=list(id=pdIdent(~1+t_sonda+fuente)), data=datos, control=control.lme, method="REML")

## modelo de FIXED INTERCEPT & RANDOM SLOPES WITH CORRELATED RANDOM SLOPES FOR DIFFERENT COVARIATES
##   lmer.5 "Satterthwaite" es aproximadamente como lme.5a pdSymm
lmer.5 <- lmer(a_por_b ~ t_sonda+fuente*especie + (t_sonda+fuente-1|id), data=datos, control=control.lmer, REML=TRUE)
##
lme.5 <- lme(a_por_b ~ t_sonda+fuente*especie, random=~t_sonda+fuente-1|id, data=datos, control=control.lme, method="REML")   ## identico a la linea siguiente
lme.5 <- lme(a_por_b ~ t_sonda+fuente*especie, random=list(id=~t_sonda+fuente-1), data=datos, control=control.lme, method="REML")
lme.5a <- lme(a_por_b ~ t_sonda+fuente*especie, random=list(id=pdSymm(~t_sonda+fuente-1)), data=datos, control=control.lme, method="REML")
lme.5b <- lme(a_por_b ~ t_sonda+fuente*especie, random=list(id=pdCompSymm(~t_sonda+fuente-1)), data=datos, control=control.lme, method="REML")



## modelo de EFECTOS ALEATORIOS ANIDADOS: RAMDOM INTERCEPT & FIXED SLOPES
lmer.anid.1 <- lmer(a_por_b ~ t_sonda+fuente + (1|especie/id), data=datos, control=control.lmer, REML=TRUE)
##    lme no puede gestionar, de manera directa, dos efectos aleatorios cruzados
##    en los modelos lme importa el orden de los efectos aleatorios: especie/id implica -> especie, id; especie/id no implica -> id, especie
##      los dos modelos lme.anid.1 son equivalentes; los modelos 1 y 1a proporcionan identicos resultados
##    lmer.anid.1 "Satterthwaite" es aproximadamente como los dos modelos lme(...)
lme.anid.1 <- lme(a_por_b ~ t_sonda+fuente, random=~1|especie/id, data=datos, control=control.lme, method="REML")
lme.anid.1 <- lme(a_por_b ~ t_sonda+fuente, random=list(especie=~1, id=~1), data=datos, control=control.lme, method="REML")
lme.anid.1a <- lme(a_por_b ~ t_sonda+fuente, random=list(especie=pdSymm(~1), id=pdSymm(~1)), data=datos, control=control.lme, method="REML")

## modelo de EFECTOS ALEATORIOS ANIDADOS: CORRELATED RANDOM INTERCEPTS & SLOPES FOR TWO CORRELATED FIXED EFFECTS
##    lmer.anid.2 "Satterthwaite" es aproximadamente como lme(...) lme.anid.2a pdSymm
lmer.anid.2 <- lmer(a_por_b ~ t_sonda+fuente + (1+t_sonda+fuente|especie/id), data=datos, control=control.lmer, REML=TRUE)
## 
lme.anid.2a <- lme(a_por_b ~ t_sonda+fuente, random=list(especie=pdSymm(~1+t_sonda+fuente), id=pdSymm(~1+t_sonda+fuente)), data=datos, control=control.lme, method="REML")
lme.anid.2b <- lme(a_por_b ~ t_sonda+fuente, random=list(especie=pdCompSymm(~1+t_sonda+fuente), id=pdCompSymm(~1+t_sonda+fuente)), data=datos, control=control.lme, method="REML")

## random=list(especie=pdSymm(~1+t_sonda+fuente), id=pdSymm(~1+t_sonda+fuente)) especifica una estructura 
##   de correlacion simetrica (pdSymm) para los efectos aleatorios "especie" e "id". 
##   Se permite que las pendientes aleatorias para t_sonda y fuente dentro de cada "especie" o "id" esten correlacionadas. 
##   El efecto de t_sonda puede depender del nivel de fuente, y viceversa, dentro de cada "especie" o "id". 
##   Tambien permite la correlacion entre el intercepto aleatorio y las pendientes dentro de cada "especie" o "id".
## random=list(especie=~1+t_sonda+fuente, id=~1+t_sonda+fuente) especifica una estructura de efectos aleatorios general y no estructurada, 
##   siendo mas flexible porque permite una estructura de covarianza más general para los efectos aleatorios.
##   Las pendientes aleatorias para t_sonda y fuente dentro de cada "especie" o "id" no están obligadas a estar simetricamente correlacionadas, 
##   aunque permite la correlación entre el intercepto aleatorio y las pendientes dentro de cada "especie" o "id", 
##   pero la correlación entre las pendientes t_sonda y fuente puede no ser necesariamente simetrica.
## Si no tenemos razones de peso para suponer una correlacion simetrica, el modelo mas flexible (lme.anid.2) podria ser una opcion mas segura.
##   Pero tambien podria ser mas complejo y dificil de ajustar, pudiendo dar errores (system is computationally singular),
##   o podria conllevar el riesgo de sobreajustarse si no disponemos de muchos datos.

## modelo de EFECTOS ALEATORIOS ANIDADOS: CORRELATED RANDOM SLOPES AND INTERCEPTS FOR TWO UNCORRELATED FIXED EFFECTS
##    lmer.anid.3 "Satterthwaite" es aproximadamente como lme(...) lme.anid.3a pdSymm
lmer.anid.3 <- lmer(a_por_b ~ t_sonda+fuente + (1+t_sonda|especie/id)+(1+fuente|especie/id), data=datos, control=control.lmer, REML=TRUE)
## 
lme.anid.3a <- lme(a_por_b ~ t_sonda+fuente, random=list(especie=pdSymm(~1+t_sonda), especie=pdSymm(~1+fuente), id=pdSymm(~1+t_sonda), id=pdSymm(~1+fuente)), data=datos, control=control.lme, method="REML")
lme.anid.3b <- lme(a_por_b ~ t_sonda+fuente, random=list(especie=pdCompSymm(~1+t_sonda), especie=pdCompSymm(~1+fuente), id=pdSymm(~1+t_sonda), id=pdSymm(~1+fuente)), data=datos, control=control.lme, method="REML")

## modelo de EFECTOS ALEATORIOS ANIDADOS: UNCORRELATED RANDOM SLOPES AND INTERCEPTS FOR TWO UNCORRELATED FIXED EFFECTS
##   todos los modelos lme(...) son casi equivalentes en sus resultados, y convergentes con lmer.anid.4 "Satterthwaite"
lmer.anid.4 <- lmer(a_por_b ~ t_sonda+fuente + (1|especie/id)+(t_sonda-1|especie/id)+(fuente-1|especie/id), data=datos, control=control.lmer, REML=TRUE)
## 
lme.anid.4 <- lme(a_por_b ~ t_sonda+fuente, random=list(especie=~1, especie=~t_sonda-1, especie=~fuente-1, id=~1, id=~t_sonda-1, id=~fuente-1), data=datos, control=control.lme, method="REML")
lme.anid.4a <- lme(a_por_b ~ t_sonda+fuente, random=list(especie=pdSymm(~1), especie=pdSymm(~t_sonda-1), especie=pdSymm(~fuente-1), id=pdSymm(~1), id=pdSymm(~t_sonda-1), id=pdSymm(~fuente-1)), data = datos, control=control.lme, method="REML")
lme.anid.4b <- lme(a_por_b ~ t_sonda+fuente, random=list(especie=pdDiag(~1+t_sonda+fuente), id=pdDiag(~1+t_sonda+fuente)), data=datos, control=control.lme, method="REML")   ## difiere de los anteriores en los grados de libertad
lme.anid.4c <- lme(a_por_b ~ t_sonda+fuente, random=list(especie=pdBlocked(list(~1, ~t_sonda-1, ~fuente-1)), id= pdBlocked(list(~1, ~t_sonda-1, ~fuente-1))), data=datos, control=control.lme, method="REML")



#### comparacion de modelos con IDENTICA PARTE FIJA y diferente aleatoria ####
## A la hora de comparar diferentes modelos anidados, el metodo REML lo podremos utilizar si, y solo si,
##   se mantienen constantes los efectos fijos!!! (sean factores o covariantes) y SOLO cambia la estructura de los efectos aleatorios.
##   NUNCA deberemos utilizar el metodo REML si en los modelos anidados que se comparan cambian los efectos fijos.
## Podemos actualizar los modelos usando ML o REML de la siguiente manera:
##     update(modelo.lmer, REML=TRUE)  o  update(modlme, method="REML")
## AICc modelos con un factor aleatorio 
## pesos de los modelos con un factor aleatorio 
AICc(lmer.1, lmer.2, lmer.3, lmer.4, lmer.5)
Weights(AICc(lmer.1, lmer.2, lmer.3, lmer.4, lmer.5))
AICc(lme.1, lme.1a, lme.1b, lme.1c, lme.2, lme.2a, lme.2b, lme.3, lme.3a, lme.3b, lme.4, lme.4a, lme.4b, lme.4c, lme.5, lme.5a, lme.5b)
Weights(AICc(lme.1, lme.1a, lme.1b, lme.1c, lme.2, lme.2a, lme.2b, lme.3, lme.3a, lme.3b, lme.4, lme.4a, lme.4b, lme.4c, lme.5, lme.5a, lme.5b))

## AICc modelos con dos factores aleatorios anidados
AICc(lmer.anid.1, lmer.anid.2, lmer.anid.3, lmer.anid.4)
Weights(AICc(lmer.anid.1, lmer.anid.2, lmer.anid.3, lmer.anid.4))
AICc(lme.anid.1, lme.anid.1a, lme.anid.2a, lme.anid.2b, lme.anid.3a, lme.anid.3b, lme.anid.4, lme.anid.4a, lme.anid.4b, lme.anid.4c)
Weights(AICc(lme.anid.1, lme.anid.1a, lme.anid.2a, lme.anid.2b, lme.anid.3a, lme.anid.3b, lme.anid.4, lme.anid.4a, lme.anid.4b, lme.anid.4c))

## cuantas veces es mejor un modelo que otro? (primero el de menor valor de AICc)
exp(-0.5 * (AICc(lmer.1) - AICc(lmer.2)))   
exp(-0.5 * (AICc(lmer.anid.1) - AICc(lmer.anid.3)))

## comparamos modelos con la misma parte fija y diferente parte aleatoria
##   usando estimas ML y no-REML; con modelos lme(...) hay que hacerlo manualmente
anova(lmer.1, lmer.2)
anova(update(lme.1, method="ML"), update(lme.2, method="ML"))
anova(lmer.1, lmer.4)
anova(update(lme.1, method="ML"), update(lme.4, method="ML"))

## equivalente a utilizar likelihood ratio tests (usando estimas ML y no-REML)
lrtest(update(lmer.1, REML=F), update(lmer.2, REML=F))
lrtest(update(lme.1, method="ML"), update(lme.2, method="ML"))



#### VISUALIZACION DE MATRICES DE VARIANZA-COVARIANZA PARA LOS EFECTOS ALEATORIOS ####
heatmap(getVarCov(lme.2), margins=c(10,10), revC=TRUE)
print(getVarCov(lme.2))

## VARIANZAS de los terminos aleatorios
summary(lmer.2)$varcor  ## modelo lmer
VarCorr(lme.2)          ## modelo lme



## elegimos un modelo
modelo <- lmer.1
modlme <- lme.1

## podemos crear el mismo modelo seleccionado utilizando la funcion glmmTMB
##   usa la misma nomenclatura para la formula que los modelos lmer
control.tmb <- glmmTMBControl(optCtrl=list(iter.max=1000, eval.max=1000), parallel=10)
modtmb <- glmmTMB(formula(modelo), data=datos, dispformula= ~1, REML=TRUE, control=control.tmb)

## para sacar la formula de lo que hemos hecho:
modelo@call; formula(modelo)
modlme$call; formula(modlme)
modtmb$call; formula(modtmb)

## los datos con los que hemos trabajado estan ahora en (respuesta y luego las columnas de las predictoras):
modelo@frame
modlme$data   ## con modelos lme aparecen todos los datos del dataframe usado
modtmb$frame



#### RESIDUOS DEL MODELO ####
## lecturas: https://arxiv.org/pdf/1502.06988.pdf
##           https://besjournals.onlinelibrary.wiley.com/doi/pdfdirect/10.1111/2041-210X.13434
##           https://peerj.com/articles/9522.pdf
##           https://link.springer.com/article/10.3758/s13428-021-01587-5

## Exploracion visual de la normalidad de los residuos clasicos del modelo
##   histograma modelo lmer
hist(residuals(modelo), density=10, freq=FALSE, main="residuos del modelo lmer", ylim=c(0,0.3), lwd=2, col="black")
curve(dnorm(x, mean=mean(residuals(modelo)), sd=sd(residuals(modelo))), col="red", lwd=2, add=TRUE, yaxt="n")
##   histograma modelo lme
hist(residuals(modlme), density=10, freq=FALSE, main="residuos del modelo lme", ylim=c(0,0.3), lwd=2, col="black")
curve(dnorm(x, mean=mean(residuals(modlme)), sd=sd(residuals(modlme))), col="red", lwd=2, add=TRUE, yaxt="n")
##   histograma modelo glmmTMB
hist(residuals(modtmb), density=10, freq=FALSE, main="residuos del modelo lme", ylim=c(0,0.3), lwd=2, col="black")
curve(dnorm(x, mean=mean(residuals(modtmb)), sd=sd(residuals(modtmb))), col="red", lwd=2, add=TRUE, yaxt="n")
##
##   qqplot modelo lmer
qqmath(modelo, id=0.05)  
qqnorm(residuals(modelo), main="residuos del modelo")
qqline(residuals(modelo), col="red", lwd=2)
##   qqplot modelo lme
qqnorm(residuals(modlme), main="residuos del modelo")
qqline(residuals(modlme), col="red", lwd=2)
##   qqplot modelo glmmTMB
qqnorm(residuals(modtmb), main="residuos del modelo")
qqline(residuals(modtmb), col="red", lwd=2)
##
## test de normalidad
## https://en.wikipedia.org/wiki/Shapiro%E2%80%93Wilk_test 
##  En vez del test de Kolmogorov-Smirnov (que asume que muestra y poblacion son coincidentes)
##  (en la mayoria de los casos tenemos una muestra que no coincide con la poblacion)
##  https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test 
shapiro.test(residuals(modelo))   ##   modelo lmer
shapiro.test(residuals(modlme))   ##   modelo lme
shapiro.test(residuals(modtmb))   ##   modelo glmmTMB
##  cuidado con este test: http://www.statisticalmisses.nl/index.php/frequently-asked-questions/77-what-is-wrong-with-tests-of-normality
##
## parametrizacion de la normalidad de los residuos con lmer, lme y glmmTMB
##  kurtosis: puntiagudez de la distribucion (Ho= 3 o 0)
##    https://en.wikipedia.org/wiki/Kurtosis
##    mayor efecto sobre las P's en la violacion de la normalidad                                     
##    si K>0 (leptokurtosis) hay un mayor error tipo II (aceptar la Ho [nula] cuando de hecho es falsa) 
##    si K<0 (platikurtosis) hay un mayor error tipo I (rechazar la Ho [nula] cuando es cierta)
kurtosis(residuals(modelo))            ## kurtosis de Pearson, Ho = 3
anscombe.test(residuals(modelo)) 
##
##  sesgo: simetria de la distribucion (Ho= 0)
##    https://en.wikipedia.org/wiki/Bias
##    poco efecto sobre las P's en la violacion de la normalidad
##    leve aumento del error de tipo I (rechazar la Ho [nula] cuando es cierta)
skewness(residuals(modelo))            ## sesgo Ho = 0
agostino.test(residuals(modelo))
##
## experimento sobre la emergencia de resultados significativos de desvio de la normalidad por puro azar
##    generamos una normal (kk) que emula los residuos del modelo
kk <- rnorm(n=dim(modelo@frame)[1], mean=modelo@frame[,1], sd=sd(modelo@frame[,1]))  
skewness(kk) 
agostino.test(kk)

## "outliers"; con {performance} solo para modelos lmer
##   method="zscore", threshold=list('zscore'=2.58)
check_outliers(modelo, method="ci", threshold=list('ci'=0.99))
check_outliers(modelo, method="zscore", threshold=list('zscore'=2.58))
plot(check_outliers(modelo, method="zscore", threshold=list('zscore'=2.58)))
plot(check_outliers(modelo, method="cook", threshold=list('cook'=1)))
##
## otra opcion mas elaborada: distancia de Cook y leverage
win.graph()
plot(hatvalues(modelo), cooks.distance(modelo))
abline(h=4/length(modelo@frame[,1]), col="red")
abline(v=2*length(fixef(modelo))/length(residuals(modelo)), col="red")
identify(hatvalues(modelo), cooks.distance(modelo))     ## si no se ven dos lineas rojas es que no se supera un umbral de preocupacion
dev.off()
##
## otra opcion visual con el paquete {car}
##   nos fijamos en las observaciones con residuos studentizados (aprox. Student-t) >2 y valores altos de leverage y distancia de Cook
influencePlot(modelo, main="distancia de Cook: tamannio del circulo")   ## con un listado de las observaciones potencialmente problematicas

## heterocedasticidad de los residuos
##  con lmer
plot(fitted(modelo), residuals(modelo), xlab="valores predichos", ylab="residuos", main="HAY HETEROCEDASTICIDAD?")
abline(h=0, lty=2, lwd=2, col="red")
lines(smooth.spline(fitted(modelo), residuals(modelo), spar=0.95), col="blue")
##  con glmmTMB
plot(fitted(modtmb), residuals(modelo), xlab="valores predichos", ylab="residuos", main="HAY HETEROCEDASTICIDAD?")
abline(h=0, lty=2, lwd=2, col="red")
lines(smooth.spline(fitted(modtmb), residuals(modtmb), spar=0.95), col="blue")
##
## otro modo mas directo
plot(modelo, id=0.05)   ## con lmer
plot(modlme, id=0.05)   ## con lme 

## ANALISIS GLOBAL DE RESIDUOS de manera sintetica para la revision de los modelos lmer:
##   residuos clasicos
##     con {LMERConvenienceFunctions}
mcp.fnc(modelo)                 
##     con {performance}   *** abrir mucho el panel de Plots ***
check_model(modelo)             
check_model(modlme)   ## solo para modelos sencillos lme(...) de random intercepts
check_model(modtmb)

## con modelos lmer y glmmTMB es posible efectuar una diagnosis utilizando residuos ordinales de DHARMa (escalados entre 0-1)
##   valoramos la correspondencia entre residuos clasicos y los de DHARMA
##     elegimos uno de los modelos mixtos lmer o glmmTMB
modelo.mixto <- modelo
residuos.mixto <- simulateResiduals(modelo.mixto)
plot(residuos.mixto$scaledResiduals, residuals(modelo.mixto))
##
## PREDICCIONES DEL MODELO lmer y EXPLICACION de los residuos de DHARMA
plot(fitted(modelo), modelo@frame[,1]); abline(a=0, b=1, col="red", lwd=2)
predict.lmer <- predictInterval(modelo, newdata=datos, level=0.999)
predicciones <- data.frame(modelo@frame[1], predict.lmer[3], predict.lmer[2], predict.lmer[1])
predicciones$residuos <- predicciones[1] - predicciones[4]
predicciones$residsDHARMA <- residuos.mixto$scaledResiduals
predicciones

##   normalidad, heterocedasticidad de los residuos y existencia de outliers
testUniformity(modelo.mixto)
testQuantiles(modelo.mixto)
testOutliers(modelo.mixto)
residuos_ordinales <- residuos.mixto$scaledResiduals
which(residuos_ordinales<0.005 | residuos_ordinales>0.995)  ## con valores predefinidos escalados extremos (aqui 1% mas extremos)
##   sinteticamente
plot(simulateResiduals(modelo.mixto))

## HOMOGENEIDAD DE VARIANZAS ENTRE LOS GRUPOS DE LOS PREDICTORES NOMINALES (factores)
## con el paquete {DHARMa} para los residuos escalados de modelos lmer y glmmTMB
##     elegimos uno de los modelos mixtos lmer o glmmTMB
modelo.mixto <- modelo
testCategorical(modelo.mixto, datos$fuente)                   ## *** ojo !!  poned aqui VUESTROS DATOS ***
testCategorical(modelo.mixto, datos$especie)                  ## *** ojo !!  poned aqui VUESTROS DATOS ***
testCategorical(modelo.mixto, datos$fuente:datos$especie)     ## *** ojo !!  poned aqui VUESTROS DATOS ***

## analisis externo de la HOMOGENEIDAD DE VARIANZAS ENTRE LOS GRUPOS DE LOS PREDICTORES NOMINALES
##   (utilizando las variables de interes)
##   varianzas de las diferentes celdas de interaccion entre factores
##   *** ojo !!  poned aqui VUESTROS DATOS ***
miformula <- as.formula(a_por_b ~ fuente*especie)
aggregate(miformula, data=datos, FUN=var)     ## varianza
aggregate(miformula, data=datos, FUN=mean)    ## media
aggregate(miformula, data=datos, FUN=length)  ## tamannio muestral
##
## disponemos de diferentes tests
##   The Levene test works well in the ANOVA framework, providing there are small to moderate deviations from the normality.
##   In this case, it outperfoms the Bartlett test. If the distribution are nearly normal, however, the Bartlett test is better.
##   The Fligner-Killeen test is to be prefered in case of strong departure from the normality (to which the Bartlett test is sensible).
##   En el test de Levene otra opcion es center=mean (e.g., la del programa STATISTICA)
##     consultad: http://www.cookbook-r.com/Statistical_analysis/Homogeneity_of_variance/
## para un solo factor
leveneTest(a_por_b ~ especie, data=datos, center=median)     ## *** ojo !!  poned aqui VUESTROS DATOS ***
leveneTest(a_por_b ~ especie, data=datos, center=mean)       ## *** ojo !!  poned aqui VUESTROS DATOS ***
bartlett.test(a_por_b ~ especie, data=datos)                 ## *** ojo !!  poned aqui VUESTROS DATOS ***
fligner.test(a_por_b ~ especie, data=datos)                  ## *** ojo !!  poned aqui VUESTROS DATOS ***
## para dos factores o mas
leveneTest(a_por_b ~ especie:fuente, data=datos, center=median)      ## *** ojo !!  poned aqui VUESTROS DATOS ***
leveneTest(a_por_b ~ especie:fuente, data=datos, center=mean)        ## *** ojo !!  poned aqui VUESTROS DATOS ***
bartlett.test(a_por_b ~ interaction(especie, fuente), data=datos)    ## *** ojo !!  poned aqui VUESTROS DATOS ***
fligner.test(a_por_b ~ interaction(especie, fuente), data=datos)     ## *** ojo !!  poned aqui VUESTROS DATOS ***
##
## relacion entre medias y varianzas
## si hay heterocedasticidad de los residuos, que no haya una clara relacion en lo siguiente
##    si la relacion es positiva, aumenta el error de tipo I
##    si la relacion es negativa, aumenta el error de tipo II
## *** ojo !!  poned aqui VUESTROS DATOS ***
miformula <- as.formula(a_por_b ~ fuente*especie)
num.col <- dim(aggregate(miformula, data=datos, FUN=mean))[2]
medias.originales <- aggregate(miformula, data=datos, FUN=mean)[num.col]
varianzas.originales <- aggregate(miformula, data=datos, FUN=var)[num.col]
cor.media.varianza <- round(cor(medias.originales[,1], varianzas.originales[,1]), 3)
plot(medias.originales[,1], varianzas.originales[,1], ylab="VARIANZAS DE LAS CELDAS", xlab="MEDIAS DE LAS CELDAS", main=paste("r entre medias y varianzas =",cor.media.varianza))
abline(lm((varianzas.originales[,1]) ~ medias.originales[,1]), col="red", lwd=2)
par(mfcol=c(1,1))



#### MULTICOLINEARIDAD ENTRE PREDICTORES ####
## para modelos lmer
check_collinearity(modelo)
plot(check_collinearity(modelo))
## para modelos lme solo para random intercept models
check_collinearity(modlme)
plot(check_collinearity(modlme))
## para modelos lmer
check_collinearity(modtmb)
plot(check_collinearity(modtmb))



#### CONTROL DE LA HETEROCEDASTICIDAD RESIDUAL ####
## podemos abordar la heterocedasticidad de los residuos usando modelos mixtos lme(...)   no con lmer(...) !!!
##   modelo Heteroskedasticity-Corrected  (hay otras opciones: mirad en Help varPower, varExp, varIdent, varFunc)
##   https://quantdev.ssri.psu.edu/sites/qdev/files/ILD_Ch06_2017_MLMwithHeterogeneousVariance.html
plot(modlme, main="modelo original")
modlme.HC0 <- update(modlme, weights=varPower())
plot(modlme.HC0, main="modelo HC0")
modlme.HC1 <- update(modlme, weights=varExp(form=~t_sonda))
plot(modlme.HC1, main="modelo HC1")
modlme.HC2 <- update(modlme, weights=varExp(form=~t_sonda|id))
plot(modlme.HC2, main="modelo HC2")
modlme.HC3 <- update(modlme, weights=varFixed(~t_sonda))
plot(modlme.HC3, main="modelo HC3")
modlme.HC4 <- update(modlme, weights=varIdent(~1|id))
plot(modlme.HC4, main="modelo HC4")
modlme.HC5 <- update(modlme, weights=varComb(varIdent(form=~1|id), varExp(form=~t_sonda)))
plot(modlme.HC5, main="modelo HC5")
## comparamos con el modelo lme original frente a 
##   los modelos que en los plots ya no muestran heterocedasticidad residual
AICc(modlme, modlme.HC0, modlme.HC2, modlme.HC5)
##
## ahora analizamos los resultados del modelo mas eficaz y plausible
##   frente al original con el problema; vemos como se modifican los coeficientes, std.err, P's ...
##   elegimos uno de los modelos previos
modlme.HC <- modlme.HC0
##
options(scipen=999)
lrtest(update(modlme, method="ML"), update(modlme.HC, method="ML"))
##
round(summary(modlme.HC)$tTable, 3)
round(summary(modlme)$tTable, 3)
##
anova.lme(modlme.HC, type="marginal")
anova.lme(modlme, type="marginal")
## 
plot(allEffects(modlme.HC), rotx=45) 
plot(allEffects(modlme), rotx=45) 

## con modelos glmmTMB
##  aqui va un ajemplo *** PONED AQUI VUESTROS DATOS *** dentro de dispformula= ~...
modtmb.HC <- update(modtmb, dispformula=~fuente)
testQuantiles(modtmb.HC)
AICc(modtmb, modtmb.HC)
##
lrtest(update(modtmb, REML=FALSE), update(modtmb.HC, REML=FALSE))
##
round(summary(modtmb.HC)$coefficients$cond, 3)
round(summary(modtmb)$coefficients$cond, 3)
##
Anova(modtmb.HC, type=3)
Anova(modtmb, type=3)
## 
plot(allEffects(modtmb.HC), rotx=45) 
plot(allEffects(modtmb), rotx=45) 



#### OMNIBUS TESTS de la "significacion" global del modelo ####
## Utilizamos la aproximacion de los "likelihood ratio tests" (lrtest) utilizando Maximum Likelihood (ML)
##   https://en.wikipedia.org/wiki/Likelihood-ratio_test
## Se comparan modelos anidados, como son el modelo de interes y el MODELO NULO 
##   el MODELO NULO es el que tiene la misma parte aleatoria, pero en la fija no tiene ningun efecto (solo el intercepto)
## Tambien podriamos comparar dos modelos de interes que tuviesen diferentes subconjuntos de efectos fijos, y misma parte aleatoria

##   modelos lmer
modelo.ML <- update(modelo, REML=FALSE)
modelo.nulo.ML <- update(modelo.ML, .~1+(1|id), REML=FALSE)  ## *** ojo !!  poned aqui VUESTROS DATOS *** con la misma parte aleatoria
##   modelos lme
modlme.ML <- update(modlme, method="ML")
modlme.nulo.ML <- update(modlme.ML, .~1, random=list(id=~1), method="ML")  ## *** ojo !!  poned aqui VUESTROS DATOS *** con la misma parte aleatoria
##   modelos glmmTMB
modtmb.ML <- update(modtmb, REML=FALSE)
modtmb.nulo.ML <- update(modtmb.ML, .~1+(1|id), REML=FALSE)  ## *** ojo !!  poned aqui VUESTROS DATOS *** con la misma parte aleatoria

## Realizamos el likelihood ratio test entre el modelo de interes y el que carece de efectos fijos
##   ponemos antes el modelo mas sencillo (i.e., menos grados de libertad, Df); tambien podemos hacerlo con anova(..., ...)
lrtest(modelo.nulo.ML, modelo.ML)
lrtest(modlme.nulo.ML, modlme.ML)
lrtest(modtmb.nulo.ML, modtmb.ML)

## Realizamos la comparacion entre el modelo de interes y el que carece de efectos fijos utilizando Akaike AICc
AICc(modelo.nulo.ML, modelo.ML)
AICc(modlme.nulo.ML, modlme.ML)
AICc(modtmb.nulo.ML, modtmb.ML)

## Aproximacion basada en metodos de remuestreo para comparar modelos anidados. Solo para modelos lmer(...).
##   ponemos antes el modelo mas complejo y luego el mas sencillo
##   con nsim definimos el numero de procesos de bootstrapping
##   con seed podemos cambiar los procesos de aleatorizacion-remuestreo
PBmodcomp(modelo.ML, modelo.nulo.ML, nsim=1000, seed=123, details=0)  ## muy lento!!!



#### SIGNIFICACION DE EFECTOS DEL MODELO ####
## nos fijamos en los coeficientes de los efectos fijos 
##  para modelos lmer con aproximacion "Satterthwaite" (por defecto) y "Kenward-Roger"
round(summary(modelo, ddf="Kenward-Roger")$coefficients, 5)
model_parameters(modelo, ci=0.95, ci_method="kenward")
round(summary(modelo)$coefficients, 5)     ## lo mismo que anniadiendo ddf="Satterthwaite"
model_parameters(modelo, ci=0.95, ci_method="satterthwaite")
##  para los modelos lme son los grados de libertad marginales
round(summary(modlme)$tTable, 5)
##  para modelos glmmTMB con la aproximacion asintotica de Wald
round(summary(modtmb)$coefficients$cond, 5)
model_parameters(modtmb, ci=0.95, ci_method="wald")

## COEFCIENTES DE REGRESION ESTANDARIZADOS PARA LAS COVARIANTES (predictores continuos)
## zeta-standarizamos (a MEDIA=0 Y SD=1) 
## para obtener los coeficientes estandarizados del modelo (hace comparables los coeficientes para valorar su importancia relativa EN VALOR ABSOLUTO)
## modelos lmer  
standardize_parameters(modelo)         ## coeficientes de regresion para la parte fija ESTANDARIZADOS
plot(standardize_parameters(modelo))
fixef(modelo)                          ## coeficientes de regresion para la parte fija SIN ESTANDARIZAR
## modelos glmmTMB  
standardize_parameters(modtmb)         ## coeficientes de regresion para la parte fija ESTANDARIZADOS
plot(standardize_parameters(modtmb))
fixef(modtmb)                          ## coeficientes de regresion para la parte fija SIN ESTANDARIZAR


## TEST DE EFECTOS ALEATORIOS: valores de Chi2 y sus p-values para los efectos aleatorios
## modelos lmer
rand(modelo)
## modelos lme comparando el modelo SIN y CON efectos aleatorios
anova(gls(formula(modlme), data=datos, method="REML"), modlme)
## modelos glmmTMB: comparando el modelo SIN y CON efectos aleatorios
formula(modtmb)
anova(update(modtmb, a_por_b ~ t_sonda + fuente * especie), modtmb)     ## *** PONED AQUI VUESTROS DATOS ***

## VARIANZAS de los terminos aleatorios
summary(modelo)$varcor    ## modelo lmer
VarCorr(modlme)           ## modelo lme
print(getVarCov(modlme))  ## modelo lme 
summary(modtmb)$varcor    ## modelo glmmTMB

## VALORES de los efectos aleatorios
qqmath(ranef(modelo))    ## solo para modelos con un unico termino aletorio
plot(ranef(modlme))
ranef(modtmb)



## TEST DE EFECTOS FIJOS
## opcion ideal cuando tengamos como predictores a factores con tres niveles o mas
## significaciones de los Fixed Effects usando la aproximacion Kenward-Roger
##    proporciona la Sum Sq de los efectos fijos
##    los DenDf son los mismos que los obtenidos en df de summary
##    los valores de la F son identicos a los cuadrados de los t value de las salidas summary
##    los valores de significaciones (p) son identicos a los obtenidos en summary
##    no proporciona los interceptos
## Kenward-Roger modification of the F-statistic for some linear mixed models fitted with lmer
## http://web.warwick.ac.uk/statsdept/useR-2011/abstracts/290311-halekohulrich.pdf
## utilizamos una version de "anova" (anova.merModLmerTest) aplicada a objetos merMod, no equivalente a la usada en objetos lm o glm
##    modelos lmer
anova(modelo, ddf="Kenward-Roger", type=3)
anova(modelo, ddf="Satterthwaite", type=3)
##    modelos lme
anova.lme(modlme, type="marginal")[-1,]
##    modelos glmmTMB (estima asintotica)
Anova(modtmb, type=3)
##
## con el comando Anova de {car} se obtienen directamente los resultados de significacion usando la aproximacion Kenward-Roger
## proporciona las significaciones para el intercepto
## los resultados de F, p y Df son identicos a los de anova(..., ddf="Kenward-Roger", type=3)
Anova(modelo, type=3, test="F")
## no lo aplicamos a modelos lmer o lme con el test="Chisq"; produce resultados erroneos



## INTERVALOS DE CONFIANZA DE LOS COEFICIENTES
## para modelos lmer
confint(modelo, level=0.95)
confint.merMod(modelo, parm="beta_", level=0.95, method="profile")
model_parameters(modelo, ci=0.95, ci_method="wald")
model_parameters(modelo, ci=0.95, ci_method="kenward")
model_parameters(modelo, ci=0.95, ci_method="satterthwaite")
## para modelos lme
intervals(modlme, levels=0.95)
## para modelos glmmTMB
model_parameters(modtmb, ci=0.95, ci_method="wald")



#### BOOTSTRAPPING DEL MODELO ####
## Podemos efectuar parametrizaciones de nuestro modelo mixto lmer utilizando
##   procedimientos de remuestreo (bootstrapping), que tambien seran de utilidad 
##   para obtener, de otro modo, estimas de significacion (si el intervalo no incluye el valor cero).
## Para cuando nuestro modelo no proporciona significaciones.
##   cambiamos level=0.99 para p=0.01. Mas rapido con la estima asintetica de "profile": confint.merMod(modelo, parm="beta_", level=0.95, method="profile")
confint.fixed <- confint.merMod(modelo, parm="beta_", level=0.95, method="boot", boot.type="perc", nsim=500, verbose=TRUE)
print(confint.fixed, digits=3)
## para modelos lme en el caso de que haya un solo termino en random=list(...); es igual que: 
intervals(modlme, levels=0.95)$fixed

## Procedimientos de remuestreo (BOOTSTRAPPING), que tambien seran de utilidad cuando nuestro modelo no proporciona significaciones.
##    para obtener, de otro modo, estimas de significacion (si el intervalo no incluye el valor cero).
##    ci_method puede ser "ETI", "HDI", "BCI"
##    es muy conveniente hacer 1000 procesos (iterations) o mas
##  modelos lmer  
boot.lmer <- bootstrap_parameters(modelo, ci=0.95, iterations=1000, test=c("p-value", "p_direction", "p_map"), ci_method="BCI")
boot.lmer
##  modelos glmmTMB (muy lento respecto a un modelo lmer)
boot.tmb <- bootstrap_parameters(modtmb, ci=0.95, iterations=1000, test=c("p-value", "p_direction", "p_map"), ci_method="BCI")
boot.tmb
##
## para modelos nlme::lme
##  suele dar problemas con modelos con factores fijos cruzados, algunos de los cuales tengan muchos niveles con pocas replicas



#### CUANTO EXPLICA EL MODELO? ####
## VALORES OBSERVADOS Y PREDICHOS
##   con modelos lmer
predicho <- fitted(modelo)
plot(modelo@frame[,1] ~ predicho, ylab="RESPUESTA", xlab="PREDICCIONES")
abline(lm(modelo@frame[,1] ~ predicho), col="green", lwd=2)
##   con modelos lme (no efectuamos esta estima con modelos lme corregidos por heterocedasticidad)
predicho <- fitted(modlme)
plot(get.response(modlme) ~ predicho, ylab="RESPUESTA", xlab="PREDICCIONES")
abline(lm(get.response(modlme) ~ predicho), col="green", lwd=2)
##   con modelos glmmTMB
predicho <- fitted(modtmb)
plot(modtmb$frame[,1] ~ predicho, ylab="RESPUESTA", xlab="PREDICCIONES")
abline(lm(modtmb$frame[,1] ~ predicho), col="green", lwd=2)

## PARTICION DE LA VARIANZA
## R2m: marginal R2 (the proportion of variance explained by the fixed factor(s) alone)
## R2c: conditional R2 (the proportion of variance explained by both the fixed and random factors; i.e. the entire model)
##  consultad: https://jonlefcheck.net/2013/03/13/r2-for-linear-mixed-effects-models/comment-page-1/
##             http://onlinelibrary.wiley.com/doi/10.1111/j.2041-210x.2012.00261.x/epdf
##             http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12225/epdf
r.squaredGLMM(modelo)   ## modelo lmer
r.squaredGLMM(modlme)   ## modelo lme  (no efectuamos esta estima con modelos lme corregidos por heterocedasticidad)
r.squaredGLMM(modtmb)   ## modelo glmmTMB   (no efectuamos esta estima con modelos lme corregidos por heterocedasticidad)

## MAGNITUD DE LOS EFECTOS FIJOS
## partial eta-squared se refiere a una magnitud relativa para cada predictor, estimada del siguiente modo:
##   Eta2 parcial = SS_efecto / (SS_efecto + SS_residual);  https://easystats.github.io/effectsize/articles/anovaES.html
##   SS_residual = sum(residuals(modelo)^2)
## tambien contamos con los metodos insesgados omega y epsilon squared; https://easystats.github.io/effectsize/articles/anovaES.html#other-measures-of-effect-size
##   Omega squared = [SS_efecto - (df_efecto * MS_residual)] / (SS_total + MS_residual)
##   Epsilon squared es como omega_squared pero asociado a adjusted-R2; denominado Adjusted Eta Squared
## con modelos lme4::lmer
eta_squared(modelo)    ## lo mismo que el siguiente
eta_squared(anova(modelo, ddf="Satterthwaite", type=3))
eta_squared(anova(modelo, ddf="Kenward-Roger", type=3))
omega_squared(anova(modelo, ddf="Kenward-Roger", type=3))
epsilon_squared(anova(modelo, ddf="Kenward-Roger", type=3))  
## con modelos nlme::lme
eta_squared(anova.lme(modlme, type="marginal"))
omega_squared(anova.lme(modlme, type="marginal"))
epsilon_squared(anova.lme(modlme, type="marginal"))
## con modelos glmmTMB no es posible abordarlo

## tabla sintetica de los resultados con la PARTICION DE LA VARIANZA DE LA RESPUESTA (y magnitudes de efectos parciales) en modelos lme4::lmer
##   utilizando la estima de grados de libertad de Kenward-Roger
##   prop_varianza es el tanto por uno de la variacion en la respuesta explicado exclusivamente por ese efecto fijo
SStotal <- sum((modelo@frame[,1]-mean(modelo@frame[,1]))^2)      
SSerror <- SStotal*(1-as.numeric(r.squaredGLMM(modelo)[2]))
tabla.SS <- as.data.frame(anova(modelo, ddf="Kenward-Roger", type=3))
prop_varianza <- tabla.SS[,1]/SStotal
tabla.SS <- data.frame(prop_varianza, tabla.SS)
names(tabla.SS)[7] <- "Pr.F"
print(tabla.SS, digits=4)
var.compartida <- (as.numeric(r.squaredGLMM(modelo)[1])-sum(tabla.SS$prop_varianza))
print(c("proporcion de la varianza atribuible a los efectos fijos =", round(as.numeric(r.squaredGLMM(modelo)[1]), 4)), quote=FALSE)
print(c("proporcion de la varianza compartida por los efectos fijos =", round(var.compartida, 4)), quote=FALSE)



#### TESTS POST-HOC ####
library(emmeans)     ## para tests post hoc
library(phia)        ## para medias de interacciones de factores
## otros metodos de ajuste son: "scheffe", "fdr", "holm"
##  *** PONED AQUI VUESTROS DATOS *** (nombre del factor en c("..."))
## modelo lmer
pairs(emmeans(modelo, c("especie")), adjust="tukey")
pairs(emmeans(modelo, c("fuente", "especie")), adjust="tukey")
## modelo lme
pairs(emmeans(modlme, c("especie")), adjust="tukey")
pairs(emmeans(modlme, c("fuente", "especie")), adjust="tukey")
## modelo glmmTMB
pairs(emmeans(modtmb, c("especie")), adjust="tukey")
pairs(emmeans(modtmb, c("fuente", "especie")), adjust="tukey")

## MEDIAS MARGINALES DE LOS NIVELES DE LOS FACTORES
##   (controlando por el efecto de las covariantes)
interactionMeans(modelo)   ## modelo lmer
interactionMeans(modlme)   ## modelo lme
## no es posible con un modelo lme
##
## otra posibilidad
##  *** PONED AQUI VUESTROS DATOS ***
emmeans(modelo, ~fuente*especie)   ## modelo lmer
emmeans(modlme, ~fuente*especie)   ## modelo lme
emmeans(modtmb, ~fuente*especie)   ## modelo glmmTMB

## tests post-hoc fuera del contexto del modelo mixto construido
library(ggstatsplot)
## unir dos factores en un factor interaccion
datos$interaccion <- with(data=datos, interaction(fuente,  posicion, sep="_"))
##
## plot
ggbetweenstats(
  data = datos,
  x = especie,                       ## factor cuyos niveles se desean comparar
  y = a_por_b,                       ## respuesta
  type = "parametric",               ## tipo de test tambien puede ser "nonparametric"
  plot.type = "box",
  pairwise.comparisons = TRUE,
  var.equal = FALSE,                 ## efectua la estima haciendo uso de los grados de libertad de Satterthwaite
  pairwise.display = "significant",
  p.adjust.method = "fdr",           ## tipo de correccion de multiples estimas de P
  effsize.type = "eta",
  centrality.plotting = FALSE,       ## proporciona media para tests parametricos
  bf.message = FALSE,                ## en tests parametricos proporciona el factor de Bayes en favor de Ho
  package = "pals",                  ## "ggsci" - "default_igv", "palettesForR" - "MATLAB"
  palette = "polychrome",            ## ved posibilidades con View(paletteer::palettes_d_names)
  title = "DIFERENCIA ENTRE ESPECIES",
  caption = "comentarios al grafico:"
)



#### VISUALIZACION DE EFECTOS ####
## valores medios para los niveles de los factores ajustados por las covariantes
##  *** PONED AQUI VUESTROS DATOS ***
## modelo lmer
plot(emmeans(modelo, ~fuente*especie))
## modelo lme
plot(emmeans(modlme, ~fuente*especie))
## modelo glmmTMB
plot(emmeans(modtmb, ~fuente*especie))
detach("package:emmeans")
detach("package:phia")

## Partial effects plots con {effects}
library(effects)
## elegimos un modelo lmer de interes que contenga a los predictores cuyos efectos parciales sobre la respuesta deseamos estimar
plot(allEffects(modelo), rotx=45, rug=F)    ## intervalos de confianza del 95%
##  para plots de efectos particulares de factores (si queremos ver los residuos parciales: partial.residuals=T)
plot(Effect("especie", partial.residuals=F, modelo), lwd=3, lty=1, rotx=45,
     partial.residuals=list(col="blue", pch=1, cex=1.25), rug=F)
plot(Effect(c("fuente", "especie"), partial.residuals=F, modelo), lwd=3, lty=1, rotx=45,
     partial.residuals=list(col="blue", pch=1, cex=1.25), rug=F)
##   para covariantes
plot(Effect("t_sonda", partial.residuals=F, modelo), lwd=3, lty=1, rotx=0,
     partial.residuals=list(col="blue", pch=1, cex=1.25), rug=F)
##
##   elegimos un modelo lme(...)
plot(allEffects(modlme), rotx=45, rug=F)    ## intervalos de confianza del 95%
##   para plots de efectos particulares
plot(Effect("especie", partial.residuals=F, modlme), lwd=3, lty=1, rotx=45,
     partial.residuals=list(col="blue", pch=1, cex=1.25), rug=F)
plot(Effect(c("fuente", "especie"), partial.residuals=F, modlme), lwd=3, lty=1, rotx=45,
     partial.residuals=list(col="blue", pch=1, cex=1.25), rug=F)
plot(Effect("t_sonda", partial.residuals=F, modlme), lwd=3, lty=1, rotx=0,
     partial.residuals=list(col="blue", pch=1, cex=1.25), rug=F)
## mirad la ayuda de "plot.effects" en Detailed Argument Descriptions
##
## elegimos un modelo glmmTMB de interes que contenga a los predictores cuyos efectos parciales sobre la respuesta deseamos estimar
plot(allEffects(modtmb), rotx=45, rug=F)    ## intervalos de confianza del 95%
##  para plots de efectos particulares de factores (si queremos ver los residuos parciales: partial.residuals=T)
plot(Effect("especie", partial.residuals=F, modtmb), lwd=3, lty=1, rotx=45,
     partial.residuals=list(col="blue", pch=1, cex=1.25), rug=F)
plot(Effect(c("fuente", "especie"), partial.residuals=F, modtmb), lwd=3, lty=1, rotx=45,
     partial.residuals=list(col="blue", pch=1, cex=1.25), rug=F)
##   para covariantes
plot(Effect("t_sonda", partial.residuals=F, modtmb), lwd=3, lty=1, rotx=0,
     partial.residuals=list(col="blue", pch=1, cex=1.25), rug=F)

## plot de valores parciales mediante el paquete {sjPlot}
##    type="int" para representar solo las interacciones en modelos factoriales
##    para mostrar ca. un error estandard en vez del intervalo al 95%: ci.lvl=0.68
##    "std2" muestra los efectos estandarizados tanto para covariantes como para los predictores binarios de los factores
## 
## plot de coeficientes estandarizados (comparables entre si en su magnitud); modelos lmer y glmmTMB (no lme)
plot_model(modelo, show.data=FALSE, type="std2", ci.lvl=0.95, terms=NULL)
plot_model(modelo, show.data=FALSE, type="std2", ci.lvl=0.95, terms=NULL, show.values=TRUE, value.offset=0.4)
##    ejemplos con variables concretas; modelos lmer, lme y glmmTMB
plot_model(modelo, show.data=FALSE, type="eff", ci.lvl=0.95, terms="t_sonda")
plot_model(modelo, show.data=FALSE, type="eff", ci.lvl=0.95, terms="especie")
plot_model(modelo, show.data=FALSE, type="eff", ci.lvl=0.95, terms=c("t_sonda", "fuente"))
plot_model(modelo, show.data=FALSE, type="eff", ci.lvl=0.95, terms=c("fuente", "especie"))
plot_model(modelo, show.data=FALSE, type="int", ci.lvl=0.95, terms=NULL)
plot_model(modelo, show.data=FALSE, type="est", ci.lvl=0.95, terms=NULL)
## para diagnosis del modelo; modelos lmer, lme y glmmTMB
plot_model(modelo, show.data=TRUE, type="resid", ci.lvl=0.95, terms=NULL)
plot_model(modelo, show.data=TRUE, type="diag", ci.lvl=0.95, terms=NULL)
##
detach(package:sjPlot)





####################################################################
#### ANALISIS PERMULACIONAL DEL MODELO BAJO AUSENCIA DE EFECTOS ####
####################################################################

## HIPOTESIS NULAS para los coeficientes con las que contrastar nuestro modelo mixto de interes lmer(...)  o  lme(...)
##   Randomization (Permutation) test under the null hypothesis that coefficients are cero.
##     Then we compare the observed coefficients in the model of interest with the 1000 permuted coefficient values.
##   The justification of the randomization test derives from the fact that under the null hypothesis of no treatment effect, 
##     the random assignment procedure produces a random shuffle of the responses. 
##   The justification of the permutation test derives from the fact that under the null hypothesis of identical distributions, 
##     all permutations of the responses are equally likely.
##   Randomization tests consider every possible permutation of the labels, permutation tests take a random sample of permutations of the labels
##   para ello aleatorizamos los valores de la variable respuesta por remuestreo-sin-reemplazo: sample(respuesta, replace=FALSE)
##   conveniente cuando ha habido desvios de los supuestos canonicos de los residuos del modelo
##     https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2687965/
##     https://link.springer.com/content/pdf/10.1007%2F978-1-4614-1365-3_10.pdf
##     https://lmcarrascal.eu/cursos/randtests.pdf
##     https://www.sciencedirect.com/topics/mathematics/randomization-test/pdf
##     https://en.wikipedia.org/wiki/Permutation_test
##     https://www.researchgate.net/publication/321385778_Randomization_Tests_or_Permutation_Tests_A_Historical_and_Terminological_Clarification
##     https://stats.stackexchange.com/questions/55742/difference-between-randomization-test-and-permutation-test
##     https://www.r-bloggers.com/2021/03/randomization-tests-make-fewer-assumptions-and-seem-pretty-intuitive/
##     https://uoftcoders.github.io/rcourse/lec12-randomization-tests.html
##     https://towardsdatascience.com/how-to-use-permutation-tests-bacc79f45749
##     https://measuringu.com/randomization-test/
##     https://deepblue.lib.umich.edu/bitstream/handle/2027.42/95955/oel_1.pdf

## El analisis de permutacion en estadistica es un tipo de metodo de remuestreo que implica reordenar (o permutar)
##   los datos observados para generar una distribucion de un estadistico dado bajo la hipotesis nula. 
##   La idea principal es evaluar los datos observados frente a una distribucion de lo que se esperaria por puro azar 
##   (hipotesis nula, Ho, de ausencia de efectos), suponiendo que Hoa es verdadera.
## Se utiliza a menudo para evaluar la significacion de parametros (e.g., coeficientes de regresión, valores F de los efectos fijos o R^2 del modelo)
##   derivados de un modelo estadistico especifico de interes. 
## El modo simplificado de como funciona es el siguiente:
##  1. Ajustamos tu modelo a los datos observados y calculamos los estadisticos de interes (por ejemplo, coeficientes, F de efectos o R^2).
##  2. Permutamos los datos originales reorganizando el orden de la variable de respuesta, mientras dejamos las variables predictoras sin cambios, 
##     y corremos el mismo modelo a esos datos permutados, calculando los estadisticos de interes.
##  3. Repetimos el paso 2 muchas veces para crear una distribucion de los estadisticos aleatorios bajo la hipotesis nula.
##  4. Comparamos nuestros estadisticos del modelo obtenidos con los datos reales con esas distribuciones de parametros asumiendo la Ho 
##     para calcular un valor frecuentista de P. Los valores de P seran las proporciones de la distribucion de permutacion que 
##     son mas extremas que los estadisticos observados en nuestro modelo original.
## Al comparar los parametros del modelo con los generados bajo la hipotesis nula Ho, podemos probar si nuestros resultados 
##     son probablemente debidos al azar, o si podrian indicar un efecto estadistico realmente existente que es muy imporbable 
##     que pudieran haberse obtenido por puro azar. 
## Las permutaciones estan destinadas a simular la aleatoriedad que esperariamos encontrar si la hipotesis nula Ho fuera verdadera. 
## Este es un metodo no parametrico, lo que significa que no se basa en supuestos sobre la distribucion de la poblacion subyacente.
## El analisis de permutacion es una buena alternativa para las situaciones en las que:
##  Los residuos violan algunas de las suposiciones de las pruebas parametricas (normalidad, homocedasticidad). 
##     El analisis permutacional no requiere que se cumplan estas suposiciones, porque genera una distribucion de estadisticos de interes
##     bajo la hipotesis nula basada directamente en nuestros datos a traves de permutaciones.
##  Los analisis permutacionales son tambien son bastante robustos a los valores atipicos (por ser outliers o tener altos valores de leverage), 
##     ya que toda la distribucion de permutaciones se genera a partir de esos datos observados reales, con lo que llevan implicitos esos valores extremos. 
##     Esto es una ventaja frente a las pruebas parametricas que pueden ser sensibles a los valores atipicos.
## Sin embargo, los analisis permutacionales no resuelven el problema de la multicolinealidad (alto VIF). 
##   Esto se debe a que el problema de la multicolinealidad involucra las relaciones entre nuestras variables predictoras, 
##   que no se ven modificadas por la permutacion de solamente la variable respuesta. 
##   Si la multicolinealidad es preocupante, es posible que necesitemos abordarla por separado (e.g., eliminando predictores altamente correlacionados),
##   pero el analisis permutacional no resolvera el problema.
## El analisis de permutacion no aborda directamente los problemas con la estructura de efectos aleatorios en un modelo mixto,
##   ya que solo implica reorganizar los datos de la variable respuesta para generar una distribucion nula. 
##   Problemas como la singularidad, la varianza muy pequennia (o nula) en los efectos aleatorios, o las altas correlaciones entre los terminos aleatorios, 
##   apuntan a problemas en la especificacion del modelo, con lo que no son algo que el analisis de permutacion pueda resolver.
##   Esos problemaso requieren nuestras decisiones para evaluar la especificacion del modelo, reducir la complejidad de la estructura de efectos aleatorios, 
##   o indagar en los datos para posibles problemas (e.g., insuficiente variabilidad y pocos datos dentro de los niveles de un efecto aleatorio).
## En resumen, el analisis de permutacion es una tecnica de remuestreo no parametrico que es robusta a muchas violaciones 
##   de las suposiciones relacionadas con la distribucion de los datos, como son los valores atipicos o las observaciones
##   con grandes valores de leverage, heterocedasticidad y normalidad residual. 
##   Funciona creando una distribucion de los estadisticos de interes derivados de nuestro modelo bajo la asuncion de que la hipotesis nula es cierta,
##   a partir de los datos realmente observados, por lo que no depende de las suposiciones de las pruebas parametricas.
##   Sin embargo, el analisis permutacional no aborda directamente el problema de la multicolinealidad. 
##   Tampoco resolvera los problemas vinculados con la singularidad, baja varianza en los efectos aleatorios y gran correlacion entre ellos.
## El numero de permutaciones necesarias para un analisis permutacional esta relacionado con la precision que necesitamos para el valor de P. 
##   Una regla general es usar un numero de permutaciones que sea al menos 10 veces el reciproco de la precision deseada para el valor de P.
##   Por ejemplo, para una precision de valor P de 0.01, tendriamos al menos 1 / 0.01 = 100 permutaciones * 10 = 1000 permutaciones.
##   Con fines practicos muchos estudios utilizan alrededor de 1000 a 5000 permutaciones. Sin embargo, con mas permutaciones podemos proporcionar valores P mas precisos.
## Uno de los desafios con los modelos mixtos es la estimacion de los grados de libertad del termino error para los efectos fijos. 
##   Diferentes metodos (como los utilizados por nlme y lme4, o los metodos de Satterthwaite o Kenward-Roger de lme4) 
##     pueden producir diferentes resultados, especialmente con tamannios de muestra pequenos, lo que a su vez puede llevar 
##     a diferencias en los valores de significacion P de los efectos fijos.
##   Los analisis basados en la permulatcion de la variable respuesta proporcionan una solucion a este problema. 
##     Debido a que crean una distribucion de los estadisticos de interes a partir de los datos originales bajo la asuncion de la hipotesis nula (Ho),
##     no requieren la estimacion de grados de libertad de la misma manera que las pruebas parametricas tradicionales. 
##     Por tanto, los analisis permutacionales pueden proporcionar la significacion de los efectos fijos en modelos mixtos 
##     sin tener que preocuparse por el calculo de los grados de libertad.

##   *** MUY LENTO ***
library(bayestestR)

## MODELOS lmer
N_boots <- 1000
N_coefs <- length(summary(modelo)$coefficients[,1])
N_Fs <- dim(anova(modelo, ddf="Satterthwaite", type=3))[1]
null_coeficientes <- as.data.frame(matrix(9999, nrow=N_boots, ncol=N_coefs))
null_Fs <- as.data.frame(matrix(9999, nrow=N_boots, ncol=N_Fs))
null_R2 <- as.data.frame(matrix(9999, nrow=N_boots, ncol=2))
respuesta <- get.response(modelo)
## tambien podemos usar ddf="Kenward-Roger", aunque es mas lento que con ddf="Satterthwaite"
for (i in 1:N_boots) {
  print(i)
  datos$respuesta.null <- sample(respuesta, replace=FALSE)
  lmer_null <- update(modelo, respuesta.null~., REML=TRUE)
  null_coeficientes[i,] <- summary(lmer_null, ddf="Satterthwaite")$coefficients[,1]
  null_Fs[i,] <- anova(lmer_null, ddf="Satterthwaite", type=3)[,5]
  null_R2[i,] <- r.squaredGLMM(lmer_null)
}
## comprobad que en la matriz de null_coeficientes no hay valores 9999 indicativos de errores por problemas de estima
colnames(null_coeficientes) <- rownames(summary(lmer_null)$coefficients)
colnames(null_Fs) <- rownames(anova(modelo, ddf="Satterthwaite", type=3))
colnames(null_R2) <- c("R2m", "R2c")
tabla.null <- describe(null_coeficientes, quant=c(0.025, 0.975, 0.005, 0.995))[, c(-1, -5:-10, -13)]
tabla.null$Estimate <- summary(modelo)$coefficients[,1]
colnames(tabla.null)[3] <- "std_error"
tabla.null$Z <- (tabla.null$Estimate - tabla.null$mean) / tabla.null$std_err
tabla.null$P <- round(pt(abs(tabla.null$Z), df=999999, lower.tail=F), 3)
F.null <- describe(null_Fs, quant=c(0.90, 0.95, 0.99, 0.995))[,c(2, 14:17)]
F.null$Estimate <- anova(modelo, ddf="Satterthwaite", type=3)[,5]
for (i in 1:N_Fs) {F.null$P[i] <- 1-ecdf(null_Fs[,i])(F.null$Estimate[i])}
R2.null <- describe(null_R2, quant=c(0.90, 0.95, 0.99, 0.995))[,c(2, 14:17)]
R2.null$Estimate <- r.squaredGLMM(modelo)[1,]
for (i in 1:2) {R2.null$P[i] <- 1-ecdf(null_R2[,i])(R2.null$Estimate[i])}

## COEFICIENTES
## que el Estimate del modelo mixto de interes quede fuera de los intervalos al 95% (Q0.025 - Q0.975) o 99% (Q0.005 - Q0.995)
print(tabla.null, digits=3)
plot(ecdf(null_coeficientes$t_sonda)); abline(h=c(0.025, 0.975), col="red")     ## *** PONED AQUI VUESTROS DATOS ***; para sacar una p para un coeficiente observado en el modelo mixto de interes

## que las F y R2 queden a la derecha de los cuantiles correspondientes a alpha = 1 - quantile
##  TABLA DE ANOVA
print(F.null, digits=3)
1-ecdf(null_Fs$especie)(2.290)  ## *** PONED AQUI VUESTROS DATOS ***; para sacar una p para una F observada en el modelo mixto de interes

##  R2 DEL MODELO
print(R2.null, digits=3)
1-ecdf(null_R2$R2m)(0.475)      ## *** PONED AQUI VUESTROS DATOS ***;  para sacar una p para una R2 observada en el modelo mixto de interes
plot(density(null_R2$R2c), lwd=2, col="red", xlab="R2 total del analisis permulacional", main="", xlim=c(0, 1)); abline(v=R2.null[2,6], col="blue")
plot(density(null_R2$R2m), lwd=2, col="red", xlab="R2 efectos fijos del analisis permulacional", main="", xlim=c(0, 1)); abline(v=R2.null[1,6], col="blue")

## INTERVALOS con los metodos ETI (Equal-Tailed Interval; por defecto) y HDI (Highest Density Interval)
##  usando el metodo ETI
plot(eti(null_R2, ci=0.95))
print(eti(null_coeficientes, ci=0.95), digits=3)  ## igual que la salida previa
plot(eti(null_coeficientes, ci=0.95))
##  usando el metodo HDI
plot(hdi(null_R2, ci=0.95))
print(hdi(null_coeficientes, ci=0.95), digits=3)
plot(hdi(null_coeficientes, ci=0.95))


## MODELOS lme
N_boots <- 1000
N_coefs <-  length(summary(modlme)$coefficients$fixed)
N_Fs <- dim(anova.lme(modlme, type="marginal"))[1]
null_coeficientes_lme <- as.data.frame(matrix(9999, nrow=N_boots, ncol=N_coefs))
null_Fs_lme <- as.data.frame(matrix(9999, nrow=N_boots, ncol=N_Fs))
null_R2_lme <- as.data.frame(matrix(9999, nrow=N_boots, ncol=2))
respuesta <- get.response(modlme)
for (i in 1:N_boots) {
  print(i)
  datos$respuesta.null <- sample(respuesta, replace=FALSE)
  lme_null <- update(modlme, respuesta.null~., method="REML")
  null_coeficientes_lme[i,] <- summary(lme_null)$coefficients$fixed
  null_Fs_lme[i,] <- anova.lme(lme_null, type="marginal")[,3]
  null_R2_lme[i,] <- r.squaredGLMM(lme_null)
}
## comprobad que en la matriz de null_coeficientes no hay valores 9999 indicativos de errores por problemas de estima
colnames(null_coeficientes_lme) <- names(summary(lme_null)$coefficients$fixed)
colnames(null_Fs_lme) <- rownames(anova.lme(modlme, type="marginal"))
colnames(null_R2_lme) <- c("R2m", "R2c")
tabla.null.lme <- describe(null_coeficientes_lme, quant=c(0.025, 0.975, 0.005, 0.995))[, c(-1, -5:-10, -13)]
tabla.null.lme$Estimate <- summary(modlme)$coefficients$fixed
colnames(tabla.null.lme)[3] <- "std_error"
F.null_lme <- describe(null_Fs_lme, quant=c(0.90, 0.95, 0.99, 0.995))[,c(2, 14:17)]
F.null_lme$Estimate <- anova.lme(modlme, type="marginal")[,3]
for (i in 1:N_coefs) {F.null_lme$P[i] <- 1-ecdf(null_Fs_lme[,i])(null_Fs_lme$Estimate[i])}
R2.null_lme <- describe(null_R2_lme, quant=c(0.90, 0.95, 0.99, 0.995))[,c(2, 14:17)]
R2.null_lme$Estimate <- r.squaredGLMM(modlme)[1,]
for (i in 1:2) {null_R2_lme$P[i] <- 1-ecdf(null_R2_lme[,i])(null_R2_lme$Estimate[i])}
## que el Estimate del modelo mixto de interes quede fuera de los intervalos al 95% (Q0.025 - Q0.975) o 99% (Q0.005 - Q0.995)
print(tabla.null.lme, digits=3)
## *** PONED AQUI VUESTROS DATOS *** eligiendo una variable predictora
plot(ecdf(null_coeficientes_lme$t_sonda)); abline(h=c(0.025, 0.975), col="red")
## que las F y R2 queden a la derecha de los cuantiles correspondientes a alpha = 1 - quantile
print(F.null_lme[-1,], digits=3)
1-ecdf(null_Fs_lme$especie)(2.290)  ## *** PONED AQUI VUESTROS DATOS ***; para sacar una p para una F observada en el modelo mixto de interes
print(R2.null_lme, digits=3)
1-ecdf(null_R2_lme$R2m)(0.475)      ## *** PONED AQUI VUESTROS DATOS ***;  para sacar una p para una R2 observada en el modelo mixto de interes
plot(density(null_R2_lme$R2c), lwd=2, col="red", xlab="R2 total del analisis permulacional", main="", xlim=c(0, 1)); abline(v=R2.null_lme[2,6], col="blue")
plot(density(null_R2_lme$R2m), lwd=2, col="red", xlab="R2 efectos fijos del analisis permulacional", main="", xlim=c(0, 1)); abline(v=R2.null_lme[1,6], col="blue")
##
## INTERVALOS con los metodos ETI (Equal-Tailed Interval; por defecto) y HDI (Highest Density Interval)
##  que los valores bajo "Estimate" queden fuera de los intervalos CI_low - CI_high
##  usando el metodo ETI
tabla.eti <- as.data.frame(eti(null_coeficientes_lme, ci=0.95))
tabla.eti$Estimate <- tabla.null.lme$Estimate
tabla.eti
plot(eti(null_coeficientes_lme, ci=0.95))
##  usando el metodo HDI
tabla.hdi <- as.data.frame(hdi(null_coeficientes_lme, ci=0.95))
tabla.hdi$Estimate <- tabla.null.lme$Estimate
tabla.hdi
plot(hdi(null_coeficientes_lme, ci=0.95))





########################################
#### Estima ROBUSTA del modelo lmer ####
########################################

library(robustlmm)
## https://cran.r-project.org/web/packages/robustlmm/vignettes/rlmer.pdf
## metodos:
## DAStau: (default) For this method, the consistency factors are computed using numerical
##         quadrature. This is slower but yields more accurate results. This is the direct analogue to
##         the DAS-estimate in robust linear regression.
## DASvar: This method computes the consistency factors using a direct approximation
##         which is faster but less accurate. For complex models with correlated random effects
##         with more than one correlation term, this is the only method available.
eqt <- formula(modelo)
## modelos rlmer
##   para modelos mixtos complejos con terminos aleatorios correlacionados el metodo DASvar es el unico posible
mm.rvar2 <- rlmer(eqt, data=datos, rel.tol=1e-08, max.iter=1000, method="DASvar")
mm.rtau2 <- rlmer(eqt, data=datos, rel.tol=1e-08, max.iter=1000, method="DAStau")
mm.rvar <- update(mm.rvar2, rho.sigma.e=psi2propII(smoothPsi, k=2.28), rho.sigma.b=psi2propII(smoothPsi, k=2.28))
mm.rtau <- update(mm.rtau2, rho.sigma.e=psi2propII(smoothPsi, k=2.28), rho.sigma.b=psi2propII(smoothPsi, k=2.28))

## los PESOS DE LAS UNIDADES MUESTRALES en los modelos se obtienen para las observaciones ("w_e") y para los efectos aleatorios ("w_b")
getME(mm.rvar, "w_e")
getME(mm.rvar, "w_b")
plot(getME(mm.rvar, "w_e"), getME(mm.rtau, "w_e"))
plot(getME(mm.rtau, "w_e"), getME(mm.rtau2, "w_e"))
par(mfcol=c(1,1))    ## fija un solo panel grafico
plot(residuals(modelo, type="response"), getME(mm.rvar, "w_e"))
plot(scale(residuals(modelo, type="response")), getME(mm.rvar, "w_e"))
plot(cooks.distance(modelo), getME(mm.rvar, "w_e"))
plot(hatvalues(modelo), getME(mm.rvar, "w_e"))

## RESIDUOS del modelo
mcp.fnc(mm.rvar)
mcp.fnc(mm.rtau)
mcp.fnc(modelo)
## plots de los residuos y pesos del modelo
##   demanda en la consola: Hit <Return> to see next plot
plot(mm.rvar, which=c(1:3))
plot(mm.rtau, which=c(1:3))

## RESULTADOS del modelo
##   varianza de los terminos aleatorios
summary(mm.rvar)$varcor
summary(mm.rtau)$varcor
summary(modelo)$varcor
##
##   los modelos rlmer no proporcionan significaciones
summary(mm.rvar)$coefficients
summary(mm.rtau)$coefficients
summary(modelo, ddf="lme4")$coefficients

## SIGNIFICACIONES de los coeficientes del modelo
## como rlmer no proporciona significaciones, 
## podriamos utilizar los grados de libertad proporcionados por el modelo lmer de interes
##   pero esta aproximacion no nos brinda obtener tablas de anova, tamannios de efectos, visualizacion de efectos parciales
##
##    modelo robusto rlmer DASvar con significaciones asociadas a los grados de libertad de Kenward-Roger
tabla.mm.rvar <- as.data.frame(summary(mm.rvar)$coefficients)
tabla.mm.rvar$dfKR <- round(summary(modelo, ddf="Kenward-Roger")$coefficients[,3], 1)
tabla.mm.rvar$P_dfKR <- (1-pt(abs(tabla.mm.rvar[,3]), df=tabla.mm.rvar[,4]))
round(tabla.mm.rvar, 4)
##
##    modelo robusto rlmer DAStau con significaciones asociadas a los grados de libertad de Kenward-Roger
tabla.mm.rtau <- as.data.frame(summary(mm.rtau)$coefficients)
tabla.mm.rtau$dfKR <- round(summary(modelo, ddf="Kenward-Roger")$coefficients[,3], 1)
tabla.mm.rtau$P_dfKR <- (1-pt(abs(tabla.mm.rtau[,3]), df=tabla.mm.rtau[,4]))
round(tabla.mm.rtau, 4)
##
##
## tambien podriamos utilizar los pesos de las observaciones e incluirlos en el modelo lmer original
modelo.w.var <- update(modelo, weights=getME(mm.rvar, "w_e"))
round(summary(modelo.w.var, ddf="Kenward-Roger")$coefficients,4)
Anova(modelo.w.var, type=3, test="F")
epsilon_squared(anova(modelo.w.var, ddf="Kenward-Roger", type=3))  
##
modelo.w.tau <- update(modelo, weights=getME(mm.rtau, "w_e"))
round(summary(modelo.w.tau, ddf="Kenward-Roger")$coefficients,4)
Anova(modelo.w.tau, type=3, test="F")
epsilon_squared(anova(modelo.w.tau, ddf="Kenward-Roger", type=3))  
##
##
## modelo original lmer sin pesar las observaciones:
round(summary(modelo, ddf="Kenward-Roger")$coefficients, 4)
Anova(modelo, type=3, test="F")
epsilon_squared(anova(modelo, ddf="Kenward-Roger", type=3))  







########################################
#### OTRO JUEGO DE DATOS Y EJEMPLOS ####
########################################

## datos de internet
meconecto <- url("http://www.lmcarrascal.eu/cursos/nwayMIXED_ANOVA_1.RData")  ## abro una conexion con internet
load(meconecto)                                                               ## cargo el archivo de la conexion
close(meconecto); rm(meconecto)                                               ## cierro la conexion

options(contrasts=c(factor="contr.sum", ordered="contr.poly"))

control.lme <- lmeControl(maxIter=200, msMaxIter=200, msVerbose=TRUE, opt="optim")    ## el optimizador  opt="nlminb"  hay veces que da problemas
control.lmer <- lmerControl(check.conv.grad=.makeCC(action ="ignore", tol=1e-6, relTol=NULL), optimizer="bobyqa", optCtrl=list(maxfun=100000))

## Intraclass correlation: the "unconditional means model" (or "null model"); Random intercept with fixed mean
##   diferencias entre los niveles del factor aleatorio (individuo)
##   random=~1|individuo,   random=list(individuo=~1)   y random=list(individuo=pdSymm(~1))    son equivalentes
lmer.0 <- lmer(ln_uso10h ~ 1 + (1|individuo), data=datos, control=control.lmer, REML=TRUE)
##
lme.0 <- lme(ln_uso10h ~ 1, random=list(individuo=pdSymm(~1)), data=datos, control=control.lme, method="REML")

## modelo de RANDOM INTERCEPT & FIXED SLOPES
lmer.1 <- lmer(ln_uso10h ~ temperatura+spp*distancia + (1|individuo), data=datos, control=control.lmer, REML=TRUE)
##
lme.1.Symm <- lme(ln_uso10h ~ temperatura+spp*distancia, random=list(individuo=pdSymm(~1)), data=datos, control=control.lme, method="REML")
lme.1.Chol <- lme(ln_uso10h ~ temperatura+spp*distancia, random=~1|individuo, data=datos, control=control.lme, method="REML")    ## este modelo y los dos siguientes son equivalentes
lme.1.Chol <- lme(ln_uso10h ~ temperatura+spp*distancia, random=list(individuo=~1), data=datos, control=control.lme, method="REML")
lme.1.Chol <- lme(ln_uso10h ~ temperatura+spp*distancia, random=list(individuo=pdLogChol(~1)), data=datos, control=control.lme, method="REML")
lme.1.Ident <- lme(ln_uso10h ~ temperatura+spp*distancia, random=list(individuo=pdIdent(~1)), data=datos, control=control.lme, method="REML")

## modelo de CORRELATED RANDOM INTERCEPT & RANDOM SLOPES, with CORRELATED random slopes for different covariates
lmer.2 <- lmer(ln_uso10h ~ temperatura+spp*distancia + (1+temperatura+distancia|individuo), data=datos, control=control.lmer, REML=TRUE)
##
lme.2.Symm <- lme(ln_uso10h ~ temperatura+spp*distancia, random=list(individuo=pdSymm(~1+temperatura+distancia)), data=datos, control=control.lme, method="REML")
lme.2.Chol <- lme(ln_uso10h ~ temperatura+spp*distancia, random=~(1+temperatura+distancia)|individuo, data=datos, control=control.lme, method="REML")    ## este modelo y los dos siguientes son equivalentes (no lo comento mas veces)
lme.2.Chol <- lme(ln_uso10h ~ temperatura+spp*distancia, random=list(individuo=~1+temperatura+distancia), data=datos, control=control.lme, method="REML")   
lme.2.Chol <- lme(ln_uso10h ~ temperatura+spp*distancia, random=list(individuo=pdLogChol(~1+temperatura+distancia)), data=datos, control=control.lme, method="REML")
lme.2.SymComp <- lme(ln_uso10h ~ temperatura+spp*distancia, random=list(individuo=pdCompSymm(~1+temperatura+distancia)), data=datos, control=control.lme, method="REML")
## estos modelos que siguen no tienen sentido al considerar independencia entre terminos aleatorios
lme.2.Diag <- lme(ln_uso10h ~ temperatura+spp*distancia, random=list(individuo=pdDiag(~1+temperatura+distancia)), data=datos, control=control.lme, method="REML")
lme.2.Ident <- lme(ln_uso10h ~ temperatura+spp*distancia, random=list(individuo=pdIdent(~1+temperatura+distancia)), data=datos, control=control.lme, method="REML")

## modelo de UNCORRELATED RANDOM INTERCEPT & RANDOM SLOPES for different covariates
lmer.3 <- lmer(ln_uso10h ~ temperatura+spp*distancia + (1|individuo)+(temperatura|individuo)+(distancia|individuo), data=datos, control=control.lmer, REML=TRUE)
##
lme.3.Symm <- lme(ln_uso10h ~ temperatura+spp*distancia, random=list(individuo=pdSymm(~1), individuo=pdSymm(~temperatura), individuo=pdSymm(~distancia)), data=datos, control=control.lme, method="REML")
lme.3.Chol <- lme(ln_uso10h ~ temperatura+spp*distancia, random=list(individuo=pdLogChol(~1), individuo=pdLogChol(~temperatura), individuo=pdLogChol(~distancia)), data=datos, control=control.lme, method="REML")

## modelo de UNCORRELATED RANDOM INTERCEPT & RANDOM SLOPES; elimina la covarianza interceptos-pendientes
lmer.4 <- lmer(ln_uso10h ~ temperatura+spp*distancia + (1|individuo)+(temperatura-1|individuo)+(distancia-1|individuo), data=datos, control=control.lmer, REML=TRUE)
##
lme.4.Symm <- lme(ln_uso10h ~ temperatura+spp*distancia, random=list(individuo=pdSymm(~1), individuo=pdSymm(~temperatura-1), individuo=pdSymm(~distancia-1)), data=datos, control=control.lme, method="REML")
lme.4.Chol <- lme(ln_uso10h ~ temperatura+spp*distancia, random=list(individuo=pdLogChol(~1), individuo=pdLogChol(~temperatura-1), individuo=pdLogChol(~distancia-1)), data=datos, control=control.lme, method="REML")
## estos modelos que siguen consideran independencia entre terminos aleatorios
lme.4.Diag <- lme(ln_uso10h ~ temperatura+spp*distancia, random=list(individuo=pdDiag(~1+temperatura+distancia)), data=datos, control=control.lme, method="REML")
lme.4.Ident <- lme(ln_uso10h ~ temperatura+spp*distancia, random=list(individuo=pdIdent(~1), individuo=pdIdent(~temperatura-1), individuo=pdIdent(~distancia-1)), data=datos, control=control.lme, method="REML")

## modelo de FIXED INTERCEPT & RANDOM SLOPES, with UNCORRELATED random slopes for different covariates
lmer.5 <- lmer(ln_uso10h~temperatura+spp*distancia + (temperatura-1|individuo)+(distancia-1|individuo), data=datos, control=control.lmer, REML=TRUE)
##
lme.5.Symm <- lme(ln_uso10h ~ temperatura+spp*distancia, random=list(individuo=pdSymm(~temperatura-1), individuo=pdSymm(~distancia-1)), data=datos, control=control.lme, method="REML")
lme.5.Chol <- lme(ln_uso10h ~ temperatura+spp*distancia, random=list(individuo=pdLogChol(~temperatura-1), individuo=pdLogChol(~distancia-1)), data=datos, control=control.lme, method="REML")
## estos modelos que siguen consideran independencia entre terminos aleatorios
lme.5.Diag <- lme(ln_uso10h ~ temperatura+spp*distancia, random=list(individuo=pdDiag(~temperatura+distancia-1)), data=datos, control=control.lme, method="REML")
lme.5.Ident <- lme(ln_uso10h ~ temperatura+spp*distancia, random=list(individuo=pdIdent(~temperatura-1), individuo=pdIdent(~distancia-1)), data=datos, control=control.lme, method="REML")

## modelo de CORRELATED RANDOM INTERCEPT & RANDOM SLOPES, with UNCORRELATED random slopes for different covariates
lmer.6 <- lmer(ln_uso10h ~ temperatura+spp*distancia + (1+temperatura|individuo)+(1+distancia|individuo), data=datos, control=control.lmer, REML=TRUE)
##
lme.6.Symm <- lme(ln_uso10h ~ temperatura+spp*distancia, random=list(individuo=pdSymm(~1+temperatura), individuo=pdSymm(~1+distancia)), data=datos, control=control.lme, method="REML")
lme.6.Chol <- lme(ln_uso10h ~ temperatura+spp*distancia, random=list(individuo=pdLogChol(~1+temperatura), individuo=pdLogChol(~1+distancia)), data=datos, control=control.lme, method="REML")


## COMPARACION DE MODELOS: CORRELATED RANDOM INTERCEPTS AND SLOPES MIXED MODELS
##   valores de Akaike AICc y df de los efectos del modelo
AICc(lmer.1, lme.1.Symm, lme.1.Chol, lme.1.Ident)
AICc(lmer.2, lme.2.Symm, lme.2.Chol, lme.2.SymComp, lme.2.Diag, lme.2.Ident)
AICc(lmer.3, lme.3.Symm, lme.3.Chol)
AICc(lmer.4, lme.4.Symm, lme.4.Chol, lme.4.Diag, lme.4.Ident)
AICc(lmer.5, lme.5.Symm, lme.5.Chol, lme.5.Diag, lme.5.Ident)
AICc(lmer.6, lme.6.Symm, lme.6.Chol)


##  matrices de varianzas-correlaciones
summary(lmer.2)$varcor
VarCorr(lme.2.Symm)
VarCorr(lme.2.Chol)
VarCorr(lme.2.SymComp)
VarCorr(lme.2.Diag)
VarCorr(lme.2.Ident)
##
summary(lmer.4)$varcor
VarCorr(lme.4.Symm)
VarCorr(lme.4.Chol)
VarCorr(lme.4.Diag)
VarCorr(lme.4.Ident)


## RESULTADOS
##  seleccionamos un modelo lme4::lmer y otro nlme::lme
modelo <- lmer.2
modlme <- lme.2.Symm

plot(simulateResiduals(modelo))
plot(modlme, id=0.05)

check_collinearity(modelo)
check_collinearity(modlme)

round(summary(modelo, ddf="Kenward-Roger")$coefficients, 5)
round(summary(modlme)$tTable, 5)
round(summary(modelo, ddf="Satterthwaite")$coefficients, 5)

anova(modelo, ddf="Kenward-Roger", type=3)
anova.lme(modlme, type="marginal")[-1,]
anova(modelo, ddf="Satterthwaite", type=3)

r.squaredGLMM(modelo)
r.squaredGLMM(modlme)

epsilon_squared(anova(modelo, ddf="Kenward-Roger", type=3))  
epsilon_squared(anova.lme(modlme, type="marginal"))
epsilon_squared(anova(modelo, ddf="Satterthwaite", type=3))  

rand(modelo)
anova(gls(formula(modlme), data=datos, method="REML"), modlme)

plot(allEffects(modelo), rotx=45, rug=F)
plot(allEffects(modlme), rotx=45, rug=F)







######################################################
####         TRANSFORMACION DE BOX-COX            ####
## busqueda automatizada del valor optimo de lambda ##
######################################################

## TRANSFORMACION BOX-COX
##   funcion boxCox del paquete {car}
## anniadir +1 a la variable respuesta si tiene ceros, en cuyo caso creamos una ecuacion de modelo nueva
## la transformacion es: Y' = (Y^lambda - 1)/lambda      ,o,     Y' = ((Y+1)^lambda - 1)/lambda 
## la siguiente linea de codigo SOLO TRANSFORMA la respuesta de nuestra ecuacion eqt
## podemos afinar y restringir aun mas el rango de lambda; por ejemplo: lambda=seq(0,1, 1/1000))

## corremos un modelo general lineal frecuentista (lm)
formula(modelo)
## eliminamos los terminos aleatorios
##   *** PONED AQUI VUESTROS DATOS ***
eqt <- as.formula(I(uso10h+1) ~ temperatura + spp * distancia)
mod_bc <- lm(eqt, data=datos)

## comenzamos el proceso automatizado de busqueda de la lamba para la transformacion optima
##   si no observamos un maximo, modificamos el rango de lambdas en lambda=seq(-2,2, 1/100)
print(boxCox(mod_bc, lambda=seq(-2,2, 1/100)))
lambda.bc <- with(boxCox(mod_bc, lambda=seq(-2,2, 1/1000)), x[which.max(y)])
print(c("parametro lambda de la transformacion Box-Cox =",round(lambda.bc, 3)), quote=FALSE)

## transformamos la variable a continuacion ... que llamamos respuesta.bc incluida en datos
datos$respuesta.bc <- ((mod_bc$model[,1]^lambda.bc)-1)/(lambda.bc)
## para explorar la relacion entre la respuesta original y su transformacion Box-Cox
plot(datos$respuesta.bc, mod_bc$model[,1], xlab="variable respuesta transformada", ylab="variable respuesta original")

## renombramos y guardamos los modelos previos
modelo.previo <- modelo
modlme.previo <- modlme

## ahora rehacemos el modelo 
##   sustituyendo la variable respuesta del modelo por datos$respuesta.bc
modelo <- update(modelo, respuesta.bc~., na.action=na.omit)
modlme <- update(modlme, respuesta.bc~., na.action=na.omit)






###################################################################
#### COMENTARIOS ACERCA DE ESTRUCTURAS ALEATORIAS EN nlme::lme ####
###################################################################

## Descripcion de una estructura de efectos aleatorios de modelos mixtos lme4::lmer
##  ejemplo: (1+temperatura+distancia|individuo)
##  La estructura por defecto para la matriz de varianza-covarianza de efectos aleatorios asume lo siguiente:
##  Definicion positiva (Positive-definiteness): todas las varianzas (elementos diagonales de la matriz) son positivas.
##  Simetria (Symmetry): La matriz de varianza-covarianza es simetrica 
##    (la covarianza entre dos efectos aleatorios cualesquiera es la misma independientemente del orden en que se consideren
##    la covarianza del efecto aleatorio de `x` sobre `y` es la misma que la de `y` sobre `x`
##  Heterogeneidad (Heterogeneity): se asume que las varianzas de los diferentes efectos aleatorios pueden ser diferentes. 
##    Por tanto, los elementos diagonales de la matriz de covarianza de los efectos aleatorios (varianzas), no estan obligados a ser iguales.
##  Estructura de covarianzas: se estiman por defecto todas las covarianzas entre pares de efectos aleatorios, 
##    permitiendo que sean diferentes (i.e., no estan obligados a ser cero o iguales entre si).
##  La estructura de los efectos aleatorios en lme4::lmer es bastante flexible, permitiendo la heterogeneidad en las varianzas 
##    y covarianzas de los efectos aleatorios. 
##    No hace suposiciones restrictivas sobre la estructura de la matriz de varianza-covarianza de los efectos aleatorios, 
##    a diferencia de lo que ocurre algunas de las estructuras que podemos especificar en nlme::lme 
##    (e.g., con pdDiag, pdCompSymm o pdIdent).

## Si solo abordamos un modelo de intercepto aleatorio con nlme::lme (i.e., random = list(id = ~1)), 
##   entonces las estructuras de covarianza mas complejas no seran necesarias (e.g., pdDiag, pdSymm, pdCompSymm, pdIdent, etc).

## Descripcion de una estructura de efectos aleatorios de modelos mixtos nlme::lme mediante pdSymm
##  ejemplo: random=list(individuo=pdSymm(~1+temperatura+distancia))
##  La funcion pdSymm, abreviatura de "positive-definite symmetric", especifica un tipo particular de estructura de efectos aleatorios. 
##  La estructura pdSymm proporciona bastante flexibilidad para modelizar la estructura de efectos aleatorios en los datos, 
##   permitiendo que cada efecto aleatorio tenga su propia varianza y que cada par de efectos aleatorios tenga su propia covarianza.
##  Permite una estructura de efectos aleatorios multivariante en la que la matriz de covarianza de los efectos aleatorios es simetrica 
##   y definida positivamente. Mas en detalle:
##  Definicion positiva (Positive-definiteness): Todas las varianzas (que son los elementos diagonales de la matriz de covarianza) son positivas. 
##  Simetria (Symmetry): La matriz de covarianza es simetrica, lo que significa que la covarianza entre dos efectos aleatorios diferentes 
##   es la misma independientemente del orden (e.g., la covarianza del efecto A con el efecto B es la misma que la covarianza del efecto B con el efecto A).
##  Heterogeneidad (Heterogeneity): Diferentes efectos aleatorios pueden tener diferentes varianzas, es decir, no se requiere que sean homogeneas.
##  Las covarianzas pueden ser diferentes: La covarianza entre un par de efectos aleatorios puede ser diferente de la covarianza entre otro par de efectos aleatorios. 
##   Solo tienen que ser simetricas con respecto a cada par, pero no tienen por que ser iguales en los distintos pares de efectos aleatorios.
##  En el ejemplo se permite que el intercepto aleatorio y las pendientes aleatorias esten correlacionadas.

## Descripcion de una estructura de efectos aleatorios de modelos mixtos nlme::lme mediante la parametrizacion "log-Cholesky".
##  ejemplo: random=list(individuo=~1+temperatura+distancia)
##  Esta estructura, de escritura sencilla, es equivalente a:  random=list(individuo=pdLogChol(~1+temperatura+distancia))
##  La descomposicion de Cholesky es un metodo utilizado para la descomposicion de matrices de covarianza hermitianas definidas positivamente. 
##   La descomposicion de Cholesky de una matriz de covarianza es una matriz triangular inferior.
##  Las varianzas (elementos diagonales) se transforman mediante un logaritmo, lo que garantiza que sigan siendo positivas.
##  Los elementos no diagonales se estiman libremente, sin ninguna restriccion de simetria.
##  La principal diferencia entre pdLogChol y pdSymm radica en la forma en que tratan los elementos no diagonales (covarianzas). 
##   La parametrizacion log-Cholesky permite una estructura de covarianza sin restricciones, mientras que pdSymm asume explicitamente una estructura simetrica. 
##   Sin embargo, ambos metodos garantizan que las varianzas sigan siendo positivas.
##  random=list(id = ~1 + x + z) permite interceptos y pendientes aleatorios pero no asume explicitamente una estructura particular de covarianzas.
##  Por tanto, especifica una estructura de correlacion no estructurada, con lo que no se asume que haya correlacion entre dos efectos aleatorios. 
##   Puede ser una buena opcion si se desconoce como de probable es que se correlacionen los efectos aleatorios.
##  La estructura de correlacion simetrica pdSymm dara lugar normalmente a estimaciones mas bajas de los componentes de la varianza
##    que la estructura de correlacion no estructurada con parametrizacion "log-Cholesky".
##  En resumen, la estructura logCholesky supone que los elementos diagonales de la matriz triangular inferior son 
##    varianzas transformadas logaritmicamente de los efectos aleatorios, y los elementos no diagonales de la matriz triangular inferior son 
##    correlaciones entre los efectos aleatorios.
##    Esta estructura logCholesky es mas flexible que las estructuras diagonales tradicionales (pdDiag) o de simetria compuesta (pdCompSymm) 
##    porque permite diferentes varianzas y correlaciones entre los efectos aleatorios.
##    Lo que puede ser util si sospechamos que los efectos aleatorios no tienen varianza homogenea y/o no estan perfectamente correlacionados entre si.
##    (importante cuando existe la posibilidad de que los efectos aleatorios pueden tener una estructura de correlacion compleja)

## Descripcion de una estructura de efectos aleatorios de modelos mixtos nlme::lme mediante pdDiag
##  ejemplo: random=list(individuo=pdDiag(~1+temperatura+distancia))
##  La funcion pdDiag especifica una estructura de covarianza diagonal para los efectos aleatorios (las celdas extra-diagonales valen cero-0)
##  En el ejemplo estamos especificando un modelo con interceptos aleatorios (~1) y pendientes aleatorias para las variables x y z (x + z), 
##   cada uno de los cuales es un efecto aleatorio independiente. Asumimos que esos efectos aleatorios varian en los distintos niveles del factor aleatorio id.
##  pdDiag especifica que los efectos aleatorios tienen sus propias varianzas, pero no estan correlacionados entre si. 
##   Por tanto, la matriz de covarianza de los efectos aleatorios es una matriz diagonal, con todos los elementos no diagonales (las covarianzas) iguales a cero. 
##   Esto significa que cualquier aumento o disminucion en un efecto aleatorio (como la pendiente aleatoria para x) no proporciona ninguna informacion 
##    sobre aumentos o disminuciones probables en cualquier otro efecto aleatorio (como el intercepto aleatorio o la pendiente aleatoria para z).
##  La caracteristica propia de la estructura pdDiag es que se supone que las covarianzas entre cualquier par de efectos aleatorios es cero, 
##   asumiendo que los efectos aleatorios no estan correlacionados entre si. 
##  Permite varianzas diferentes para cada efecto aleatorio asumiendo que las covarianzas entre los efectos aleatorios son = cero;
##      es decir, los efectos aleatorios no estan correlacionados.
##  Se trata de una estructura mas restrictiva en comparacion con pdSymm o la parametrizacion log-Cholesky.

## Descripcion de una estructura de efectos aleatorios de modelos mixtos nlme::lme mediante pdIdent
##  ejemplo: random=list(id=pdIdent(~1 + x + z))
## pdIdent del paquete nlme se utiliza para especificar una estructura de efectos aleatorios en la que 
##   las varianzas de los efectos aleatorios son todas iguales (homocedasticidad) y todas las covarianzas son cero 
##   (i.e., los efectos aleatorios de los interceptos (~1) y pendientes para las variables x y z no estan correlacionados).
## La matriz de varianzas-covarianzas de los efectos aleatorios es una matriz de identidad, 
##   en la que todos los elementos diagonales (varianzas) son iguales y positivos, 
##   y todos los elementos no diagonales (covarianzas) son cero.
## Por tanto, es una forma de modelizar interceptos aleatorios y pendientes aleatorias que varian a traves 
##   de diferentes niveles de un factor aleatorio, pero asumiendo correlacion cero entre los diferentes efectos aleatorios, 
##   y que todos los efectos aleatorios tienen la misma varianza. 
## Se trata de una estructura aun mas restrictiva que pdDiag.

## nlme::lme permite varias estructuras mas para la matriz de covarianza de efectos aleatorios.
##  pdCompSymm: Esta funcion especifica una estructura de "simetria compuesta". 
##    En esta estructura, todas las varianzas son iguales (homocedasticidad) y todas las covarianzas tambien son iguales.
##    Puede ser apropiada en diversos escenarios analiticos y experimentales:
##      Medidas repetidas: Si tenemos medidas repetidas en los mismos sujetos y esperamos que la correlacion entre las medidas 
##        sea constante a lo largo del tiempo. Sin embargo, tambien supone que la varianza de los efectos aleatorios es constante  
##        en el tiempo, lo que puede no ser cierto en todos los casos.
##      Datos transversales: Si tenemos datos recogidos en el mismo punto en el tiempo a traves de multiples sujetos, 
##        y esperamos que la correlacion entre cualquier par de sujetos sea la misma, esta estructura de simetria compuesta 
##        podria ser una buena eleccion.
##      Datos balanceados: La simetria compuesta suele ser un buen supuesto cuando se tienen datos balanceados, 
##        es decir, cuando cada sujeto tiene el mismo numero de observaciones. 
##  pdCorSymm: Esta funcion especifica una estructura de "simetria de correlacion". 
##    En esta estructura, todas las varianzas son iguales y todas las correlaciones entre diferentes efectos aleatorios son iguales. 
##    Difiere de pdCompSymm en que escala las covarianzas por las varianzas para obtener correlaciones.
##  pdBlocqued: Es util cuando se tienen diferentes "bloques" de variables que pueden tener diferentes estructuras de covarianza. 
##    Por ejemplo, si algunos de sus efectos aleatorios son medidas tomadas al mismo tiempo y otras se toman en momentos diferentes, 
##     es posible que deseemos permitir diferentes estructuras de covarianza dentro de cada bloque.

## random = list(id=~1, id=~x-1, id=~z-1)
##  equivalencia lmer: (1|id) + (0+x|id) + (0+z|id)    o     (1|id) + (x-1|id) + (z-1|id)
##  lmer(...) y lme(...) NO producen exactamente los mismos resultados
##  glmmTMBr(...) y lme(...) SI producen (casi) los mismos resultados 
##  Se asume que los interceptos y las pendientes aleatorios son independientes dentro de cada variable (x, z) para cada nivel de id.
##  Este modelo lme especifica un intercepto aleatorio y dos pendientes aleatorias no correlacionadas (para x y z) dentro del factor aleatorio id.
##  Permite que cada nivel de id tenga su propio intercepto, y su propia pendiente para x y z, 
##   pero estas pendientes no incluyen un termino de intercepto y se especifica que no estan correlacionadas con el intercepto aleatorio ni entre si.

## Diferencias entre estas dos estructuras de logCholesky:
##   random=list(id=~1+x+z)
##   random=list(id=~1, id=~x-1, id=~z-1)
##  list(id=~1+x+z) modeliza un intercepto aleatorio (~1) y pendientes aleatorias para x y z dentro de cada nivel del factor aleatorio id,
##     permitiendo que todos esos efectos aleatorios esten correlacionados. 
##     Se trata de un modelo multivariante de efectos aleatorios, con un intercepto aleatorio ~1 y las pendientes aleatorias para x y z dentro de cada nivel id. 
##  list(id=~1, id=~x-1, id=~z-1) modeliza un intercepto aleatorio (~1) y dos pendientes aleatorias para x y z como efectos aleatorios
##     separados dentro de cada nivel del factor de agrupacion id. 
##     El -1 en las formulas para x y z significa que los efectos aleatorios de x y z no estan correlacionados con el intercepto ni entre si. 
##     Implica modelizar tres efectos aleatorios univariantes separados para cada id.
##  Por tanto, la principal diferencia radica en los supuestos sobre las covarianzas entre los efectos aleatorios: 
##     list(id=~1+x+z) permite la covarianza entre efectos aleatorios, mientras que list(id=~1, id=~x-1, id=~z-1) no.

## Hay varias razones por las que los resultados de nlme::lme y lme4::lmer pueden no ser exactamente los mismos, 
##   incluso cuando se especifican modelos equivalentes:
##  Metodos de estimacion: Ambos paquetes utilizan tecnicas basadas en la maxima verosimilitud para la estimacion de parametros, 
##    pero implementan metodos especificos diferentes. 
##    La funcion lme utiliza un algoritmo de optimizacion de funciones denominado "algoritmo de Lindstrom-Bates" o "algoritmo de Pinheiro-Bates". 
##      Se trata de una combinacion del metodo Newton-Raphson y el metodo Fisher Scoring.
##      El metodo Newton-Raphson es un algoritmo de busqueda de raices que produce aproximaciones sucesivamente mejores 
##        a las raices de los efectos fijos y la matriz de varianza-covarianza de los efectos aleatorios.
##      El metodo de Fisher se usa para aproximar la matriz hessiana en las iteraciones. Es una tecnica iterativa para encontrar la estimacion 
##        de maxima verosimilitud de los parametros de un modelo estadistico.
##    La funcion lmer utiliza la "aproximacion de Laplace" o la "cuadratura adaptativa de Gauss-Hermite" para aproximar la funcion 
##        de log-verosimilitud en modelos mixtos. Tambien utiliza algoritmos de optimizacion como el metodo de Nelder-Mead, bobyqa, 
##        o el metodo de Gradiente Conjugado (dependiendo de los ajustes especificos) para encontrar los valores optimos de los 
##        parametros que maximizan la verosimilitud (ML) o la maxima verosimilitud restringida (REML).
##      La aproximacion de Laplace proporciona una forma adecuada de aproximar integrales complejas que aparecen en la funcion de verosimilitud de los modelos mixtos. 
##      La cuadratura adaptativa de Gauss-Hermite es un metodo mas preciso, pero computacionalmente mas intensivo, para aproximar estas integrales.
##    Las diferencias entre estos metodos pueden dar lugar a ligeras diferencias en los parametros estimados.
##    Si los algoritmos de de optimizacion convergen a soluciones ligeramente diferentes debido a cuestiones numericas, 
##      entonces dan lugar a estimaciones diferentes.
##  Los detalles reales de como se implementan y ajustan los efectos aleatorios pueden diferir entre lmer y lme.
##    https://www.flutterbys.com.au/stats/tut/tut9.1.html

## EXPLICACIONES, RESUMEN Y SUGERENCIAS
##  Los modelos mixtos lme(...) tienen algunas limitaciones a la hora de especificar terminos aleatorios complejos (efectos aleatorios cruzados).
##  Los modelos mixtos lmer(...) gestionan de manera sencilla esos efectos aleatorios cruzados (e.g.: 1|factor_aleatorio_1 + 1|factor_aleatorio_2)
##  Los modelos mixtos lme(...) utilizan una aproximacion clasica de grados de libertad de disennios balanceados ANOVA n-way
##  Los modelos mixtos lmer(...) utilizan dos aproximaciones diferentes: Satterthwaite y Kenward-Roger.
##   La de Kenward-Roger es la mas restrictiva respecto a los grados de libertad del denominador de los efectos fijos al usar "random slopes"
##  Los modelos mixtos lme(...) proporcionan resultados similares a los lmer con Satterthwaite (aunque los grados de libertad difieren).
##  Los modelos mixtos lmer(...) no gestionan la heterocedasticidad de los residuos. 
##  Los modelos mixtos lme(...) gestionan la heterocedasticidad de los residuos a traves del argumento weights=varPower(), weights=varExp(), varFixed(), varIdent(), ...
##  Los modelos mixtos lme(...) tambien gestionan situaciones de autocorrelacion (e.g., espacial y temporal) mediante el argumento cor=corAR1(), corARMA(), corEXp, corSpher, ...

## Debemos evitar situaciones en las cuales la varianza asociada a los terminos aleatorios sea muy pequennia, tendente a cero.
##   Este asunto conlleva estimas singulares, al ser la estructura de efectos aleatorios demasiado compleja para ser apoyada por los datos.
##   En modelos lmer(...) frecuentemente aparece esta alarma: boundary (singular) fit
##   en este caso se obtiene un ajuste "singular", indicativo de que el modelo esta sobreajustado 
##     o que los efectos aleatorios son muy pequennios (<1.0e-3)
##     (i.e., la estructura de efectos aleatorios es demasiado compleja para ser apoyada por los datos)
##   esto sugiere eliminar la parte mas compleja de la estructura de los efectos aleatorios (generalmente "random slopes"),
##     que conduce a un modelo mas parsimonioso, no sobreajustado.
##   Esta situacion poco deseable podria ocurrir en varios casos:
##     Sobreajuste: el modelo tiene demasiados terminos de efectos aleatorios dada la cantidad de datos y su estructura, 
##                  o se estima el mismo componente de varianza mas de una vez, puede ocurrir un ajuste singular.
##                  Como consecuencia, algunas de las pendientes aleatorias están perfectamente correlacionadas (+1 o -1) entre sí o con los interceptos aleatorios
##     "Sparse data": Cuando tenemos muy pocas observaciones por cada nivel de los factores aleatorios o poca variacion dentro de esos niveles.
##     Datos no balanceados: es decir, diferentes numeros de observaciones a traves de los niveles de los efectos aleatorios.
##   Como gestionar el problema de singularidad?
##     Simplificando la parte aleatoria: Si tenemos pendientes aleatorias en el modelo, podriamos querer eliminarlas o reducirlas.
##     Inspeccionando los datos: verificamos si los niveles en la estructura de efectos aleatorios estan balanceados y/o 
##        si hay suficiente variacion dentro de estos niveles. Si hay niveles con muy pocas observaciones, consideraremos eliminarlos o reunirlos con otros niveles.
##   Los ajustes singulares no son inherentemente problematicos y no invalidan nuestro modelo. 
##     De hecho, son menos frecuentes con nlme::lme que con lme4::lmer.
##     El problema simplemente sugiere que nuestro modelo podria ser demasiado complejo dado los datos, y deberiamos considerar un modelo mas sencillo. 

## El mensaje de advertencia "Model failed to converge with 1 negative eigenvalue" indica que el procedimiento de ajuste del modelo 
##   no logro encontrar exitosamente un conjunto de estimaciones de parametros que maximicen la verosimilitud de los datos 
##   dado el modelo. Esto se conoce como un problema de convergencia.
## Los problemas de convergencia en modelos mixtos generalmente ocurren cuando el modelo es complejo 
##   (es decir, tiene muchos parametros, especialmente parametros de pendientes aleatorias) en relacion con la 
##   cantidad de datos disponibles, los datos son escasos, o el modelo esta mal especificado.
## Los eigenvalues negativos son particularmente problematicos, ya que sugieren que la matriz de covarianza estimada 
##   para los efectos aleatorios no es definida positiva. Una matriz de covarianza debe ser semi-definida positiva 
##   (lo que significa que todos sus eigenvalues son mayores o iguales a cero) para ser una matriz de covarianza valida. 
##   Si una matriz de covarianza tiene un valor propio negativo, entonces algunos valores de la matriz
##   de varianzas-covarianzas es negativo, lo cual es imposible en la realidad.
## Esta alerta surge cuando:
##   Hay estructuras complejas de efectos aleatorios: Si el modelo incluye muchas pendientes aleatorias y/o correlaciones 
##     entre efectos aleatorios, y especialmente si estas no estan justificadas por los datos 
##     (i.e., no hay suficientes datos para estimarlos con precision).
##   Hay escasez de datos: si tenemos muchos niveles de un efecto aleatorio con pocas observaciones por nivel.
##   Se han estimado correlaciones entre efectos aleatorios cercanas a -1 o 1.
##     Una correlacion perfecta o casi perfecta entre efectos aleatorios (e.g., > 0.95 o < -0.95) puede llevar a 
##     la inestabilidad numerica, resultando en autovalores negativos.
## Como afrontar este problema?
##   Simplificando el modelo, reduciendo la complejidad de la estructura de efectos aleatorios. 
##     Se puede lograr eliminando algunos (o todos) efectos de pendientes aleatorias, 
##     o asumir su independencia no estimando las correlaciones entre esos efectos aleatorios (sean entre pendientes y/o interceptos).
##   Revisando los datos para valorar si los niveles de un efecto aleatorio tienen muy pocas observaciones, 
##     considerando si esos niveles pueden ser eliminados (o combinados con otros).
##   Aumentando el numero de las iteraciones: A veces, aumentar el numero maximo de iteraciones permitidas para el algoritmo de optimizacion 
##     (control = lmerControl(optCtrl = list(maxfun = 1e5)) en lme4::lmer) puede resolver el problema.
##   Usando diferentes optimizadores, que pueden especificarse en la funcion lmerControl
##     usando "Nelder-Mead", "BFGS", "nloptwrap", "optimx". 

## Escalamiento/Centrado de las predictoras, para lograr estabilidad numerica en las estimas de los predictores.
## Con las funciones `lme4::lmer` y `nlme::lme` estandarizar las variables a veces puede hacer que la estimacion del modelo sea numericamente mas estable, 
##   particularmente en modelos mixtos con muchos terminos aleatorios y de orden superior (como interacciones entre factores aleatorios),
##   o cuando las escalas de los predictores son muy diferentes (e.g., pH y altitud en m).
## Se suele necesitar la estandarizacion cuando aparecen algunos problemas potenciales como:
##   1. Advertencias de convergencia: Tanto `lme4::lmer` como `nlme::lme` utilizan metodos iterativos para estimar los parametros del modelo. 
##      Si estos metodos no convergen, la estimacion de los parametros no se ha estabilizado, lo cual puede ser una señal de que
##      los predictores necesitan ser reescalados.
##   2. Coeficientes muy grandes o muy pequeños: Si las escalas de los predictores son muy diferentes (en decimas, en unidades o en miles), 
##      algunos coeficientes pueden terminar siendo muy grandes o muy pequeños. Esto puede ser una señal de que la estandarizacion podria ser de utilidad.
##   3. Advertencias de alta correlacion entre terminos aleatorios: Si estamos utilizando interacciones, los predictores no estandarizados a veces pueden conducir 
##      a una alta multicolinealidad en modelos de pendientes aleatorias. Estandarizar los predictores a veces puede aliviar este problema.
## En resumen, puede ser una buena practica estandarizar las variables predictoras en modelos mixtos, ya que puede facilitar la interpretacion de los coeficientes, 
##   reducir la inestabilidad numerica en la busqueda de convergencia y evitar problemas con la multicolinealidad. 
## Ejemplos de covariantes a incluir en la formula: X --> I(scale(X));    X+X^2+X^3 --> I(poly(X, 3))
plot(datos$t_sonda, I(scale(datos$t_sonda)))
pairs(I(poly(datos$t_sonda, 3)), main="TODOS ESCALADOS    1: lineal, 2: cuadratico, 3: cubico")
plot(I(poly(datos$t_sonda, 3))[,1], I(poly(datos$t_sonda, 3))[,2], xlab="variable X escalada a media=0", ylab="variable X^2 escalada a media=0")
plot(I(poly(datos$t_sonda, 3))[,1], I(poly(datos$t_sonda, 3))[,3], xlab="variable X escalada a media=0", ylab="variable X^3 escalada a media=0")
## En el caso de poly(X, 3), R genera tres terminos polinomiales: x, x^2, y x^3, cada uno de los cuales esta ortogonalizado con respecto a los otros. 
##   La media de esos terminos es cero. En cuanto a la desviacion tipica (sd), no se lleva a un valor especifico como en el caso de la estandarización (sd=1). 
##   La ortogonalización garantiza que la varianza de cada termino polinomial es constante.
##   Por tanto, aunque los terminos polinomiales generados por poly(X, 3) estan ortogonalizados (y por lo tanto, no estan correlacionados entre si), 
##   esto no significa que esten estandarizados (con media cero y desviacion estandar, sd, igual a uno).

##   Para un modelo mixto con interceptos y pendientes aleatorias, el tamannio de la muestra deberia ser idealmente lo suficientemente grande
##     como para que cada grupo (definido por los niveles de los factores aleatorios) tenga un numero razonable de observaciones. 
##     Esto se debe a que los interceptos y pendientes aleatorias se estiman en funcion de la variacion dentro de estos grupos.
##     Una regla general podria ser tener al menos 5-10 observaciones por grupo, y al menos 10-20 grupos, siendo esos grupos los
##     niveles de los factores aleatorios. No obstante esos números dependen de las especificidades del estudio, incluyendo 
##     la magnitud de los efectos y la variabilidad en los datos.

## Los modelos mixtos glmmTMB y lme producen resultados muy similares.

## Si no hemos tenido un buen ajuste a los supuestos canonicos de los modelos mixtos (normalidad, homocedasticidad de los residuos, puntos outliers)
##   * construiremos hipotesis nulas para examinar la magnitud de los coeficientes de los efectos fijos (cuyos valores deben estar fuera de los intervalos de confianza nulos)
##   * obtendremos los parametros de los modelos mixtos mediante bootstraps (cuyos intervalos de confianza no deben incluir el valor cero)

## Anova{car}: (Type III Wald F tests with Kenward-Roger df)
## http://ir.library.oregonstate.edu/xmlui/bitstream/handle/1957/5262/mydissertation.pdf
## Instead of approximating the Wald-type statistic itself by an F distribution as in the Satterthwaite approximation, 
## Kenward and Roger considered a scaled form of the statistic that makes the approximation more flexible since
## they approximate the denominator degrees of freedom and the scale as well. Second, Kenward and Roger
## modified the procedure in such a way that they get the right values for both the denominator degrees of freedom
## and the scale for the two special cases considered in chapter 3. Third, the adjusted estimator of the
## variance-covariance matrix of the fixed effects estimator, which is less biased than the conventional estimator, 
## as we saw in chapter 2, is used in constructing the Wald-type statistic.
## The K-R, the Satterthwaite, and the proposed methods produce the same estimate of the denominator degrees of 
## freedom, and the scale estimate is one. Even though the estimates are the same, the Satterthaite and the K-R tests
## are not necessarily identical, and this is true because in the K-R statistic, we use the adjusted estimator of the
## variance-covariance matrix of the fixed effects estimator, whereas the conventional estimator is used in the Satterthwaite test.
## We conducted a simulation study for three kinds of block designs; partially balanced block designs, balanced incomplete block
## designs, and complete block designs with missing data. The K-R, the proposed, the Satterthwaite, and the Containment methods
## were compared in the study. Three factors were considered in the study: the variance components ratio, the efficiency factor,
## and the sample size. In most cases, we found the K-R and the proposed methods performed as well as or better than the other
## methods. For designs with small values for all the three factors, the Satterthwaite method performed poorly.

## Kenward-Roger modification of the F-statistic for some linear mixed models fitted with lmer
## http://web.warwick.ac.uk/statsdept/useR-2011/abstracts/290311-halekohulrich.pdf





##############################
## BOOTSTRAPPING DEL MODELO ##
##############################
library(bayestestR)
## Model-based (Semi-)Parametric Bootstrap for Mixed Models; para modelos lmer
mcmc.fixed <- bootMer(modelo, FUN=fixef, nsim=500, type="parametric", use.u=FALSE, seed=1, verbose=TRUE)
coef.mcmc.fixed <- as.data.frame(mcmc.fixed$t)
## resumen con coeficientes del modelo de interes (original), la mediana de los coeficientes del bootstrap (bootMed),
##   su error estandard (bootSE) y el bias (bootBias)
summary(mcmc.fixed)
## tabla de resultados con cuantiles para alfa 0.05 y 0.01
tabla.mixedmodel.simulada <- describe(coef.mcmc.fixed, quant=c(0.025, 0.975, 0.005, 0.995))[,c(1:4, 14:17)]
colnames(tabla.mixedmodel.simulada)[4] <- "std_error"
tabla.mixedmodel.simulada$coeficiente.modelo <- round(fixef(modelo), 5)
tabla.mixedmodel.simulada$std_error.modelo <- round(summary(modelo, ddf="Kenward-Roger")$coefficients[,2], 5)
tabla.mixedmodel.simulada$P.modelo <- round(summary(modelo, ddf="Kenward-Roger")$coefficients[,5], 5)
## comparar los coeficientes obtenidos por bootstrapping (en mean) con los originales del modelo (coeficiente.modelo)
##   si los valores de kurtosis y sesgo no son cercanos a cero eso es indicativo de la existencia de valores extremos (influyentes, perdidos, ...)
## utilizaremos los intervalos de confianza derivados de cuantiles (95%: Q0.025 y Q0.975; 99%: Q0.005 y Q0.995)
##   que no incluyan el valor "cero"
print(tabla.mixedmodel.simulada, digits=5)
## intervalos de confianza ETI y HDI
##   que los intervalos no incluyan el valor CERO
eti(coef.mcmc.fixed, ci=0.95)
plot(eti(coef.mcmc.fixed, ci=0.95))
hdi(coef.mcmc.fixed, ci=0.95)
plot(hdi(coef.mcmc.fixed, ci=0.95))
## Valoracion del grado de correlacion entre los coeficientes del modelo
##  los valores de r obtenidos se deberian aproximar a "cero" si los efectos fijos fuesen independientes entre si
correlaciones <- cor(coef.mcmc.fixed)
corrplot(correlaciones, type="upper", method="circle")   ## para valores numericos de las correlaciones incluimos: method="number"
