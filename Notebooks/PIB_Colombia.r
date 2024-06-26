# Ejemplo de identificaci�n, estimaci�n,  diagn�sticos y pron�sticos en un modelo ARIMA
#  
# ----------------------------------------------------------------------------
rm(list=ls(all=TRUE))  # elimina todos los objetos creados que todav�a est�n en R

# Cargar librer�as
library(fpp2)
library(TSA)
library(FitAR)
library(urca)            
library(fBasics)

# lectura de datos desde un archivo de texto. Los datos corresponden
# al Producto Interno Bruto de Colombia observado trimestralmente y desestacionalizado
#  desde el primer trimestre del a�o 2000 al cuarto trimestre del a�o 2015.

#(z=ts(scan("D:/Curso_Series_Financieras/datos/PIB_real_2000_1_2015_4.txt")))
(z=ts(scan(file.choose())))
length(z)   # la serie es corta. Esto puede dificultar la identificaci�n del modelo.

##### ETAPA DE IDENTIFICACI�N #####

# gr�fica de la serie
plot.ts(z, type="o", cex=0.6)
#
### An�lisis de estabilidad de la varianza: Es necesario usar la ###
### transformaci�n Box-Cox? ###

# An�lisis gr�fico
par(mfrow=c(2,1))
plot.ts(z, type="o", cex=0.6)
plot.ts(diff(z), type="o", cex=0.6)
dev.off()

# usando la librer�a forecast: supone que los datos son dependientes
BoxCox.lambda(z, method = c("loglik"), lower = -2, upper = 2)  # maximiza log likelihood
# sugiere una transformaci�n logaritmica

# comparaci�n de la serie original y la transformada
par(mfrow=c(2,1))
plot.ts(z, type="o", cex=0.6)
plot.ts(log(z), type="o", cex=0.6)
dev.off()

# comparaci�n de las primeras diferencia de serie original y la transformada
par(mfrow=c(2,1))
plot.ts(diff(z), type="o", cex=0.6)
plot.ts(diff(log(z)), type="o", cex=0.6)
dev.off()

# se observa que la transformaci�n logaritmica no parece tener efecto.

# otra alernativa usando la librer�a forecast
BoxCox.lambda(z, method = c("guerrero"), lower = -2, upper = 2)  #minimiza coef. de variaci�n

# usando la librer�a FitAR: Enc encuentra lambda en un modelo AR(p) aproximado y diferenciado si se requiere.
(p=round(length(z)^(1/3)))
mod1=arima(z, c(p, 1, 0))
BoxCox(mod1)

# Conclusi�n: se propone usar la transformaci�n lambda=1, es decir, no se transforma la serie.  
# y se sigue trabajando con la serie original z

### determinaci�n del valor de d ###
# correlogramas muestrales de la serie original
(max_rezag=round(length(z)/4))
par(mfrow=c(2,1))
Acf(z, lag.max=max_rezag, ci=0,ylim=c(-1,1))    # ACF decrece lentamente
Pacf(z, lag.max=max_rezag, ci=0, ylim=c(-1,1))  # PACF concentrada en el primer rezago
dev.off()

# del gr�fico de z y sus correlogramas se concluye que la serie es no estacionaria:
# Qu� tipo de serie es: TS o DS?

### prueba de ra�ces unitarias ###

# usando la prueba aumentada de Dickey-Fuller de la librer�a urca
(maxlag=floor(12*(length(z)/100)^(1/4)))

## Bajo el modelo de selecci�n AIC ##
ru_z=ur.df(z, type = c("trend"), lags=maxlag, selectlags = c("AIC")) # Usando AIC
summary(ru_z)
# valide el modelo elegido usando AIC

## Bajo el criterio de selecci�n BIC ##
ru_z=ur.df(z, type = c("trend"), lags=maxlag, selectlags = c("BIC")) # Usando BIC
summary(ru_z)    # parece que el coeficiente estimado de la primera diferencia rezagada no
                 # es estad�sticamente significativo.

# redefinici�n del modelo basado en el criterio BIC
ru_z=ur.df(z, type = c("trend"), lags=0)   # redefinici�n del modelo
summary(ru_z)
# valide el modelo elegido usando BIC

## Validaci�n de los modelos ##
# se chequea cada modelo y se usa la prueba usando el modelo sugerido por el criterio AIC
  # validaci�n de la ecuaci�n ADF
    resid=ru_z@testreg$residuals             # residuales del modelo ADF
    plot(ru_z)                               # diagn�sticos 
    auto.arima(resid, max.p=5, max.q=5)      # busqueda "autom�tica" 
    cheq=Arima(resid, c(0,0,0), include.constant=TRUE) 
    summary(cheq)    
    tsdiag(cheq, gof=12) # Prueba conjunta para Rho_k
    
  # Verificacion de normalidad en los residuales
    resid_sd=resid/ru_z@testreg$sigma        # estandarizaci�n de los residuales
    qqnorm(resid_sd,  xlab = "Cuantiles Te�ricos", ylab = "Cuantiles Muestrales")
    abline(a=0, b=1)

    shapiro.test(resid_sd)             #prueba de Shapiro-Wilks
    jarqueberaTest(resid_sd )          # prueba Jarque-Bera 

# el modelo elegido, basado en el criterio AIC, para la prueba ADF satisface los supuestos b�sicos.
# Bajo el criterio BIC la prueba de Ljung-Box dice que hay al menos un coeficiente de correlaci�n
# que no es cero

# Conclusi�n: Parece que la serie z contiene al menos una ra�z unitaria.

### Examen de m�s ra�ces unitarias ###
# gr�fica de la serie diferenciada una vez
plot.ts(diff(z), type="o")

# prueba de si hay ra�z unitarias en z diferenciada una vez.
ru_dif_z=ur.df(diff(z), type = c("drift"), lags=maxlag, selectlags = c("AIC"))
summary(ru_dif_z)  

# redefinici�n del modelo para la prueba
ru_dif_z=ur.df(diff(z), type = c("drift"), lags=0)
summary(ru_dif_z)  

# se validan y se elige el primer modelo
  # validaci�n de la ecuaci�n ADF
    resid1=ru_dif_z@testreg$residuals             # residuales del modelo ADF
    plot(ru_dif_z)                                # diagn�sticos 
    auto.arima(resid1, max.p=5, max.q=5)          # busqueda "autom�tica" 
    cheq1=Arima(resid1, c(0,0,0), include.constant=TRUE) 
    summary(cheq1)    
    tsdiag(cheq1, gof=12)
  # Verificacion de normalidad en los residuales
    resid1_sd=resid/ru_dif_z@testreg$sigma        # estandarizaci�n de los residuales
    qqnorm(resid1,  xlab = "Cuantiles Te�ricos", ylab = "Cuantiles Muestrales")
    qqline(resid1)

    shapiro.test(resid1)                #prueba de Shapiro-Wilks
    jarqueberaTest(resid1)             # prueba Jarque-Bera 

# Conclusi�n: se cumplen los supuestos de la prueba. 
# La serie z diferenciada no contiene r�iz unitaria.
# Por tanto, la serie z pertenece a la clase de modelos ARIMA con deriva y d=1:
#     phi(B)Wt=theta0+theta(B)at  con thetao>0  

### determinaci�n de los valores de (p, q) del modelo ARMA para Wt=(1-B)zt ###
### correlogramas muestrales para z diferenciada una vez ###
par(mfrow=c(2,1))
Acf(diff(z), lag.max=12, ci=0, ylim=c(-1,1))
Pacf(diff(z), lag.max=12, ci=0, ylim=c(-1,1))

lag.plot(diff(z), lags=4, do.lines=F, diag=F)

# correlogramas muestrales  para (1-B)z con bandas
par(mfrow=c(2,1))
Acf(diff(z), lag.max=16)
Pacf(diff(z), lag.max=16)

# Parece que el proceso que gener� la serie z es un ARIMA(0,1,0).
# selecci�n usando criterios de informaci�n: use script criter_info_PIB_Col.r
# usando el criterio de informaci�n AIC indica el modelo ARIMA(0,1,3). 
# usando el criterio de informaci�n BIC indica el modelo ARIMA(0,1,0). 

# selecci�n "autom�tica" del modelo usando librer�a foreast
auto.arima(z, d=1, max.p=5, max.q=5, ic=c("aic"))
auto.arima(z, d=1, max.p=5, max.q=5, ic=c("bic"))
# Usando el criterio de informaci�n AIC parece que el proceso que gener� 
# la serie es un ARIMA(0,1,3) con deriva.
# Usando el criterio de informaci�n BIC parece que el proceso que gener� 
# la serie es un ARIMA(0,1,0) con deriva.

# selecci�n del modelo usando la eacf: se elige el modelo que se�ala el v�rtice 
# de un tri�ngulo de ceros en la tabla                      
eacf(diff(z))    # se encuentra en el paquete TSA
# Hay evidencia de que el proceso que gener� la serie es un ARIMA(0,1,0) o un
# ARIMA(1,1,0).

##### ETAPA DE ESTIMACI�N #####
# ARIMA(0,1,0): estimaci�n ML exacta del con valores iniciales dados por la estimaci�n condicional
mod1_CSS_ML=Arima(z, c(0, 1, 0), include.drift=TRUE, lambda=1, method = c("CSS-ML"))
summary(mod1_CSS_ML) 
(res1_CSS_ML=residuals(mod1_CSS_ML))
res1_est=res1_CSS_ML/(mod1_CSS_ML$sigma2^.5)  # estandarizaci�n de los residuales

# ARIMA(1,1,0): estimaci�n ML exacta del con valores iniciales dados por la estimaci�n condicional
mod3_CSS_ML=Arima(z, c(1, 1, 0), include.drift=TRUE, lambda=1, method = c("CSS-ML"))
summary(mod3_CSS_ML) 
(res3_CSS_ML=residuals(mod3_CSS_ML))
res3_est=res3_CSS_ML/(mod3_CSS_ML$sigma2^.5)  # estandarizaci�n de los residuales

# ARIMA(0,1,3): estimaci�n ML exacta del con valores iniciales dados por la estimaci�n condicional
mod4_CSS_ML=Arima(z, c(0, 1, 3), include.drift=TRUE, lambda=1, method = c("CSS-ML"))
summary(mod4_CSS_ML) 
(res4_CSS_ML=residuals(mod4_CSS_ML))
res4_est=res4_CSS_ML/(mod4_CSS_ML$sigma2^.5)  # estandarizaci�n de los residuales

###### ETAPA DE DIAGN�STICOS #####

## ra�ces de los polinomios: ##
# ARIMA(0,1,0) NO HAY POLINOMIOS AR ni MA.

# ARIMA(1,1,0)HAY POLINOMIO AR.
#ra�ces del polinomio AR
autoplot(mod3_CSS_ML)  

# ARIMA(0,1,3)HAY POLINOMIO MA.
#ra�ces del polinomio MA
autoplot(mod4_CSS_ML)  

## An�lisis de los residuales ##
# ARIMA(0,1,0)
checkresiduals(mod1_CSS_ML, lag=15)   # de la librer�a forecast
tsdiag(mod1_CSS_ML, gof=15)           # de la librer�a stats
# parece que el modelo deja dudas en explicar bien toda la estructura de dependencia
# aqu� se detendr�an los diagn�sticos para ese modelo y se tratar�a de identificar un nuevo modelo.

# ARIMA(1,1,0)
tsdiag(mod3_CSS_ML, gof=15)
# parece que el modelo deja dudas en explicar bien toda la estructura de dependencia
# aqu� se detendr�an los diagn�sticos para ese modelo y se tratar�a de identificar un nuevo modelo.

# ARIMA(0,1,3)
tsdiag(mod4_CSS_ML, gof=15)
# No se rechaza que los residuales del modelo son ruido blanco 
# En adelante se continua con los otros diagn�sticos para este mdelo

# chequeo de normalidad
# gr�fico cuantil-cuantil
qqnorm(res4_est,  xlab = "Cuantiles Te�ricos", ylab = "Cuantiles Muestrales")
abline(a=0, b=1)
 
# histograma, densidad kernel y gr�fico normal
hist(res4_CSS_ML, prob=T)
lines(density(res4_CSS_ML), lwd=3, col="orange")
curve(dnorm(x, mean=mean(res4_CSS_ML), sd=sd(res4_CSS_ML)), col="darkblue", lwd=2, add=TRUE)

# pruebas de normalidad
shapiro.test(res4_est)                 # prueba de Shapiro-Wilks

# se usa la librer�a fBasics para la descripci�n de la distribuci�n de los residuales, incluyendo
# los coeficientes de asimetr�a y curtosis muestrales y para la prueba de ljung_Box

basicStats(res4_est)  # kurtosis corresponde al "exceso de curtosis"=curtosis-3 
normalTest(res4_est, method=("jb"))  # librer�a fBasics: puede realizar otras pruebas
                                     # "ks" for the Kolmogorov-Smirnov one�sample test, 
                                     # "sw" for the Shapiro-Wilk test, 
                                     # "jb" for the Jarque-Bera Test, 
                                     # "da" for the D'Agostino Test. The default value is "ks"

# conclusi�n: No se rechaza la normalidad, a un nivel de significancia de alpha=0.05.
# Sin embargo, si se usa alpha>=0.09, la prueba de Shapiro-Wilk rechazar�a normalidad. 

# detecci�n de observaciones at�picas distantes
# chequeo de observaciones at�picas extremas (no es un an�lisis completo de outliers)
plot.ts(res4_est, type="o", ylim=c(-4,4))
abline(a=-3, b=0, col="red", lty=2)
abline(a=3, b=0, col="red", lty=2)

ind=(abs(res4_est)>3)
sum(ind)					# n�mero de observaciones at�picas
(grupo=cbind(seq(z),res4_est, ind))

# valores ajustados del modelo para z
(ajust=mod4_CSS_ML$fitted)     
# gr�fico para los valores ajustados y los valores observados
ts.plot(z,ajust)   # gr�fico de las series contra el tiempo
lines(z, col="black", type="o", cex=.5)
lines(ajust, col="red", type="o", cex=.5)

# gr�fico de dispersi�n de valores observados vs valores ajustados
plot(as.vector(ajust), as.vector(z), type="p")   # gr�fico de dispersi�n de la serie observada 
abline(0,1, col="red")                          # contra la serie ajustada

# Evaluaci�n de la significancia estad�stica de los coeficientes estimados
# Construcci�n de los estad�sticos de prueba, t
# estad�sticos t
summary(mod4_CSS_ML) 
(t=coef(mod4_CSS_ML)/(diag(vcov(mod4_CSS_ML)))^.5)
# valores P
(val_p=2*pnorm(t, lower.tail=FALSE))
round(val_p, 4)
# a un nivel de significancia(aproximado) de alpha=0.10, el coficiente ma1 no es sigificativo;
# pero los coeficientes del ma2  y ma3 s� son significativos.

# estimaci�n del ARIMA(0,1,3) con el coeficiente ma1 restringido a 0.
#-------------------------------------------------------------------
# ARIMA(0,1,3): estimaci�n ML exacta del con valores iniciales dados por la estimaci�n condicional
mod4_CSS_ML=Arima(z, c(0, 1, 3), include.drift=TRUE, lambda=1, method = c("CSS-ML"),
fixed=c(0,NA,NA,NA))
summary(mod4_CSS_ML) 

autoplot(mod4_CSS_ML)
tsdiag(mod4_CSS_ML)
(res4_CSS_ML=residuals(mod4_CSS_ML))
res4_est=res4_CSS_ML/(mod4_CSS_ML$sigma2^.5)  # estandarizaci�n de los residuales

# chequeo de normalidad
# gr�fico cuantil-cuantil
qqnorm(res4_est,  xlab = "Cuantiles Te�ricos", ylab = "Cuantiles Muestrales")
abline(a=0, b=1)
 
# pruebas de normalidad
shapiro.test(res4_est)                 # prueba de Shapiro-Wilks
normalTest(res4_est, method=("jb"))    # prueba de Jarque-Bera

# estad�sticos t
(t=coef(mod4_CSS_ML)[2:4]/(diag(vcov(mod4_CSS_ML)))^.5)
# valores P
(val_p=2*pnorm(t, lower.tail=FALSE))
round(val_p, 4)
# todos los coeficientes son significativos a un nivel de significancia=0.09.
qnorm(.045, lower.tail=F)

# valores ajustados del modelo para z transformada
(ajust=mod4_CSS_ML$fitted)     
# gr�fico para los valores ajustados y los valores observados
ts.plot(z,ajust)   # gr�fico de las series contra el tiempo
lines(ajust, col="red", type="o", cex=.5)

# gr�fico de dispersi�n de valores observados vs valores ajustados
plot(as.vector(z),as.vector(ajust), type="p")   # gr�fico de dispersi�n de la serie observada 
abline(0,1, col="red")                          # contra la serie ajustada

# Conclusi�n: el modelo pasa los diagn�sticos.
# modelo elegido: un IMA(1,3) dado por  (1-B)Zt=theta0+(1-theta2*B^2-theta3*B^3)at
# El PIB evoluciona con una tendencia aleatoria mezclada con una tendencia determin�stica lineal.
# mas una componente estacionaria ma(3) restringida con theta1=0.

# Es el modelo �til para pronosticar?
#
# ======================= EVALUACI�N DE LOS PRON�STICOS =======================
#
# se contruye el modelo usando 56 datos (conjunto de entrenamiento) y se dejan
# los �ltimos 8 (conjunto de prueba) para evaluar la capacidad de pron�stico del modelo 
# estimaci�n ML exacta con valores iniciales dados por la estimaci�n condicional
# modelo sin restricci�n
(mod_Evalpron=Arima(z[1:56], c(0, 1, 3), include.drift=T, method = c("CSS-ML")))
summary(mod_Evalpron) # compare con el modelo con todos los datos

# modelo restringido
mod_Evalpron=Arima(z[1:56], c(0, 1, 3), include.drift=TRUE, lambda=1, method = c("CSS-ML"),
fixed=c(0,NA,NA,NA))
summary(mod_Evalpron)

# Diagn�stico sobre el modelo estimado #
autoplot(mod_Evalpron)
res_eval=residuals(mod_Evalpron)
tsdiag(mod_Evalpron, gof=14)
checkresiduals(mod_Evalpron)
qqnorm(residuals(mod_Evalpron))
qqline(residuals(mod_Evalpron))
shapiro.test(residuals(mod_Evalpron))

# Pron�sticos #
(z_pred=forecast(mod_Evalpron, h=8, level=c(80, 95), fan=FALSE))

# gr�fico de los pron�sticos
plot(z_pred)

# o, alternativamente,
autoplot(z_pred)

# Evaluaci�n de los pr�nosticos
# lista de los valores reales y los pron�sticos
cbind(z[57:64], z_pred$mean)

ts.plot(z[57:64], z_pred$mean, z_pred$lower[,2], z_pred$upper[,2], type="o")
lines(z_pred$mean, col="red")
lines(z_pred$lower[,2], col="blue")
lines(z_pred$upper[,2], col="blue")

# precisi�n de los pron�sticos
# --------------------------------------------------
(recm=(mean((z[57:64]-ts(z_pred$mean))^2))^.5)            # ra�z cuadrada error cuadr�tico medio 
(recmp=100*(mean(((z[57:64]-ts(z_pred$mean))/z[57:64])^2))^.5) # ra�z cuadrada error cuadr�tico medio en porcentaje

(eam=mean(abs(z[57:64]-ts(z_pred$mean))))                 # error absoluto medio
(eamp=100*mean(abs((z[57:64]-ts(z_pred$mean))/z[57:64]))) # error absoluto medio en porcentaje

# descomposici�n del error cuadr�tico medio
# --------------------------------------------------
# error cuadr�tico medio
(ecm=mean((z[57:64]-ts(z_pred$mean))^2)) 

# proporci�n de sesgo en media
(prop_sesgom=(mean(z_pred$mean)-mean(z[57:64]))^2/ecm) 

# c�lculo de varianza sesgadas
(num_pron=length(z_pred$mean))
(sigmap=(((num_pron-1)/num_pron)*var(z_pred$mean))^.5) # desv. est�nd de los pron�sticosdatos 
(sigmar=(((num_pron-1)/num_pron)*var(z[57:64]))^.5)  # desv. est�nd de los datos

# proporci�n de sesgo en varianza
(prop_var=(sigmap-sigmar)^2/ecm)         

# proporci�n de covarianza
(prop_covar=2*(1-cor(z_pred$mean, z[57:64]))*sigmap*sigmar/ecm)

(prop_sesgom+prop_var+prop_covar)  # Chequeo

# EJERCICIO: Construya un modelo para el PIB usando la transformaci�n logar�tmica y compare
# sus resultados con el modelo anterior.

# Verdaderos Pron�sticos
# --------------------------------------------------
par(mfrow=c(1,2))
(z_pron<-forecast(mod4_CSS_ML, h=8, level=c(80,95),  fan=F))
plot(z_pron)

(z_pron1<-forecast(mod4_CSS_ML, h=8, level=c(80,95),  fan=T))
plot(forecast(z_pron))

# o, alternativamente
autoplot(forecast(z_pron))
autoplot(forecast(z_pron1))

