library(forecast)

# Entrada de datos desde un archivo de texto
(z=ts(scan("D:/Curso_exten_series_UdeA/datos/PIB_real_2000_1_2015_4.txt")))

# selección del modelo usando criterios de información: AIC o BIC
# defina d, p y q
d=1
p=4
q=4

(num_modelos=(p+1)*(q+1))
aic=matrix(rep(-99, times=num_modelos*3), nrow=num_modelos, ncol=3)
bic=matrix(rep(-99, times=num_modelos*3), nrow=num_modelos, ncol=3)

k=1
for(i in 0:p) {
for(j in 0:q) {
mod=Arima(z, c(i,d,j), include.drift=T, method=c("CSS-ML"))
aic[k, 1]=i
aic[k, 2]=j
aic[k, 3]=mod$aic

bic[k, 1]=i
bic[k, 2]=j
bic[k, 3]=mod$bic
k=k+1
}
}

aic=data.frame(aic)
attach(aic)
p=X1
q=X2
Aic=X3
(AIC=data.frame(cbind(p, d, q, Aic)))

bic=data.frame(bic)
attach(bic)
p=X1
q=X2
Bic=X3
(BIC=data.frame(cbind(p, d, q, Bic)))

cbind(AIC[order(Aic),], "     ", BIC[order(Bic),])

which.min(AIC$Aic)#Posición del Valor mínimo segun AIC
which.min(BIC$Bic)#Posición del Valor mínimo segun BIC
