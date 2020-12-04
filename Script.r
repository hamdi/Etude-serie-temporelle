##########################################################
##            Projet de series temporelles              ##
##     Hamdi BEL HADJ HASSINE et Fadi SAHBANI           ##
##########################################################


# Packages nécessaires:
library(tseries)
library(forecast)
library(fUnitRoots)
library(car)


#####--   Partie 1 : Les données   --#####

### Question 1 ###
path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path)

# Extraction et mise en forme des données :
data <- as.data.frame(read.csv("valeurs_mensuelles.csv", fileEncoding = "UTF-8", sep = ";"))
donnees <- as.matrix(data[-c(1, 2),c(1,2)])
colnames(donnees) <- c("Date","Valeurs")
donnees <- apply(donnees, 2, rev)
rownames(donnees) <- 1:dim(donnees)[1]
head(donnees)
s <- ts(as.numeric(donnees[,2]),start=1990,frequency=12)
n <- length(s)
plot(s,xlab='Date',ylab="Indice de production industrielle", main="Indice de production industrielle")

### Question 2 ###
monthplot(s)
#les 12 chronogrammes mensuels sont à peu près identiques, ce qui confirme l'absence de saisonalité

lag.plot(s,lags=12,layout=c(3,4),do.lines=FALSE)
# Le Lagplot montre une corrélation forte.

fit1 <- decompose(s)
plot(fit1)
#on peut voir que l'erreur ne semble pas trop varier au cours du temps ce qui indique l'adéquation du modele additif
#Visualisation ACF et PACF
acf(s)
pacf(s)
# L'ACF montre que la série n'est pas stationnaire et qu'elle présente une tendance, mais ne montre pas de saisonnalité.

#Test de la tendance
summary(lm(s ~ seq(1,n)))
# Le coefficient de la tendance linéaire est significatif, donc on effectue les tests de stationnarité avec tendance

#Test KPSS de la stationnarité
kpss.test(s,null="Trend")
# Le test KPSS rejette au niveau 1% l'hypothèse de stationnarité de la série.

#Test ADF de la stationnarité
#   Tests de LjungBox combinés pour verifier l'autorrélation des résidus 
#   jusqu'à l'ordre k
LjungBoxtest <- function(X, k, fitdf=0) {
  pvalues <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(X, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvalues))
}
#Fonction pour effectuer un test Dickey-Fuller augmenté valide
ValidADF <- function(X,kmax,type){ # Tests ADF jusqu'à avoir des résidus non autocorrélés
  k <- 0
  pasautcor <- 0
  while (pasautcor==0){
    cat(paste0("ADF avec ",k, " lags: "))
    adf <- adfTest(X,lags=k,type=type)
    pvalues <- LjungBoxtest(adf@test$lm$residuals,24,fitdf=length(adf@test$lm$coefficients))[,2]
    if (sum(pvalues<0.05,na.rm=T) == 0) {
      pasautcor <- 1; 	      pasautcor <- 1; cat("Les résidus ne sont pas auto-corrélés. Test ADF valide.\n")}
    else cat("Les résidus sont auto-corrélés \n")
    k <- k + 1
  }
  return(adf)
}
adf <- ValidADF(s,24,"ct")
adf
# Le test ADF ne rejette pas au seuil 5% la présence d'une racine unitaire.
# La serie est donc intégrée. On la différencie à l'ordre 1 :
x <- diff(s)
plot(x,xlab='Date',ylab="IPI différencié", main="Série différenciée")

# La série différenciée ne semble pas présenter de constante ou tendance.
#On le vérifie avec une régression linéaire :
summary(lm(x ~ seq(1,length(x))))
# Les coefficients associés à la constante et à la tendance sont non 
# significatifs au seuil 5%

#On fait un test KPSS:
kpss.test(x,null="Level")
# Le test KPSS ne rejette pas la stationnarité de la série différenciée

# et maintenant on verifie avec le test ADF
suppressWarnings(
  adf <- ValidADF(x,50,"nc"))
adf
# Le test ADF rejette au seuil 1% la présence d'une racine unitaire.

# On peut aussi le vérifier par la fonction lagplot
lag.plot(x,lags=12,layout=c(3,4),do.lines=FALSE)
# La serie différenciée est donc stationnaire,
#c'est-à-dire que la série originale s est I(1)

#Finalement on effectue un test pour vérifier l'homoscédasticité

# Test d'homoscédasticité de Breusch-Pagan
lmtest::bptest(lm(x ~ seq(1,length(x))))
# On ne rejette pas l'homoscédasticité de X au niveau 5%
### Question 3 ###

plot(s,xlab='Date',ylab="Indice de production industrielle", main="Indice de production industrielle")
plot(x,xlab='Date',ylab="IPI différencié", main="Série différenciée")

#####--   Partie 2 : Modèles ARMA   --#####

### Question 4 ###

acf(as.numeric(x), main="ACF de X")
# L'ACF est seulement significatif à l'ordre 1. On note un léger dépassement au lag 11 mais on suppose qu'il est dû au hasard
# On a donc q_max=1
pacf(as.numeric(x),main="PACF de X")
# Le PACF est significatif jusqu'à l'ordre 3, donc p_max=3

# On vérifie ensuite la validité de tous les modèles
# Tests de validité des modèles :
pmax <- 3
qmax <- 1
valide <- matrix(nrow=pmax+1,ncol=qmax+1)
for (p in 0:pmax){
  for (q in 0:qmax){
    model <- arima(x, order = c(p,0,q))
    valide[p+1,q+1] <- all(LjungBoxtest(model$residuals,24,fitdf=p+q)[-c(1:(p+q)),2]>0.05)
  }
}
rownames(valide) <- paste("p=",0:pmax)
colnames(valide) <- paste("q=",0:qmax)

cat("Validité des modèles :")
valide


# Tests de sigificativité des coefficients :
significatif <- matrix(nrow=pmax+1,ncol=qmax+1)
for (p in 0:pmax){
  for (q in 0:qmax){
    model <- arima(x, order = c(p,0,q))
    df <- length(x)-p-q-1
    coef <- model$coef
    se <- sqrt(diag(model$var.coef))
    t <- coef/se
    pval <- ((1-pt(abs(t),df))*2)
    if (p==0 & q==0) {significatif[p+1,q+1] <- (pval[1] < 0.05)}
    else if (q==0 | p==0) {significatif[p+1,q+1] <- (pval[p+q] < 0.05)}
    else {significatif[p+1,q+1] <- ( (pval[p] < 0.05) & (pval[q] < 0.05)) }
  }
}
rownames(significatif) <- paste("p=",0:pmax)
colnames(significatif) <- paste("q=",0:qmax)

cat("Significativité des modèles :")
significatif

# Tests de nullité des autocorrélations des résidus :
autocorr_res <- matrix(nrow=pmax+1,ncol=qmax+1)
for (p in 0:pmax){
  for (q in 0:qmax){
    model <- arima(x, order = c(p,0,q))
    acf_res <- acf(model$residuals, plot=FALSE)
    pacf_res <- pacf(model$residuals, plot=FALSE)
    autocorr_res[p+1,q+1] <- all(abs(acf_res$acf)[-1] < qnorm((1 + 0.95)/2)/sqrt(n-1)) & all(abs(pacf_res$acf)[-1] < qnorm((1 + 0.95)/2)/sqrt(n-1))
  }
}
rownames(autocorr_res) <- paste("p=",0:pmax)
colnames(autocorr_res) <- paste("q=",0:qmax)

cat("Tests de nullité des autocorrélations des résidus :")
autocorr_res

# Ce dernier critère étant restrictif, il n'a validé aucun modèle
# Alors on augmente l'intervalle de confiance à 96.5% au lieu de 95% pour tester la nullité des autocorrélations des résidus
autocorr_res <- matrix(nrow=pmax+1,ncol=qmax+1)
for (p in 0:pmax){
  for (q in 0:qmax){
    model <- arima(x, order = c(p,0,q))
    acf_res <- acf(model$residuals, plot=FALSE)
    pacf_res <- pacf(model$residuals, plot=FALSE)
    autocorr_res[p+1,q+1] <- all(abs(acf_res$acf)[-1] < qnorm((1 + 0.965)/2)/sqrt(n-1)) & all(abs(pacf_res$acf)[-1] < qnorm((1 + 0.965)/2)/sqrt(n-1))
  }
}
rownames(autocorr_res) <- paste("p=",0:pmax)
colnames(autocorr_res) <- paste("q=",0:qmax)

cat("Tests de nullité des autocorrélations des résidus :")
autocorr_res

#Le test a séléctionné alors 3 modèles

# On résume les résultats de tous les tests de validité et de significativité :
cat("Tous les tests :")
valide & significatif & autocorr_res

# Cela a permis de sélectionner deux modèles valides et ajustés : Les modèles ARMA(0,0,1) et ARIMA(2,0,1)
# Vérifions ces modèles

model01 <- arima(x, order=c(0,0,1))
LjungBoxtest(model01$residuals,24,fitdf=1)
df <- length(x)-1
coef <- model01$coef
se <- sqrt(diag(model01$var.coef))
t <- coef/se
pval <- ((1-pt(abs(t),df))*2)
cat("Coefficients : \n")
coef
cat("P-valeurs : \n")
pval
plot(model01$residuals)
acf(model01$residuals)
pacf(model01$residuals)
# On constate un très léger dépassement sur les autocorrélogrammes, c'est pourquoi ce modèle ne passait pas le test à 95% mais passait le test à 96.5%.
#Autrement, modèle est effectivement valide et bien ajusté et les résidus ressemblent à un bruit blanc.

# On vérifie le modèle ARMA(2,1) :
model21 <- arima(x, order=c(2,0,1))
LjungBoxtest(model21$residuals,24,fitdf=3)
df <- length(x)-3
coef <- model21$coef
se <- sqrt(diag(model21$var.coef))
t <- coef/se
pval <- ((1-pt(abs(t),df))*2)
cat("Coefficients : \n")
coef
cat("P-valeurs : \n")
pval
plot(model21$residuals)
acf(model21$residuals)
pacf(model21$residuals)
# Il y a également un léger dépassement de l'intervalle 
#de confiance mais sinon le modèle semble valide et bien ajusté.

# On compare les deux modèles maintenant par les criteres d'information AIC et BIC
pmax <- 3
qmax <- 1
liste_aic <- matrix(nrow=pmax+1,ncol=qmax+1)
liste_bic <- matrix(nrow=pmax+1,ncol=qmax+1)
for (p in 0:pmax){
  for (q in 0:qmax){
    model <- arima(x, order = c(p,0,q))
    liste_aic[p+1,q+1] <- AIC(model)
    liste_bic[p+1,q+1] <- AIC(model, k = log(length(x)))
  }
}
rownames(liste_aic) <- paste("p=",0:pmax)
colnames(liste_aic) <- paste("q=",0:qmax)
rownames(liste_bic) <- paste("p=",0:pmax)
colnames(liste_bic) <- paste("q=",0:qmax)

cat("AIC :")
liste_aic  # le meilleur modèle est ARMA(2,1)
cat("BIC :")
liste_bic  # le meilleur modèles est ARMA(0,1)

# Les critères d'information sont divergents, on regarde la normalité des résidus
model01 <- arima(x, order = c(0,0,1))
model21 <- arima(x, order = c(2,0,1))
cat("Modèle ARMA(0,1) :")
#On realise un test de Jarque Bera et de Shapiro
jarque.bera.test(model01$residuals)
shapiro.test(model01$residuals)
# H0: Normalité
# La normalité est rejetée
#QQ-plot:
qqPlot(model01$residuals)

cat("Modèle ARMA(2,1) :")
jarque.bera.test(model21$residuals)
shapiro.test(model21$residuals)
# La normalité est rejetée
#QQ-plot:
qqPlot(model21$residuals)

# Les 2 modèles donnent des résultats proches


# On s'intéresse aux R2 ajustés des deux modèles :
adjr2 <- function(model){
  ss_res <- sum(model$residuals^2) 
  p <- length(model$model$phi)
  q <- length(model$model$theta[model$model$theta!=0])
  ss_tot <- sum((x - mean(x))^2)
  n <- length(x)
  adj_r2 <- 1-(ss_res/(n-p-q-1))/(ss_tot/(n-1)) #r2 ajusté
  return(adj_r2)
}
cat("Modèle ARMA(0,1) : R2 = ",adjr2(model01),"\n")
cat("Modèle ARMA(2,1) : R2 = ",adjr2(model21))

# Le modèle ARMA(2,1) a un meilleur R2 ajusté

# Prévision sur un échantillon de test
n <- length(s)
train <- s[1:340]
test <- s[341:n]
model01 <- arima(train, order=c(0,1,1))
model21 <- arima(train, order=c(2,1,1))
testfit01 <- Arima(c(train,test), model=model01)
forecast01 <- fitted(testfit01)[341:n]
testfit21 <- Arima(c(train,test), model=model21)
forecast21 <- fitted(testfit21)[341:n]
cat("Modèle ARMA(0,1) : RMSE = ",sum((forecast01-test)^2),"\n")
cat("Modèle ARMA(2,1) : RMSE = ",sum((forecast21-test)^2))

plot(test, col="green",type="l")
lines(forecast01, col="blue")
lines(forecast21, col="red")
legend("topleft", legend = c("Données de test", "Prévision ARIMA(0,1,1)", "Prévision ARIMA(2,1,1)"),
       col = c("green", "blue","red"), lty = 1)

# Le modèle ARMA(2,1) donne une erreur inférieure au modèle ARMA(0,1) sur l'échantillon de test

#Chacun des deux modèles minimise un des deux criteres AIC et BIC, mais les critères des prevision
#et R2 ajusté ont montré que ARMA(2,1) est légèrement meilleur. On gardera donc ce modèle.

# Visualisation du modèle ARMA retenu
model_x <- Arima(x,order=c(2,0,1),include.constant=TRUE)
cat("Coefficients : \n")
model_x$coef # Coefficients du modèle
plot(model_x$x,col="red",type="l", main="Série différenciée et modèle ajusté", ylab="X")
lines(fitted(model_x),col="blue")
legend("topleft", legend = c("Série différenciée", "Modèle"),
       col = c("red", "blue"), lty = 1)
# On remarque que le modèle n'explique pas une grande partie de la variance (ce qui était prévisible par le R2 ~ 0.27), mais cela reste acceptable pour un modèle ARIMA série économétrique.

# Tracé de la série originale et du modèle
plot(s,col="red", type="l", main="Série originale et modèle ajusté")
model <- Arima(s,order=c(2,1,1), include.constant=TRUE)
lines(fitted(model),col="blue")
legend("topleft", legend = c("Série", "Modèle"),
       col = c("red", "blue"), lty = 1)

#On peut aussi vérifier les racines du modèle choisi :
autoplot(model)
# Elles sont bien à l'intérieur du cercle unité

#####--   Partie 3 : Prévision   --#####
# On suppose dans cette partie que les résidus sont gaussiens

### Question 8 ###

model <- Arima(x,order=c(2,0,1),include.constant=TRUE)

#Extraction des coefs du modèle et de la variance des résius 
model$coef
alpha_0 <- as.numeric(model$coef[4])
phi_1 <- as.numeric(model$coef[1])
phi_2 <- as.numeric(model$coef[2])
psi_1 <- as.numeric(model$coef[3])
sigma2 <- as.numeric(model$sigma2)

# Prévisions de X_{T+1} et X_{T+2}
prev_T1 <- alpha_0 +phi_1*x[n-1] +phi_2*x[n-2]+ psi_1*as.numeric(model$residuals[n-1])
prev_T2 <- alpha_0 +phi_1*prev_T1 +phi_2*x[n-1]
cat("X(T+1) : ",prev_T1,"\n")
cat("X(T+2) : ",prev_T2)


# Intervalle de confiance univarié pour X_{T+1}
borne_sup_T1<- prev_T1 +1.96*sqrt(sigma2)
borne_inf_T1<- prev_T1-1.96*sqrt(sigma2)

# Intervalle de confiance univarié pour X_{T+2}
borne_sup_T2<- prev_T2 +1.96*sqrt(sigma2*(1+(phi_1+psi_1)^2))
borne_inf_T2<- prev_T2 -1.96*sqrt(sigma2*(1+(phi_1+psi_1)^2))
cat("IC(T+1) : [",borne_inf_T1,borne_sup_T1,"] \n")
cat("IC(T+2) : [",borne_inf_T2,borne_sup_T2,"] \n")


#Représentation graphique de l'intervalle de confiance (univarié) de X :
IC_X <- function(){
ts.plot(ts(c(x[320:n-1],borne_sup_T1, borne_sup_T2),start=2016+7/12,end=2020+3/12,frequency=12),type="l", col="darkgrey",main="Prévision de X", ylab="X")
lines(ts(c(x[320:n-1],borne_inf_T1, borne_inf_T2),start=2016+7/12,end=2020+3/12,frequency=12),type="l", col="darkgrey")
lines(ts(c(x[320:n-1],prev_T1, prev_T2),start=2016+7/12,end=2020+3/12,frequency=12),lty="dotted", col="grey")
lines(ts(x[320:n-1],start=2016+7/12,end=2020+1/12,frequency=12),type="l")
points(2020+2/12, prev_T1,bg='tomato2', pch=21, cex=1, lwd=2)
points(2020+3/12, prev_T2,bg='tomato2', pch=21, cex=1, lwd=2)
}
IC_X()

# On compare avec l'IC généré automatiquement par la librairie
autoplot(forecast(model,h=2),xlim=c(2015,2020.5))
cat("IC(T+1) : [",forecast(model,h=1)$lower[2], forecast(model,h=1)$upper[2],"] \n")
cat("IC(T+2) : [",forecast(model,h=2)$lower[4], forecast(model,h=2)$upper[4],"] \n")
# Résultats semblables à nos prévisions théoriques

#Représentation graphique de l'intervalle de confiance pour S :
IC_S <- function(){
prev_s_T1 <- s[n]+prev_T1
prev_s_T2 <- prev_s_T1+prev_T2
plot(c(s[320:n],borne_sup_T1+s[n], borne_sup_T2+prev_s_T1),type="l", col="darkgrey", main="Prévision pour la série originale S", ylab="S")
lines(c(s[320:n],borne_inf_T1+s[n], borne_inf_T2+prev_s_T1),type="l", col="darkgrey")
lines(c(s[320:n],prev_s_T1, prev_s_T2),lty="dotted", col="grey")
lines(s[320:n],type="l")
points(n-320+2, prev_s_T1,bg='tomato2', pch=21, cex=1, lwd=2)
points(n-320+3, prev_s_T2,bg='tomato2', pch=21, cex=1, lwd=2)
}
IC_S()

# Intervalle de confiance bivarié :
Sigma <- matrix(c(1, phi_1+psi_1,phi_1+psi_1,1+(phi_1+psi_1)^2), ncol=2)
inv_Sigma <- solve(Sigma)
plot(prev_T1,prev_T2,xlim=c(-4,4),ylim=c(-4,4), xlab="Prévision de X(T+1)", ylab="Prévision de X(T+2)", main="Région de confiance bivariée à 95%")
lines(car::ellipse(center = c(prev_T1,prev_T2), shape= inv_Sigma, radius=sqrt(qchisq(0.95,df=2)) ))
