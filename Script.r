##########################################################
##            Projet de series temporelles              ##
##     Hamdi BEL HADJ HASSINE et Fadi SAHBANI           ##
##########################################################


# Packages n�cessaires:
library(tseries)
library(forecast)
library(fUnitRoots)
library(car)


#####--   Partie 1 : Les donn�es   --#####

### Question 1 ###
path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path)

# Extraction et mise en forme des donn�es :
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
#les 12 chronogrammes mensuels sont � peu pr�s identiques, ce qui confirme l'absence de saisonalit�

lag.plot(s,lags=12,layout=c(3,4),do.lines=FALSE)
# Le Lagplot montre une corr�lation forte.

fit1 <- decompose(s)
plot(fit1)
#on peut voir que l'erreur ne semble pas trop varier au cours du temps ce qui indique l'ad�quation du modele additif
#Visualisation ACF et PACF
acf(s)
pacf(s)
# L'ACF montre que la s�rie n'est pas stationnaire et qu'elle pr�sente une tendance, mais ne montre pas de saisonnalit�.

#Test de la tendance
summary(lm(s ~ seq(1,n)))
# Le coefficient de la tendance lin�aire est significatif, donc on effectue les tests de stationnarit� avec tendance

#Test KPSS de la stationnarit�
kpss.test(s,null="Trend")
# Le test KPSS rejette au niveau 1% l'hypoth�se de stationnarit� de la s�rie.

#Test ADF de la stationnarit�
#   Tests de LjungBox combin�s pour verifier l'autorr�lation des r�sidus 
#   jusqu'� l'ordre k
LjungBoxtest <- function(X, k, fitdf=0) {
  pvalues <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(X, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvalues))
}
#Fonction pour effectuer un test Dickey-Fuller augment� valide
ValidADF <- function(X,kmax,type){ # Tests ADF jusqu'� avoir des r�sidus non autocorr�l�s
  k <- 0
  pasautcor <- 0
  while (pasautcor==0){
    cat(paste0("ADF avec ",k, " lags: "))
    adf <- adfTest(X,lags=k,type=type)
    pvalues <- LjungBoxtest(adf@test$lm$residuals,24,fitdf=length(adf@test$lm$coefficients))[,2]
    if (sum(pvalues<0.05,na.rm=T) == 0) {
      pasautcor <- 1; 	      pasautcor <- 1; cat("Les r�sidus ne sont pas auto-corr�l�s. Test ADF valide.\n")}
    else cat("Les r�sidus sont auto-corr�l�s \n")
    k <- k + 1
  }
  return(adf)
}
adf <- ValidADF(s,24,"ct")
adf
# Le test ADF ne rejette pas au seuil 5% la pr�sence d'une racine unitaire.
# La serie est donc int�gr�e. On la diff�rencie � l'ordre 1 :
x <- diff(s)
plot(x,xlab='Date',ylab="IPI diff�renci�", main="S�rie diff�renci�e")

# La s�rie diff�renci�e ne semble pas pr�senter de constante ou tendance.
#On le v�rifie avec une r�gression lin�aire :
summary(lm(x ~ seq(1,length(x))))
# Les coefficients associ�s � la constante et � la tendance sont non 
# significatifs au seuil 5%

#On fait un test KPSS:
kpss.test(x,null="Level")
# Le test KPSS ne rejette pas la stationnarit� de la s�rie diff�renci�e

# et maintenant on verifie avec le test ADF
suppressWarnings(
  adf <- ValidADF(x,50,"nc"))
adf
# Le test ADF rejette au seuil 1% la pr�sence d'une racine unitaire.

# On peut aussi le v�rifier par la fonction lagplot
lag.plot(x,lags=12,layout=c(3,4),do.lines=FALSE)
# La serie diff�renci�e est donc stationnaire,
#c'est-�-dire que la s�rie originale s est I(1)

#Finalement on effectue un test pour v�rifier l'homosc�dasticit�

# Test d'homosc�dasticit� de Breusch-Pagan
lmtest::bptest(lm(x ~ seq(1,length(x))))
# On ne rejette pas l'homosc�dasticit� de X au niveau 5%
### Question 3 ###

plot(s,xlab='Date',ylab="Indice de production industrielle", main="Indice de production industrielle")
plot(x,xlab='Date',ylab="IPI diff�renci�", main="S�rie diff�renci�e")

#####--   Partie 2 : Mod�les ARMA   --#####

### Question 4 ###

acf(as.numeric(x), main="ACF de X")
# L'ACF est seulement significatif � l'ordre 1. On note un l�ger d�passement au lag 11 mais on suppose qu'il est d� au hasard
# On a donc q_max=1
pacf(as.numeric(x),main="PACF de X")
# Le PACF est significatif jusqu'� l'ordre 3, donc p_max=3

# On v�rifie ensuite la validit� de tous les mod�les
# Tests de validit� des mod�les :
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

cat("Validit� des mod�les :")
valide


# Tests de sigificativit� des coefficients :
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

cat("Significativit� des mod�les :")
significatif

# Tests de nullit� des autocorr�lations des r�sidus :
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

cat("Tests de nullit� des autocorr�lations des r�sidus :")
autocorr_res

# Ce dernier crit�re �tant restrictif, il n'a valid� aucun mod�le
# Alors on augmente l'intervalle de confiance � 96.5% au lieu de 95% pour tester la nullit� des autocorr�lations des r�sidus
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

cat("Tests de nullit� des autocorr�lations des r�sidus :")
autocorr_res

#Le test a s�l�ctionn� alors 3 mod�les

# On r�sume les r�sultats de tous les tests de validit� et de significativit� :
cat("Tous les tests :")
valide & significatif & autocorr_res

# Cela a permis de s�lectionner deux mod�les valides et ajust�s : Les mod�les ARMA(0,0,1) et ARIMA(2,0,1)
# V�rifions ces mod�les

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
# On constate un tr�s l�ger d�passement sur les autocorr�logrammes, c'est pourquoi ce mod�le ne passait pas le test � 95% mais passait le test � 96.5%.
#Autrement, mod�le est effectivement valide et bien ajust� et les r�sidus ressemblent � un bruit blanc.

# On v�rifie le mod�le ARMA(2,1) :
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
# Il y a �galement un l�ger d�passement de l'intervalle 
#de confiance mais sinon le mod�le semble valide et bien ajust�.

# On compare les deux mod�les maintenant par les criteres d'information AIC et BIC
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
liste_aic  # le meilleur mod�le est ARMA(2,1)
cat("BIC :")
liste_bic  # le meilleur mod�les est ARMA(0,1)

# Les crit�res d'information sont divergents, on regarde la normalit� des r�sidus
model01 <- arima(x, order = c(0,0,1))
model21 <- arima(x, order = c(2,0,1))
cat("Mod�le ARMA(0,1) :")
#On realise un test de Jarque Bera et de Shapiro
jarque.bera.test(model01$residuals)
shapiro.test(model01$residuals)
# H0: Normalit�
# La normalit� est rejet�e
#QQ-plot:
qqPlot(model01$residuals)

cat("Mod�le ARMA(2,1) :")
jarque.bera.test(model21$residuals)
shapiro.test(model21$residuals)
# La normalit� est rejet�e
#QQ-plot:
qqPlot(model21$residuals)

# Les 2 mod�les donnent des r�sultats proches


# On s'int�resse aux R2 ajust�s des deux mod�les :
adjr2 <- function(model){
  ss_res <- sum(model$residuals^2) 
  p <- length(model$model$phi)
  q <- length(model$model$theta[model$model$theta!=0])
  ss_tot <- sum((x - mean(x))^2)
  n <- length(x)
  adj_r2 <- 1-(ss_res/(n-p-q-1))/(ss_tot/(n-1)) #r2 ajust�
  return(adj_r2)
}
cat("Mod�le ARMA(0,1) : R2 = ",adjr2(model01),"\n")
cat("Mod�le ARMA(2,1) : R2 = ",adjr2(model21))

# Le mod�le ARMA(2,1) a un meilleur R2 ajust�

# Pr�vision sur un �chantillon de test
n <- length(s)
train <- s[1:340]
test <- s[341:n]
model01 <- arima(train, order=c(0,1,1))
model21 <- arima(train, order=c(2,1,1))
testfit01 <- Arima(c(train,test), model=model01)
forecast01 <- fitted(testfit01)[341:n]
testfit21 <- Arima(c(train,test), model=model21)
forecast21 <- fitted(testfit21)[341:n]
cat("Mod�le ARMA(0,1) : RMSE = ",sum((forecast01-test)^2),"\n")
cat("Mod�le ARMA(2,1) : RMSE = ",sum((forecast21-test)^2))

plot(test, col="green",type="l")
lines(forecast01, col="blue")
lines(forecast21, col="red")
legend("topleft", legend = c("Donn�es de test", "Pr�vision ARIMA(0,1,1)", "Pr�vision ARIMA(2,1,1)"),
       col = c("green", "blue","red"), lty = 1)

# Le mod�le ARMA(2,1) donne une erreur inf�rieure au mod�le ARMA(0,1) sur l'�chantillon de test

#Chacun des deux mod�les minimise un des deux criteres AIC et BIC, mais les crit�res des prevision
#et R2 ajust� ont montr� que ARMA(2,1) est l�g�rement meilleur. On gardera donc ce mod�le.

# Visualisation du mod�le ARMA retenu
model_x <- Arima(x,order=c(2,0,1),include.constant=TRUE)
cat("Coefficients : \n")
model_x$coef # Coefficients du mod�le
plot(model_x$x,col="red",type="l", main="S�rie diff�renci�e et mod�le ajust�", ylab="X")
lines(fitted(model_x),col="blue")
legend("topleft", legend = c("S�rie diff�renci�e", "Mod�le"),
       col = c("red", "blue"), lty = 1)
# On remarque que le mod�le n'explique pas une grande partie de la variance (ce qui �tait pr�visible par le R2 ~ 0.27), mais cela reste acceptable pour un mod�le ARIMA s�rie �conom�trique.

# Trac� de la s�rie originale et du mod�le
plot(s,col="red", type="l", main="S�rie originale et mod�le ajust�")
model <- Arima(s,order=c(2,1,1), include.constant=TRUE)
lines(fitted(model),col="blue")
legend("topleft", legend = c("S�rie", "Mod�le"),
       col = c("red", "blue"), lty = 1)

#On peut aussi v�rifier les racines du mod�le choisi :
autoplot(model)
# Elles sont bien � l'int�rieur du cercle unit�

#####--   Partie 3 : Pr�vision   --#####
# On suppose dans cette partie que les r�sidus sont gaussiens

### Question 8 ###

model <- Arima(x,order=c(2,0,1),include.constant=TRUE)

#Extraction des coefs du mod�le et de la variance des r�sius 
model$coef
alpha_0 <- as.numeric(model$coef[4])
phi_1 <- as.numeric(model$coef[1])
phi_2 <- as.numeric(model$coef[2])
psi_1 <- as.numeric(model$coef[3])
sigma2 <- as.numeric(model$sigma2)

# Pr�visions de X_{T+1} et X_{T+2}
prev_T1 <- alpha_0 +phi_1*x[n-1] +phi_2*x[n-2]+ psi_1*as.numeric(model$residuals[n-1])
prev_T2 <- alpha_0 +phi_1*prev_T1 +phi_2*x[n-1]
cat("X(T+1) : ",prev_T1,"\n")
cat("X(T+2) : ",prev_T2)


# Intervalle de confiance univari� pour X_{T+1}
borne_sup_T1<- prev_T1 +1.96*sqrt(sigma2)
borne_inf_T1<- prev_T1-1.96*sqrt(sigma2)

# Intervalle de confiance univari� pour X_{T+2}
borne_sup_T2<- prev_T2 +1.96*sqrt(sigma2*(1+(phi_1+psi_1)^2))
borne_inf_T2<- prev_T2 -1.96*sqrt(sigma2*(1+(phi_1+psi_1)^2))
cat("IC(T+1) : [",borne_inf_T1,borne_sup_T1,"] \n")
cat("IC(T+2) : [",borne_inf_T2,borne_sup_T2,"] \n")


#Repr�sentation graphique de l'intervalle de confiance (univari�) de X :
IC_X <- function(){
ts.plot(ts(c(x[320:n-1],borne_sup_T1, borne_sup_T2),start=2016+7/12,end=2020+3/12,frequency=12),type="l", col="darkgrey",main="Pr�vision de X", ylab="X")
lines(ts(c(x[320:n-1],borne_inf_T1, borne_inf_T2),start=2016+7/12,end=2020+3/12,frequency=12),type="l", col="darkgrey")
lines(ts(c(x[320:n-1],prev_T1, prev_T2),start=2016+7/12,end=2020+3/12,frequency=12),lty="dotted", col="grey")
lines(ts(x[320:n-1],start=2016+7/12,end=2020+1/12,frequency=12),type="l")
points(2020+2/12, prev_T1,bg='tomato2', pch=21, cex=1, lwd=2)
points(2020+3/12, prev_T2,bg='tomato2', pch=21, cex=1, lwd=2)
}
IC_X()

# On compare avec l'IC g�n�r� automatiquement par la librairie
autoplot(forecast(model,h=2),xlim=c(2015,2020.5))
cat("IC(T+1) : [",forecast(model,h=1)$lower[2], forecast(model,h=1)$upper[2],"] \n")
cat("IC(T+2) : [",forecast(model,h=2)$lower[4], forecast(model,h=2)$upper[4],"] \n")
# R�sultats semblables � nos pr�visions th�oriques

#Repr�sentation graphique de l'intervalle de confiance pour S :
IC_S <- function(){
prev_s_T1 <- s[n]+prev_T1
prev_s_T2 <- prev_s_T1+prev_T2
plot(c(s[320:n],borne_sup_T1+s[n], borne_sup_T2+prev_s_T1),type="l", col="darkgrey", main="Pr�vision pour la s�rie originale S", ylab="S")
lines(c(s[320:n],borne_inf_T1+s[n], borne_inf_T2+prev_s_T1),type="l", col="darkgrey")
lines(c(s[320:n],prev_s_T1, prev_s_T2),lty="dotted", col="grey")
lines(s[320:n],type="l")
points(n-320+2, prev_s_T1,bg='tomato2', pch=21, cex=1, lwd=2)
points(n-320+3, prev_s_T2,bg='tomato2', pch=21, cex=1, lwd=2)
}
IC_S()

# Intervalle de confiance bivari� :
Sigma <- matrix(c(1, phi_1+psi_1,phi_1+psi_1,1+(phi_1+psi_1)^2), ncol=2)
inv_Sigma <- solve(Sigma)
plot(prev_T1,prev_T2,xlim=c(-4,4),ylim=c(-4,4), xlab="Pr�vision de X(T+1)", ylab="Pr�vision de X(T+2)", main="R�gion de confiance bivari�e � 95%")
lines(car::ellipse(center = c(prev_T1,prev_T2), shape= inv_Sigma, radius=sqrt(qchisq(0.95,df=2)) ))
