source("C:/GIT/depth_in_FCA/Neuauflage_12_05_2026/R/functions.R", encoding = 'UTF-8')

library(haven)
setwd("C:/GIT/Datasets/Allbus_2023")
dat <- read_sav("ZA8831_v1-3-0.sav")
#Z <- na.omit(data.frame(links_rechts=dat$pa01,income=as.numeric(dat$incc),age=(dat$age),education=dat$iscd11,equality=dat$im28,need=dat$im27,equity=dat$im22,entitlement=dat$im29))

#Z <- na.omit(data.frame(age=dat$age,spirit=dat$rb14))


#rb07, rb14, ra06
#df <- na.omit(data.frame(age=as.numeric(dat$age),income=as.numeric(dat$incc),education=as.numeric(dat$iscd11),religioesitaet=as.numeric(dat$rb07) ,spiritualitaet=as.numeric(dat$rb14)))


#wichtigkeit
df <- na.omit(data.frame(age=as.numeric(dat$age),income=as.numeric(dat$incc),education=as.numeric(dat$iscd11),li01=as.numeric(dat$li01) ,li02=as.numeric(dat$li02),li03=as.numeric(dat$li03),li04=as.numeric(dat$li04),li05=as.numeric(dat$li05),li07=as.numeric(dat$li07),li08=as.numeric(dat$li08),li09=as.numeric(dat$li09),li10=as.numeric(dat$li10)   ))

#Gruende fuer Kirchenaustritt
#df <- na.omit(data.frame((dat[,c(93,(178:190))])))

#for(k in (4:12)){df[,k] <- df[,k] - mean(df[,k])}
#library(npreg)

#y=df$li02
#M <- gsm(y~df$age,spar=0)
#plot(df$age,y,col="grey")
#points(df$age,fitted.gsm(M),type="p",col="red")


CT <- oofos:::ranking_scaling(df[,-(1:3)],remove.full.columns=FALSE)
#CT <- oofos:::ranking_scaling(df[,-1],remove.full.columns=FALSE)
#CT1 <- oofos:::get_auto_conceptual_scaling(df[,-(1:3)])

#CT2 <- NULL

#for(k in (0:7)){
#CT2 <- cbind(CT2,ranking_scaling(df[,-(1:3)],remove.full.columns = TRUE,threshold=k))
#}

#CT2 <- t(unique(t(CT2)))


#CT <- cbind(CT1,1-CT1,CT2,1-CT2)
dim(CT)
#CT <- t(unique(t(CT)))
#dim(CT)
#CT <- cbind(CT2,1-CT2)

#CT1 <- oofos:::get_auto_conceptual_scaling(cbind(df$religioesitaet,df$spiritualitaet))
#CT2 <- ranking_scaling(cbind(df$religioesitaet,df$spiritualitaet),remove.full.columns = FALSE)
#context <- cbind(CT1,CT2)

ages <-sort(unique(df$age))
I <- rep(0,length(ages))
TT <- sizes <- cbp <- II <- I
condmed <- condmed1 <- condmed2 <- condmed3 <- condmed4 <- condmed5 <- condmed7 <- condmed8 <-condmed9 <- condmed10 <-I
L <- U <- I
condmean1 <- condmean2 <- condmean3 <- condmean4 <- condmean5 <- condmean7 <- condmean8 <- condmean9 <- condmean10 <- condmed1

sortages <- sort(df$age)
Td <- Tv <- list()
for(k in seq_len(length(ages))){
 for( t in (1:100)){
  i=which(  abs(df$age-ages[k])<=t)
  if(length(i)>15){break}
}
  print(k)
  
  sizes[k] <- length(i)
 print(max(rowMeans(CT)))
  Td[[k]] <-  Tukey_depth(matrix(CT[i,],nrow=length(i)))*length(i)
  Tv[[k]] <- Tverberg_depth(CT[i,])$depths
  #plot(Td[[k]],Tv[[k]])
  
  
  #D <- peeling_depth(CT[i,])
  #D <- list(depths=Tukey_depth(CT[i,]))#

  #D <- Tverberg_depth_par(CT[i,])
  
  #temp_tukey <- Tukey_depth(CT[i,])
  j=which(D$depths>=quantile(D$depths,0.8))
  
  #TT[k] <- max(temp_tukey[j])
  condmed[k] <- mean(df$li02[i])

  i <- i[j]
  
 # cbp[k] <- temp$cbp_bound
  
  
  
  condmed1[k] <- median(df$li01[i])
  condmean1[k] <- mean(df$li01[i])
  
  condmed2[k] <- median(df$li02[i])
  condmean2[k] <- mean(df$li02[i])
  
  condmed3[k] <- median(df$li03[i])
  condmean3[k] <- mean(df$li03[i])
  
  condmed4[k] <- median(df$li04[i])
  condmean4[k] <- mean(df$li04[i])
  
  condmed5[k] <- median(df$li05[i])
  condmean5[k] <- mean(df$li05[i])
  
  condmed7[k] <- median(df$li07[i])
  condmean7[k] <- mean(df$li07[i])
  
  condmed8[k] <- median(df$li08[i])
  condmean8[k] <- mean(df$li08[i])
  
  condmed9[k] <- median(df$li09[i])
  condmean9[k] <- mean(df$li09[i])
  
  condmed10[k] <- median(df$li10[i])
  condmean10[k] <- mean(df$li10[i])
  

  #D=DepthProc::depthTukey(as.matrix(df[i,-(1:3)]),as.matrix(df[i,-(1:3)]));j=which.max(D);I[k]=i[j]
 # D <-  Tukeys_depth(matrix(CT[i,],nrow=length(i)));
  #D=peeling_depth(CT[i,]);j=which.max(D);II[k]=i[j]
  j=which(D$depths>=quantile(D$depths,0.85))
  L[k] <- mean(df$li09[i[j]])
  U[k] <- max(df$li09[i[j]])


}

#plot(df$age,df$li07,col="grey")
#lines(ages,condmean1)
#lines(ages,condmean2)
#lines(ages,condmean3)
#lines(ages,condmean4)
#lines(ages,condmean5)
#lines(ages,condmean7)
#lines(ages,condmean8)
#lines(ages,condmean9)
#lines(ages,condmean10)

N <- length(ages)

df_long <-  data.frame( x=rep(ages,9),y=c(condmean1,condmean2,condmean3,condmean4,condmean5,condmean7,condmean8,condmean9,condmean10),
                   area=c(rep(as.factor("own family & children"),N),rep(as.factor("occupation & work"),N),rep(as.factor("leisure & recreation"),N),rep(as.factor("friends & acquaintances"),N),rep(as.factor("relatives"),N),rep(as.factor("politics & public life

"),N),rep(as.factor("neighborhood"),N),rep(as.factor("church"),N),rep(as.factor("religion"),N)))
                   

library(ggplot2)

okabe_ito <- c(
  "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00",
  "#CC79A7", "#000000", "#999999"
)

ggplot(df_long, aes(x, y, colour = area)) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = okabe_ito)

#########
lines(ages,condmed,col="grey")
lines(ages,L,col="red")
lines(ages,U,col="blue")
lines(ages,df$li02[I],col="blue")

i=which(df$age <=5)
M <- gsm(y~age+income+education,data=df)


library(ggplot2)

ggplot(df, aes(x = x, y = y)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", color = "blue", linewidth = 1.2) +
  labs(title = "Conditional Mean (LOESS Smooth)",
       x = "x", y = "E(y | x)") +
  theme_minimal()



library(dplyr)
library(ggplot2)

alpha <- 0.8

df_summary <- df %>%
  mutate(x_bin = cut(x, breaks = 20)) %>%
  group_by(x_bin) %>%
  summarise(
    x_mid = mean(x),
    y_q = quantile(y, probs = alpha),
    .groups = "drop"
  )

ggplot(df_summary, aes(x = x_mid, y = y_q)) +
  geom_line(color = "purple", linewidth = 1) +
  geom_point(color = "purple") +
  labs(title = paste("Binned Conditional", alpha, "Quantile"),
       x = "x", y = "Quantile") +
  theme_minimal()

library(ggplot2)

ggplot(df, aes(x = x, y = y)) +
  geom_density_2d_filled() +
  labs(title = "Conditional Density of y given x",
       x = "x", y = "y") +
  theme_minimal()


dim(Z)


mat <- as.matrix(cbind(Z$equality,Z$need,Z$equity,Z$entitlement))
colnames(mat) <- c("equality","need","equity","entitlement")
context <- ranking.scaling(mat,remove.full.columns = FALSE)
context <- cbind(context,oofos:::get_auto_conceptual_scaling(as.matrix(cbind(Z$equality,Z$need,Z$equity,Z$entitlement))))
indexs <- list()
dev.new()
par(mfrow=c(3,4))

for(k in (1:10)){
  i <- which(Z$links_rechts==k)
  depths <- peeling_depth(context[i,])
  j <- which(depths==max(depths))
  indexs[[k]] <- i[j]
  temp <- context[i[j[1]], (1:16)];dim(temp) <- c(4,4)
  colnames(temp) <- rownames(temp) <- colnames(mat)
  plot(as.relation(temp),main=k)
  print(c(k,median(Z$equality[i])))#print(length(unique(depths))/length(i))
}






tM1 <-lm(as.numeric(equality)~ log(as.numeric(income)), data=Z)
summary(M1)

M2 <-lm(as.numeric(equality)~ as.numeric(age), data=Z)
summary(M2)

M3 <-lm(as.numeric(equality)~as.numeric(age) + log(as.numeric(income)), data=Z)
summary(M3)

Mat <- (-1)* cor(as.matrix(cbind(Z$links_rechts,Z$education,Z$income,Z$equality,Z$need,Z$equity,Z$entitlement)))

M <- gsm(as.numeric(equality)~as.numeric(links_rechts),data=Z)
plot(as.numeric(Z$links_rechts),as.numeric(Z$equality),col="grey")
points(as.numeric(Z$links_rechts),fitted.gsm(M))
