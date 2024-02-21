#Neil Davies 08/11/22
#This generates the data for the advanced MR course

set.seed(12456)

data<-data.frame(matrix(NA,nrow = 40000,ncol = 1))

#Create ID variable
data$id<-seq.int(nrow(data)) 

#Drop the first column (suggestions welcome for how to get rid of this line)
data<-subset(data,select="id")

#Draw parental genotypes from the binomial distribution 
data$m1 <- rbinom(40000, 1, 0.6)
data$m2 <- rbinom(40000, 1, 0.6)
data$ma<-data$m1+data$m2
data$f1 <- rbinom(40000, 1, 0.6)
data$f2 <- rbinom(40000, 1, 0.6)
data$fa<-data$f1+data$f2

#Draw inherited allele for the each of two children:
data$mia_sib1 <- rbinom(40000, 1, 0.5)
data$fia_sib1 <- rbinom(40000, 1, 0.5)

data$fia_sib2 <- rbinom(40000, 1, 0.5)
data$mia_sib2 <- rbinom(40000, 1, 0.5)

#Define sib genotype, offspring 1 and 2 defined o1 o2
#Offspring 1
data$o11 <-data$m1*data$mia_sib1+data$m2*(1-data$mia_sib1)
data$o12 <-data$f1*data$fia_sib1+data$f2*(1-data$fia_sib1)
data$o1a<-data$o11+data$o12
#Offspring 2
data$o21 <-data$m1*data$mia_sib2+data$m2*(1-data$mia_sib2)
data$o22 <-data$f1*data$fia_sib2+data$f2*(1-data$fia_sib2)
data$o2a<-data$o21+data$o22

#Thus, we've got parents, and two sibs
#Generate error terms for the exposure and the outcomes
data$o1e_error<-rnorm(40000,0,1)
data$o2e_error<-rnorm(40000,0,1)
data$o1o_error<-rnorm(40000,0,1)
data$o2o_error<-rnorm(40000,0,1)
data$me_error<-rnorm(40000,0,1)
data$fe_error<-rnorm(40000,0,1)

#Generate parental BMI
data$m_bmi<-data$ma+data$me_error
data$f_bmi<-data$fa+data$fe_error

#Generate exposure (BMI)
data$o1_bmi<-data$o1a+data$o1e_error
data$o2_bmi<-data$o2a+data$o2e_error

#Generate outcome
#This is a function, not of BMI, but of parents' BMI
data$o1_educ<-data$m_bmi+data$f_bmi+data$o1o_error
data$o2_educ<-data$m_bmi+data$f_bmi+data$o2o_error

save(data, file = "family_sims.RData")