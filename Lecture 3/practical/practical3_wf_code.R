#Neil Davies 08/11/22
#Code for Practical 3 for the Advanced MR short course - within family MR

# 0. Clear workspace
rm(list = ls())

# 1. Load libraries

#install.packages("AER")
#install.packages("plm")
#install.packages("multiwayvcov")
#install.packages("ivpack")
#install.packages("data.table") 

library("AER")
library('data.table')
library('plm')
#library('ivpack')
library('dplyr')
library('multiwayvcov')
library('lmtest')
library('sandwich')
#library('ivreg')

#Load the data, see generate_sim.R for data generating process. 

load("family_sims.RData")

#Part 1. 

#1. Summarise the data
head(data)

#1.a
#Maternal genotypes
table(data$m1)
table(data$m2)
table(data$ma)

#paternal genotypes
table(data$f1)
table(data$f2)
table(data$fa)

#Offspring 1 alleles\
table(data$o11)
table(data$o12)
table(data$o1a)

#Offspring 2 alleles\
table(data$o21)
table(data$o22)
table(data$o2a)


#One genotype with two alleles, that have been coded in additive format. 

#2.
#We can test whether the inheritance of genetic variants within families is random, by estimating the 

#Correlation between genotypes of mother and father of offspring1

#The maternally inherited variant
summary(lm(formula = o11~m1+m2+f1+f2, data = data))
#The paternally inherited variant
summary(lm(formula = o12~m1+m2+f1+f2, data = data))

#3.
#Correlation between genotypes of mother and father of offspring2
#The maternally inherited variant
summary(lm(formula = o21~m1+m2+f1+f2, data = data))
#The paternally inherited variant
summary(lm(formula = o22~m1+m2+f1+f2, data = data))


#4. 
#Correlation between genotypes of offspring1 and offspring2
summary(lm(formula = o2a~o1a, data = data))

#5. 
#First regression of phenotypic education on bmi 
summary(lm(formula = o1_educ~o1_bmi, data = data))


#6. 
#Note there is no causal relationship in our simulated code, but we get a strong positive association via parents alleles.
#If we control for parents genotype
summary(lm(formula = o1_educ~o1_bmi+ma+fa, data = data))

#Everything is fixed!
#Note: all the confounding here is through parental genotype. If we also had confounding at the offspring level then this estimate would also be biased.


#7. 
#Classic MR on the first sib
#Use the offspring's genotype as a instrument for their BMI
classic_mr<-ivreg(formula = o1_educ~o1_bmi|o1a, data = data)
summary(classic_mr)

#Biased, b=0.98569, se=0.01506, p<2e-16. 
#We know from the data generating process, that the true coefficient on the offspring BMI is zero. 


#8. 
#What happens when we adjust for father and mothers' genotype?
wf_mr<-ivreg(formula = o1_educ~o1_bmi+ma+fa|o1a+ma+fa, data = data)
summary(wf_mr)

#Everything works! Unbiased, b=0.01194, se=0.01732, p=0.491, no evidence of an effect. 
#The WF estimate is slightly less precise. 

##Part 2.

#However, above we only used one of the two siblings. 
#How can we include both siblings, as this would increase power? 

#1
#We can reshape into long format:
data_long1=subset(data, select=-c(fia_sib2,mia_sib2,o21,o22, o2a, o2e_error,o2o_error, o2_bmi, o2_educ))
data_long1$sibling=1

colnames(data_long1)[8] ="mia_sib"
colnames(data_long1)[9] ="fia_sib"
colnames(data_long1)[10] ="o1"
colnames(data_long1)[11] ="o2"
colnames(data_long1)[12] ="oa"
colnames(data_long1)[13] ="oe_error"
colnames(data_long1)[14] ="oo_error"
colnames(data_long1)[19] ="o_bmi"
colnames(data_long1)[20] ="o_educ"

data_long2=subset(data, select=-c(fia_sib1,mia_sib1,o11,o12, o1a, o1e_error,o1o_error, o1_bmi, o1_educ))
data_long2$sibling=2

colnames(data_long2)[8] ="fia_sib"
colnames(data_long2)[9] ="mia_sib"
colnames(data_long2)[10] ="o1"
colnames(data_long2)[11] ="o2"
colnames(data_long2)[12] ="oa"
colnames(data_long2)[13] ="oe_error"
colnames(data_long2)[14] ="oo_error"
colnames(data_long2)[19] ="o_bmi"
colnames(data_long2)[20] ="o_educ"

data_long=rbind(data_long1,data_long2)

#We can run the same regression using both siblings:

#We now need to use the plm (panel data estimators) package, which allows for family fixed effects and cluster robust standard errors
#In plm we specify the group variable for family, and use the within estimator. 

#First lets check if we use two siblings and cluster across families WITHOUT adjusting for parental genotype 
#we get the same biased estimator as above for the ivreg package. 

#2. 
#Classic MR, without adjustment for parental genotype
estimates<-plm(o_educ~o_bmi|oa, data = data_long, index = "id", model = "pooling", inst.method = "bvk")
summary(estimates)

#3. 
estimates_robust<-coeftest(estimates,vcov=vcovHC(estimates,type="HC0",cluster="group"))
estimates_robust

#Without adjusting the standard errors, the classic MR is still biased, b=0.993815, se=0.010744, p<2e-16. 
#After adjustment for the clustering within families the SE increases slightly, but the estimates are still biased.
# b=0.993815, se=0.012190. 

#4.
#Within family MR: 
estimates_wf<-plm(o_educ~o_bmi+ma+fa|oa+ma+fa, data = data_long, index = "id", model = "pooling", inst.method = "bvk")
summary(estimates_wf)

#5. 
estimates_robust_wf<-coeftest(estimates_wf,vcov=vcovHC(estimates,type="HC0",cluster="group"))
estimates_robust_wf

#After adjustment for parental genotype, but without clustering the standard errors by family the estimates are now unbiased, but:
#Within family MR b=-0.0061713, se=0.0123816, p=0.6182
#After clustering the standard errors, the beta is still unbiased, but the SE increases. 
#Within family MR, clustered SEs b=-0.0061713 , se=0.0121903, p=0.6127

#In many circumstances, studies will not have sampled both parents and offspring. 
#An alternative is to use siblings. The difference in genotypes between siblings is random. 
#We can exploit this fact to control for familial or demographic effects. 

#There are a number of ways to implement this estimator. 
#Here we will run two, the within "fixed effects" estimator, and the within transformation. 
#They achieve the same ends, the within transformation requires slightly more code, but is more computationally efficient. 

#7.
#We can call the fixed effects estimator using the plm package:
estimates_fe<-plm(o_educ~o_bmi|oa, data = data_long, index = "id", model = "within", inst.method = "bvk")
summary(estimates_fe)

#Note - all we have done is changed the model to "within", which calls the fixed effects estimator, and 
#we now, DO NOT adjust for parental genotype. 

#8. 
estimates_robust_fe<-coeftest(estimates_fe,vcov=vcovHC(estimates,type="HC0",cluster="group"))
estimates_robust_fe

#Without adjusting the standard errors, the fixed effects WF MR is unbiased, b=-0.0076306, se=0.0100510, p=0.4477. 
#After adjustment for the clustering within families the SE increases slightly, but the estimates are still biased.
# b=-0.0076306, se=0.0121903. 


##Part 3.

#A second method is to use the within transformation.

#1. 
#To do this we need to calculate and adjust for the mean genotype. 
# Group by mean of multiple columns
data_long2 <- data_long %>% group_by(id) %>% 
  summarise(mean_oa=mean(oa),
            .groups = 'drop') %>%
  as.data.frame()
#Merge this back onto the main dataframe
data_long3 <- merge(data_long, data_long2, by = c("id"))   

#2. 
#Run the standard MR regressions as above, but now adjusting for the mean sibling genotype
estimates_within<-plm(o_educ~o_bmi+mean_oa|oa+mean_oa, data = data_long3, index = "id", model = "pooling", inst.method = "bvk")
estimates_robust_within<-coeftest(estimates_within,parm="o_bmi", vcov=vcovHC(estimates,type="HC0",cluster="group"))

#3. 
#Summarise IV estimates for fixed effects WF MR
summary(estimates_within)
estimates_robust_within

#Without adjusting the standard errors, the within transformation WF MR is unbiased, b=-0.0076306, se=0.0183318, p=0.6772. 
#After adjustment for the clustering within families the SE decreases slightly, but the estimates are still biased.
# b=-0.0076306, se=0.0121903. 

#Either using the fixed effects estimator, or using the within transformation across families of siblings, 
#can be used to control for familial or demographic biases. 
