####################################################################
##########  Survival analysis using Kaplan Meier Estimator  ########
####################################################################

######## Installing 'Survival' package
install.packages("survival")

######## Loading 'Survival' package
library(survival)

######## Opening lung dataset
View(lung)

#############################################################################################
####  Estimating Kaplan Meier survival probabilities for males and females step-by-step  ####
#############################################################################################

######## For males

######## Subsetting lung dataset for males
male=lung[lung$sex==1,]

######## Stratification of male dataset on the basis of 'status' column
######## A status value of 1 indicates that the patient was censored
######## A status value of 2 indicates that the patient died
######## 'km_male' dataframe contains total male patients (both censored & dead) at each time interval
km_male=aggregate(male$status,list(male$time),length)

######## 'm' dataframe contains total censored and total dead male patients at each time interval
m=data.frame(tapply(male$status,list(male$time,male$status),length))

######## Renaming columns in 'km_male' dataset with appropriate titles
colnames(km_male)[1] = "Time_Interval"
colnames(km_male)[2] = "Dead_Censored"

######## Assigning columns from 'm' data frame to 'km_male' data frame
km_male$male_dead=m$X2
km_male$male_censored=m$X1

######## Updating all the 'NA' values to 0 in km_male dataframe
km_male[is.na(km_male)] = 0

######## Calculating number of males alive (total - (censored + dead)) at each time interval
for(i in 1:nrow(km_male))
{
  km_male$male_live[i] = ifelse(i==1,nrow(male),km_male$male_live[i-1]-km_male$male_dead[i-1]-km_male$male_censored[i-1])
}

######## Calculating survival probability for males at each time interval
km_male$survival=1-(km_male$male_dead/km_male$male_live)

######## Taking the product survival probabilities at each time interval 
######## to obtain the total survival probability for all males
for(i in 1:nrow(km_male))
{
  km_male$probability[i] = ifelse(i==1,km_male$survival[i],km_male$survival[i]*km_male$probability[i-1])
}

######## Calculating standard error for the survival probability for males
for(i in 1:nrow(km_male))
{
  merror = km_male$male_dead/(km_male$male_live*(km_male$male_live-km_male$male_dead))
  km_male$standard_error[i] = round(km_male$probability[i]*(sqrt(sum(merror[1:i]))),4)
}

######## Generating 95% confidence intervals for the survival probabilities for males
km_male$Lower_CI =round(km_male$probability-1.96*km_male$standard_error,4)
km_male$Upper_CI = round(ifelse(km_male$probability+1.96*km_male$standard_error>1,1,km_male$probability+1.96*km_male$standard_error),4)


######## For females

######## Subsetting lung dataset for females
female=lung[lung$sex==2,]

######## Stratification of female dataset on the basis of 'status' column
######## A status value of 1 indicates that the patient was censored
######## A status value of 2 indicates that the patient died
######## 'km_male' dataframe contains total female patients (both censored & dead) at each time interval
km_female=aggregate(female$status,list(female$time),length)

######## 'f' dataframe contains total censored and total dead female patients at each time interval
f=data.frame(tapply(female$status,list(female$time,female$status),length))

######## Renaming columns in 'km_female' dataset with appropriate titles
colnames(km_female)[1] = "Time_Interval"
colnames(km_female)[2] = "Dead_Censored"

######## Assigning columns from 'f' data frame to 'km_female' data frame
km_female$female_dead=f$X2
km_female$female_censored=f$X1

######## Updating all the 'NA' values to 0 in 'km_female' dataframe
km_female[is.na(km_female)] = 0

######## Calculating number of females alive (total - (censored + dead)) at each time interval
for(i in 1:nrow(km_female))
{
  km_female$female_live[i] = ifelse(i==1,nrow(female),km_female$female_live[i-1]-km_female$female_dead[i-1]-km_female$female_censored[i-1])
}

######## Calculating survival probability for females at each time interval
km_female$survival=1-(km_female$female_dead/km_female$female_live)

######## Taking the product survival probabilities at each time interval 
######## to obtain the total survival probability for all females
for(i in 1:nrow(km_female))
{
  km_female$probability[i] = ifelse(i==1,km_female$survival[i],km_female$survival[i]*km_female$probability[i-1])
}

######## Calculating standard error for the survival probability for females
for(i in 1:nrow(km_female))
{
  ferror = km_female$female_dead/(km_female$female_live*(km_female$female_live-km_female$female_dead))
  km_female$standard_error[i] = round(km_female$probability[i]*(sqrt(sum(ferror[1:i]))),4)
}

######## Generating 95% confidence intervals for the survival probabilities for females
km_female$Lower_CI =round(km_female$probability-1.96*km_female$standard_error,4)
km_female$Upper_CI = round(ifelse(km_female$probability+1.96*km_female$standard_error>1,1,km_female$probability+1.96*km_female$standard_error),4)


######## Plotting survival curves and their confidence intervals generated above for males and females
plot(km_male$Time_Interval,km_male$probability, col = "blue",type = "s", lwd=4,xlab="Time_Elapsed", ylab="Survival_Probability", main="Kaplan meier-Survival curves")
lines(km_female$Time_Interval,km_female$probability, col = "red",type = "s", lwd=4)
legend(750,1,legend=c("Female","Male"),lty=c(1,1),col=c("red", "blue")) 
cm = km_male[km_male$male_censored!=0,]
cf = km_female[km_female$female_censored!=0,]
points(cm$Time_Interval,cm$probability, pch=3, lwd=2)
points(cf$Time_Interval,cf$probability, pch=3, lwd=2)
lines(km_male$Time_Interval, km_male$Lower_CI, type="s", pch=19, cex=0.1, col="blue")
lines(km_male$Time_Interval, km_male$Upper_CI, type="s", pch=19, cex=0.1, col="blue")
lines(km_female$Time_Interval, km_female$Lower_CI, type="s", pch=19, cex=0.1, col="red")
lines(km_female$Time_Interval, km_female$Upper_CI, type="s", pch=19, cex=0.1, col="red")


##########################################################################################
####  Estimating survival probabilities for males and females using survfit function  ####
##########################################################################################

######## Creating a survival object 'survobj'
lung$survobj= Surv(lung$time,lung$status)

######## Generating survival curves for males and females using survfit function
kmfitl=survfit(lung$survobj~lung$sex)
kmfit=survfit(lung$survobj~1)
######## Viewing the summary of survfit function that gives the chi-squared value and p-value
summary(kmfitl)
summary(kmfit)
######## Plotting the survival curves for male and female patients
######## Survival curves for females is plotted in red
######## Survival curves for males is plotted in blue
plot(kmfitl, col=c("blue","red"))


###################################################################
########  Using Log rank test to test the null hypothesis  ########
###################################################################

######## Null hypothesis: The probability of survival in male and female patients with lung cancer is the same


######## Creating a 'gender' column from the values in column 'sex'
######## A value of 1 is Male
######## A value of 2 is Female
lung$gender=ifelse(lung$sex==1,"male","female")

######## Stratification of lung dataset on the basis of 'status' and 'gender' columns
######## A status value of 1 indicates that the patient was censored
######## A status value of 2 indicates that the patient died
######## 'km' dataframe contains total patients (both censored & dead) at each time interval
km=aggregate(lung$status,list(lung$time),length)
######## 'ml' dataframe contains total patients (both censored & dead) at each time interval in each group (male and female)
ml=data.frame(tapply(lung$status,list(lung$time,lung$status,lung$gender),length))

######## Renaming columns in 'km' dataset with appropriate titles
colnames(km)[1] = "Time_Interval"
colnames(km)[2] = "Dead_Censored"

######## Assigning columns from 'ml' data frame to 'km' data frame
km$male_dead=ml$X2.male
km$female_dead=ml$X2.female
km$male_censored=ml$X1.male
km$female_censored=ml$X1.female

######## Updating all the 'NA' values to 0 in km dataframe
km[is.na(km)] = 0

######## Adding 'Total_censored' column which is the sum of male and female patients censored at each time interval
km$Total_censored = km$male_censored+km$female_censored

######## Calculating number of males and femles alive (total - (censored + dead)) at each time interval
for(i in 1:nrow(km))
{
  km$male_live[i] = ifelse(i==1,nrow(lung[lung$sex==1,]),km$male_live[i-1]-km$male_dead[i-1]-km$male_censored[i-1])
  km$female_live[i] = ifelse(i==1,nrow(lung[lung$sex==2,]),km$female_live[i-1]-km$female_dead[i-1]-km$female_censored[i-1])
}

######## Calculating total people alive at the beginning of each time interval 
######## Number of males alive + Number of females alive
km$Total_live=km$male_live + km$female_live

######## Calculating total people that died during each time interval
km$Total_dead=abs(km$Dead_Censored-km$Total_censored)


######## Calculating chi-squared value as part of the Long rank test
######## chi-squared value = (observed value - expected value)^2/expected value

######## Calculating expected value for males and females at each time interval
######## Expected value for males is the number of males expected to be dead for a given time interval
######## Expected value for females is the number of females expected to be dead for a given time interval
km$male_ex =round(((km$Total_dead*km$male_live)/km$Total_live),2)
km$female_ex =round(((km$Total_dead*km$female_live)/km$Total_live),2)

om = sum(km$male_dead)    ######## 'om' is the observed value for males
em = sum(km$male_ex)      ######## 'em' is the expected value for males
of = sum(km$female_dead)  ######## 'of' is the observed value for females
ef = sum(km$female_ex)    ######## 'ef' is the expected value for females
chi_square = ((om-em)^2/em)+((of-ef)^2/ef)
chi_square

######## Once we obtain the chi-square value, we looked up the p-value from chi-squared distribution table for given degrees of freedom
####As p-value(0.00131)< 0.05, we reject the null hypothesis that the probability of survival in male and female patients with lung cancer is the same

