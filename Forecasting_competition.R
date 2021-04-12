install.packages("fpp","fpp2","fpp3")
install.packages("gtrendsR")
install.packages("tidyverse","tidyquant","timetk","sweep","forecast")
install.packages("ggplot2")
library(fpp)
library(fpp2)
library(fpp3)
library(gtrendsR)
library(tidyverse)
library(tidyquant)
library(timetk)
library(sweep)
library(forecast)
library(ggplot2)
library(h2o)
library(MARSS)
library(urca)
#library(lubridate)

install.packages('prophet')
library(prophet)
library(Rcpp)
library(rlang)
library(MLmetrics)


####################################################
######                                        ######
######                                        ######
######                 PREAMBULE              ######
######             declare functions          ######
######                                        ######
####################################################



Loss <- function(err, type = 'MSFE')
{
  switch(type,
         'MSFE' = colMeans(err^2), #by default
         'RMSFE' = sqrt(colMeans(err^2)),
         'MAPE' = colMeans(abs(err))) 
}


####################################################
######                                        ######
######                  Start                 ######
######                    &                   ######
######               load dataset             ######
######                                        ######
####################################################


rm() # delete all the existing variables


## DOWNLOADING DATA  -- weekly Google searchs for "car" in the US
#Gtrend <- gtrends(c("car"),geo=c("US"),gprop ="web", time = "today+5-y")
#data <- Gtrend$interest_over_time
data <- read_csv('data3_week9.csv')
print(head(data))

car_weekly <- data %>%
  mutate(Week = yearweek(date)) %>%
  select(hits,date,Week) %>%
  as_tsibble(index = Week) 

autoplot(car_weekly, .vars=hits) + geom_line(aes(group=1), color = "steelblue")
head(car_weekly)  #the Weak variable will progressively move through time
print(tail(car_weekly))


# GENERATE A TIME SERIES OF HITS -- if you don't want to use the tsibble at first
attach(data)
Y <- ts(hits, frequency = 52, start = c(lubridate::year(car_weekly$date)[1],lubridate::week(car_weekly$date)[1]))
plot(Y)


####################################################
######                                        ######
######              cleaning function         ######
######                                        ######
####################################################




clean_time_series <- function(Y, method='None', start_index=221, end_index=235){ # Find "better" ways to clean the Series -- implement other methods ! :)
  Y_cleaned = Y 
  # Changing the dataset might not be the best way to deal with such a problem as we corrupt the signal
  
  if (method == 'None'){
    # do nothing
  }
  if (method == 'Naive'){
    # Replace values by last "valid" value
    Y_cleaned[start_index:end_index] = Y_cleaned[start_index]
  }
  if (method == 'Arima Fit'){
    # Replace values by arima forecast
    forecast_length = end_index - start_index + 1
    arima_fit = Arima(subset(Y, start=1, end=end_index), order=c(1,1,1), seasonal=c(1,1,1), method="CSS")
    forecasts = forecast(arima_fit, forecast_length)
    
    autoplot(forecasts)
    Y_cleaned[(start_index +1) : (start_index + forecast_length)] = forecasts$mean
  }
  # TRY holt winters
  
  return (Y_cleaned)
}

Y_cleaned = clean_time_series(Y, method='Arima Fit')
print(Y_cleaned)
plot(Y_cleaned)


for (i in 221:235){
  car_weekly$hits[i]=Y_cleaned[i]
}


####################################################
######                                        ######
######               Set parameters           ######
######                                        ######
####################################################


# the dataset to be considered is Y
T <- dim(car_weekly)[1] # sample size
F <- 6  # number of out-of-sample forecasts
T0 <- floor(T*0.8) # Number of in-sample observations (training)
T1 <- T-T0 # Number of out-of-sample observations (validating)

####################################################
######                                        ######
######         EXAMPLE ARIMA(1,0,0)           ######
######                                        ######
####################################################


## Example, fit an AR(1) and assess its forecasting accuracy in-sample

Forecast_SARIMA <- function(Y, T, T0, F, OOS = 1, c_ARIMA, c_SARIMA  = list(order=c(1,1,0)), windo = 'recursive', roll = W, confidence = FALSE)
{
  if(OOS == 0) #not performed by default, unless specified as OOS=0 while calling the function
  {
    Yhat_SARIMA <- matrix(NA,nrow=T,ncol=F)
    CI_SARIMA <- matrix(NA,nrow=T,ncol=F)
    loss_SARIMA <- matrix(NA,nrow=T,ncol=1)
    for(t in (T0+1-F):T) #we start before T0 to make sure we have (multistep) forecasts starting at T0+1
    {
      t0 <- switch(windo, 
                   'recursive' = 1, #by default
                   'rolling' = t-W,
                   'expanding' = max(1,floor(t-T^W)))
      
      fit_SARIMA <- Arima(Y[t0:t], order=c_ARIMA, seasonal = c_SARIMA) #by default: c_SARIMA = c(0,0,0)
      loss_SARIMA[t] <- Loss(as.matrix(fit_SARIMA$residuals))
      for(j in 1:F)
      {
        if(t+j<T+1) 
        {
          Fhat <- forecast(fit_SARIMA, h = j, level =c(95))
          Yhat_SARIMA[t+j,j] = Fhat$mean[j] # Yhat[t+j,j] is the j-step ahead forecast of t+j
          if(confidence == TRUE)
          {
            CI_SARIMA[t+j,j] = Fhat$lower[j] # confidence interval if required
            CI_SARIMA[t+j,F+j] = Fhat$upper[j] # confidence interval if required (can be used to check in-sample the quality of the confidence interval)
          }
        } 
      }
    }
    return(cbind(Yhat_SARIMA,loss_SARIMA))
  }
  if(OOS == 1) #by default
  {# now perform an out-of-sample exercise
    Yhat_SARIMA = matrix(NA,F,1)
    t0 <- switch(windo,
                 'recursive' = 1,
                 'rolling' = T-W,
                 'expanding' = max(1,floor(T-T^W)))
    fit_SARIMA <- Arima(Y[t0:T], order=c_ARIMA, seasonal = c_SARIMA) 
    Yhat_SARIMA <- ts(forecast(fit_SARIMA, h = F)$mean,frequency = 52,start = c(lubridate::year(car_weekly$date)[T],lubridate::week(car_weekly$date)[T]+1))
    lower_SARIMA <- ts(forecast(fit_SARIMA, h = F, level = c(95))$lower,frequency = 52,start = c(lubridate::year(car_weekly$date)[T],lubridate::week(car_weekly$date)[T]+1))
    upper_SARIMA <- ts(forecast(fit_SARIMA, h = F, level = c(95))$upper,frequency = 52,start = c(lubridate::year(car_weekly$date)[T],lubridate::week(car_weekly$date)[T]+1))
    return(cbind(Yhat_SARIMA, lower_SARIMA, upper_SARIMA))
  }
}


auto.arima(Y_cleaned, trace=TRUE)


FIT_SARIMA <- Forecast_SARIMA(Y_cleaned, T, T0, F, OOS=0, c_ARIMA = c(2,1,1), c_SARIMA = list(order=c(1,1,1)))
# you may want to try the different options of the "windo" parameter 

RMSFE_SARIMA <- sqrt(Loss(Y_cleaned[(T0+1):T] - FIT_SARIMA[(T0+1):T,1:F])/(1:F)) # RMSFE scaled by sqrt(horizon)
#scaling the RMSFE is not necessary -- it's just as a benchmar as in the case of the Random Walk, the MSFE increases linearly with the horizon
names(RMSFE_SARIMA) <- 1:F
# this is a row vector containing the MSFE for horizons 1 to F
# this is what you should compare across models
inLoss.ts <- ts(sqrt(FIT_SARIMA[,F+1]),frequency = 52, start = c(lubridate::year(car_weekly$date)[1],lubridate::week(car_weekly$date)[1]))

#plotting in sample one-step ahead (horizon = 1) loss function
# this graph is useful in checking the stability of the forecasting performance
autoplot(inLoss.ts) + geom_line(aes(group=1), color = "red")+ylab("In Sample RMSFE Loss")+xlab("end of estimation sample")

# and compare it with the scaled multistep ahead
RMSFE_SARIMA # the out-of-sample scaled RMSFE at each horizon




####################################################
######                                        ######
######        PROPHET FACEBOOK                ######
######                                        ######
####################################################


df=data.frame(ds= car_weekly$date, y=car_weekly$hits)


####################################################
#####                 In sample             ########
####################################################

IS_prophet <- matrix(NA,nrow=T,ncol=F)

loss_prophet <- matrix(NA,nrow=T,ncol=1)

print(T0+1-F)
t0<-1


for(t in (T0+1-F):(T)) #we start before T0 to make sure we have (multistep) forecasts starting at T0+1
{
 
  fit_prophet<-prophet(weekly.seasonality=FALSE)
  fit_prophet<-add_country_holidays(fit_prophet,country_name='US')
  data <- df[t0:t,]
  
  fit_prophet <- fit.prophet(fit_prophet,data) 
  future <- make_future_dataframe(fit_prophet, periods=F, freq= "week")
  forecast_prophet <- predict(fit_prophet, future)
  
  if ((t+F)<=T){
    loss_prophet[t] <- RMSE(df$y[(t+1):(t+F)],forecast_prophet$yhat)
  }
  for(j in 1:F)
  {
    if(t+1<T+1) 
    { 
      IS_prophet[t+1,j]<- forecast_prophet$yhat[j] # Yhat[t+j,j] is the j-step ahead forecast of t+j
    
  }
}
}

print(tail(IS_prophet))
print(loss_prophet)



###################################################
####### out sample                         ########
###################################################

m<-prophet(weekly.seasonality=FALSE)
m<-add_country_holidays(m,country_name='US')
m<-fit.prophet(m,df)
print(m)
print(Loss(as.matrix(m$residuals)))
print(m)

future <- make_future_dataframe(m, periods=F, freq= "week")
tail(future)

forecast_df <- predict(m, future)
print(tail(forecast_df))
#result_df <- forecast_df[c('yhat')]
result_df <- forecast_df[c('yhat','yhat_lower','yhat_upper')]
print(tail(result_df))

prediction_week <- result_df[261:(261+(F-1)),]

data_df <- result_df$yhat[53:260]
print(data_df)


d_df <- df[53:260, -1]
print(d_df)


rmse_normal <-RMSE(data_df,d_df)
print(rmse_normal)


dyplot.prophet(m, forecast_df)

##################################################
#####       submission week 7                #####
##################################################

IS_week <- 1:T #(Guillaume) check seq(lubridate::week(ymd(car_weekly$date[1])),lubridate::week(ymd(car_weekly$date[T])))
OOS_week <- seq(lubridate::week(ymd(car_weekly$date[T]))+1,lubridate::week(ymd(car_weekly$date[T])+1)+F)

in_sample <- result_df[1:T,]
print(prediction_week)

write.csv(cbind(OOS_week,prediction_week),"HW_week7_Group_12_OOS.csv")
write.csv(cbind(IS_week,in_sample),"HW_week7_Group_12_IS.csv")



####################################################
######                                        ######
######            FORECAST COMNINAISON        ######
######                                        ######
####################################################


#####################################
#####  out sample       #############
#####################################

OOS <- Forecast_SARIMA(Y_cleaned, T, T0, F, OOS=1, c_ARIMA = c(2,1,1), c_SARIMA = list(order=c(1,1,1),period =52))
prediction_combi = data.frame(yhat=(OOS[,1]+ prediction_week$yhat)/2, lower=(OOS[,2]+ prediction_week$yhat_lower)/2 , upper=(OOS[,3]+ prediction_week$yhat_upper)/2)
rownames(prediction_combi) <- rownames(prediction_week)
#print(prediction_combi)

#####################################
#####  in sample         ############
#####################################

IS_SARIMA <- Forecast_SARIMA(Y_cleaned, T, T0, F, OOS=0, c_ARIMA = c(1,1,1), c_SARIMA = list(order=c(1,1,0)))
#print(IS_SARIMA)

insample_combi <- matrix(NA,nrow=T,ncol=F)
for (i in 1:T){
  for (j in 1:F){
    insample_combi[i,j]<- (IS_SARIMA[i,j]+IS_prophet[i,j])/2
  }
}
print(tail(insample_combi))


loss_combi <- matrix(NA,nrow=T,ncol=1)

for(t in (T0+1-F):T){
  if ((t+F)<=T){
    loss_combi[t] <- RMSE(df$y[(t+1):(t+F)],insample_combi[t+1,1:6])
  }
}

#print((loss_combi))
#print(IS_SARIMA[,7])

IS_COMBINATION <-cbind(insample_combi,loss_combi)
print(tail(IS_COMBINATION))


##################################################
#####       submission week 9                #####
##################################################


IS_week <- 1:T #(Guillaume) check seq(lubridate::week(ymd(car_weekly$date[1])),lubridate::week(ymd(car_weekly$date[T])))
OOS_week <- seq(lubridate::week(ymd(car_weekly$date[T]))+1,lubridate::week(ymd(car_weekly$date[T])+1)+F)


print(prediction_week)

write.csv(cbind(OOS_week,OOS),"HW_week9_Group_12_OOS.csv")
write.csv(cbind(IS_week, IS_COMBINATION),"HW_week9_Group_12_IS.csv")


plot.new()
par(mar=c(4,4,3,5))
plot( ts(IS_COMBINATION[203:254,7]),col='red',axes=F,xlab="",ylab="") 
axis(2, ylim=c(0,10),col="blue",col.axis="blue",at=seq(0, 10, by=2)) 
mtext("Axe de la courbe bleue",side=2,line=2.5,col="blue")  

plot(ts(IS_SARIMA[203:254,7]), col='blue',axes=F,xlab="",ylab="") 
plot( ts(IS_COMBINATION[203:254,7]),col='red',axes=F,xlab="",ylab="") 
IS_prophet
plot( ts(loss_prophet[203:254,1]),col='red',axes=F,xlab="",ylab="") 

####################################################
######                                        ######
######             REQUIRED OUTPUT            ######
######                                        ######
####################################################



## ASSUMING THE AR(1) is your preferred choice: produce and save in-sample (IS) and out-of-sample (OOS) forecasts (T0+1-F):T  
IS_SARIMA <- Forecast_SARIMA(Y_cleaned, T, T0, F, OOS=0, c_ARIMA = c(1,1,1), c_SARIMA = list(order=c(1,1,0)))

IS_week <- 1:T #(Guillaume) check seq(lubridate::week(ymd(car_weekly$date[1])),lubridate::week(ymd(car_weekly$date[T])))
write.csv(cbind(IS_week,IS_SARIMA),"HW_week9_Group_12_IS.csv")
OOS <- Forecast_SARIMA(Y_cleaned, T, T0, F, OOS=1, c_ARIMA = c(1,1,1), c_SARIMA = list(order=c(1,1,0),period =52))
OOS_week <- seq(lubridate::week(ymd(car_weekly$date[T]))+1,lubridate::week(ymd(car_weekly$date[T])+1)+F)
write.csv(cbind(OOS_week,OOS),"HW_week6_Group_12_OOS.csv")
