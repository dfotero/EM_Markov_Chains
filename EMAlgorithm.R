
# EM-Algorithm ------------------------------------------------------------
# Author: Daniel Felipe Otero-Leon
# email: dfotero@umich.edu
# The algorithm is based on Sherlaw-johnson  AC,  Gallivan  S,  Burridge  J (1995)  Estimating  a  Markov  Transition  Matrix  from  Observational  Data.  Palgrave  Macmillan  Journals  on  behalf  of  the  Operational  Research  Society  Stable



library(dplyr)
library(sqldf)
library(expm)
library(lubridate)
library(progress)


#'@title Observation Data
#'@description Organize the observational data.
#'@param data Observational data: has to have the following structure: Date of occurance (ddate), system state (state)
#'@return A data frame called O that contains previous state (u), current state (v), the time spent between the ocurrance of these two states (Arrival), and counts the number of observations for each combination (n).
observationData<-function(data)
{
  data<-data.frame(data,u=0,v=0,Arrival=0)
  data<-data[order(data$ddate),]
  rows<-nrow(data)
  
  pb <- progress_bar$new(format = " Observational data running [:bar] :percent eta: :eta",total = rows, clear = FALSE, width= 60)
  
  for(i in 2:rows)
  {
    pb$tick()
    data$Arrival[i]<-as.numeric(data$ddate[i]-data$ddate[i-1])
    data$v[i]<-data$state[i]
    data$u[i]<-data$state[i-1]
    Sys.sleep(1 / 100)
  }
  
  data<-data[data$Arrival>0,]
  
  O<-count(data,u,v,Arrival)
  
  return(O)
}

#'@title Initial probability matrix
#'@description From the data, creates an initial Markov Chain which represents the frequency of each combination regardless of time between states.
#'@param obs Organized observational data: The organized observational data frame that contains  previous state (u), current state (v), the time spent between the ocurrance of these two states (Arrival), and counts the number of observations for each combination (n).
#'@param states Vector with the system states.
#'@return A square stochastic matrix of size length(states).
probInitial<-function(obs,states)
{
  numL<-nrow(states)
  prob<-matrix(0,nrow=numL,ncol=numL)
  
  for(i in 1:numL)
  {
    theSum<-0
    
    for(j in 1:numL)
    {
      
      if(nrow(obs[obs$u==i & obs$v==j,])>0)
      {
        prob[i,j]<-mean(obs[obs$u==i & obs$v==j,]$n)
        theSum<-theSum+mean(obs[obs$u==i & obs$v==j,]$n)
      }
      else
      {
        prob[i,j]<-1
        theSum<-theSum+1
      }
    }
    
    prob[i,]<-prob[i,]/theSum
    
  }
  return(prob)
}


#'@title Transition Matrix
#'@description A iteration of the EM-algorithm estimating a new transition matrix from the observations and the previors transition matrix.
#'@param obs Organized observational data: The organized observational data frame that contains  previous state (u), current state (v), the time spent between the ocurrance of these two states (Arrival), and counts the number of observations for each combination (n).
#'@param states Vector with the system states.
#'@param prob A square stochastic matrix of size length(states).
#'@return A square stochastic matrix of size length(states).
transitionMatrix<-function(obs,states,prob)
{
  numL<-nrow(states)
  S<-matrix(0,nrow=numL, ncol=numL)
  
  for(i in 1:numL)
  {
    for(j in 1:numL)
    {
      
      for(u in 1:numL)
      {
        for(v in 1:numL)
        {
          theObs<-obs[obs$u==u & obs$v==v,]
          numW<-nrow(theObs)
          
          if(numW>0)
          {
            for(w in 1:numW)
            {
              time<-theObs$Arrival[w]
              o<-theObs$n[w]
              for(l in 0:(time-1))
              {
                if((Prob %^% time)[u,v]>0)
                {
                  S[i,j]<-S[i,j]+o*(prob %^% l)[u,i]*prob[i,j]*(prob %^% (time-l-1))[j,v]/(prob %^% time)[u,v]
                }
              }
            }
          }
        }
      }
      
    }
  }
  return(S)
}

#'@title The EM-algorithm (main function)
#'@description Creates a Mrkov chain from observational data
#'@param data Observational data: has to have the following structure: Date of occurance (ddate), system state (state)
#'@param states Vector with the system states.
#'@param epsilon Error where the algorithm should stop. Suggest 10^(-3).
#'@param days Days between epochs.
#'@return A square stochastic matrix of size length(states).
EMAlgorithm<-function(data,states,epsilon,days)
{
  obs<-observationData(data)
  
  obs$Arrival<-ceiling(obs$Arrival/days)
  obs<-aggregate(obs$n,by=list("u"=obs$u,"v"=obs$v,"Arrival"=obs$Arrival),sum)
  names(obs)[4]<-"n"
  
  prob<-probInitial(obs,states)
  isDif<-TRUE
  numL<-nrow(states)
  ind<-1
  while (isDif)
  {
    ind<-ind+1
    
    isDif<-FALSE
    S<-transitionMatrix(obs,states,prob)
    newP<-matrix(0,nrow=numL,ncol=numL)
    themax<-0
    for(i in 1:numL)
    {
      theSum<-sum(S[i,])
      for(j in 1:numL)
      {
        if(theSum>0)
        {
          newP[i,j]<-S[i,j]/theSum
          if(abs(newP[i,j]-prob[i,j])>themax)
          {
            themax<-abs(newP[i,j]-prob[i,j])
          }
          if(abs(newP[i,j]-prob[i,j])>epsilon)
          {
            isDif<-TRUE
          }
        }
      }
    }
    prob<-newP
  }
  print(ind)
  return(prob)
}
