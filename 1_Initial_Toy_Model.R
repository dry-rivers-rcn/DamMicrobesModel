#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Title: Initial Toy Model
#Coder: Nate Jones (cnjones7@ua.edu)
#Date: 10/15/2019 
#Purpose: Begin exploring Toy Model to explore drying and its impact on Biogeochemistry
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Next Steps 
# -Add conservative solute flux
# -Add reactive solut flux

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#1.0 Setup Workspace------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Clear memory (delete this once incoporated into larger workflow)
rm(list=ls())

#download relevant packages
library(tidyverse)

#User defined variables
reach_length <- 10000       # Each reach is 10 km in length
reach_slope <- 0.0003       # Slope of each reach 
reach_bkf_depth <- 1    # Bankful depth
reach_wd_ratio <- 5       # Width to depth ratio of bankful channel
reach_roughness<-0.07      # Mannings N
simulation_length<-6       # days


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2.0 Estimate hydrologic parameters---------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2.1 Hydraulic Geometry---------------------------------------------------------
#Estimate bankful width (assume rectangular channel)
reach_bkf_width<-reach_wd_ratio*reach_bkf_depth

#Estimate bankfull flow (cms) using mannings equation
Q_bkf<- 
  (1/reach_roughness)*                             #k/n
  (reach_bkf_depth*reach_bkf_width)*               #Bankful Area
  (reach_bkf_width^(2/3))*     #Hydraulic Radius ^ (2/3)
  (reach_slope^.5)                                 #slope^0.5

#Estimate Bankful Velocity 
v_bkf<- Q_bkf/(reach_bkf_depth*reach_bkf_width)

#2.1 Inflow hydrograph ---------------------------------------------------------
#Define time step so that the Courant number (v*dt/dx) is less than or equal to 1
timestep<- floor((reach_length/v_bkf))

#Create tibble
df<-tibble(
      #Create 2 time steps
      time = seq(0,simulation_length*86400, timestep),
      #Create flow 
      Q_in = Q_bkf*sin(pi*time/86400) + Q_bkf/10
    ) %>% 
  mutate(Q_in = if_else(Q_in<=Q_bkf/10, Q_bkf/10, Q_in))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#4.0 Create function to estimate fluxes-----------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#4.1 Kinematic Wave-------------------------------------------------------------
kinematic_wave<-function(Q_in,                    #Inflow
                         Q_in_1,                  #Qout from previous time step
                         dt = timestep,           #time step (seconds)
                         dx = reach_length,       #Reach length
                         slope =  reach_slope,    #Channel Slope
                         width = reach_bkf_width, #Channel Width
                         n = reach_roughness      #mannings roughness
                        ){
  #Estiate model parameters
  beta <- 3/5  #Common Assumption (See Chow 1988)
  alpha<- (n*(width^(2/3))/(slope^0.5))^beta
  
  #Estimate Reach Outflow
  Q_out<- Q_in + (((Q_in^(1-beta))/(alpha*beta*dx))*(Q_in_1-Q_in))*dt
  
  #Export Qout
  #Q_out<-if_else(Q_out<0, 0, Q_out)
  Q_out
}

#4.2 Conservative Solute Flux---------------------------------------------------


#4.3 Reactive solute flux-------------------------------------------------------



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#5.0 Create function to estimate conservative solute flux-----------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
df<-df %>%
  #Run model
  mutate(Q_out = 
     pmap_dbl(
       #Model Inputs
       list(Q_in =   Q_in, 
            Q_in_1 = lag(Q_in, default = 0)), 
       #Function
       kinematic_wave)) %>% 
  #Add time_hr 
  mutate(time_hr = time/3600)

plot(df$time_hr, df$Q_in, type="l", xlab="Time [hr]", ylab="Flow [cms]", ps=12, cex.lab=14/12, cex.axis=10/12)
points(df$time_hr, df$Q_out, col="red", type="l")

