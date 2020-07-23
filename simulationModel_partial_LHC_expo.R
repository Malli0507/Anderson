library(deSolve)
library(tidyverse)
library(cowplot)

# set working directory
path <- "C:\\Users\\Carolin\\Documents\\Rscripts\\Covid-19\\"
setwd(path)

### the ODE model
covidmodel <- function(t, vars, parms) {
  with(as.list(c(vars, parms)), {
    foi <- 
      (beta1 * I1 +
         beta1 * I2m + beta1 * I3m +
         beta2 * b_par * I2s + b_par * beta3 * I3s +
         beta2 * m_par * Y2 + beta3 * m_par * Y3) / (popsize)
    return(
      list(
        c(
          dS = -foi * S,
          dE = foi * S - sigma * E,
          dI1 = sigma * E - gamma1 * I1,
          dI2m = gamma1 * I1 * (1 - p_par) - gamma2 * I2m,
          dI3m = gamma2 * I2m - gamma3 * I3m,
          dI2s = gamma1 * I1 * p_par - (alpha + gamma2) * I2s,
          dI3s = gamma2 * I2s - gamma3 * I3s,
          dY2 = alpha * I2s - gamma2 * Y2,
          dY3 = gamma2 * Y2 - gamma3 * Y3,
          dInf = foi * S,
          dRep = alpha * I2s
        )
      )
    )
  })
}

### the function doing the simulations
# There are default parameter values for everything, 
# so simply calling the function produces a matrix with simulated variables.
# Alternative parameters can be chosen as argument to the function, 
# e.g. calling "covidsim(alpha = 0.3)"
covidsim <- function(
  parameters = c(
    beta1 = 0.25, beta2 = 0.16, beta3 = 0.016,
    sigma = 1, gamma1 = 0.2, gamma2 = 0.14, gamma3 = 0.14,
    alpha = 0.5, p_par = 0.5, b_par = 0.8, m_par = 0.01,
    start_dist = 60, end_dist = 108, effect_dist = 1
 ),
  popsize = 60e+06, seedsize = 10, burnintime = 30,
  timewindow = 365, lockdown_factor, after_lockdown_factor,
  ...
) {
  altparms <- list(...)
  for(i in names(altparms)) {
   	parameters[i] <- altparms[[i]]
  }
  
  initialstate <- c(
    S = popsize, E = seedsize, I1 = 0,
    I2m = 0, I3m = 0,
    I2s = 0, I3s = 0,
    Y2 = 0, Y3 = 0,
    Infected = 0, Reported = 0
  )
  
  startparms <- c(popsize = popsize, parameters)
  startparms["m_par"] <- 1
  startparms["b_par"] <- 1
  simres <- ode(
    y = initialstate,
    times = seq(0, burnintime, 1),
    func = covidmodel,
    parms = startparms
  )
  
  nodistparms <- c(popsize = popsize, parameters)
  nodistparms["m_par"] <- 1
  nodistparms["b_par"] <- 1
  simres <- ode(
    y = tail(simres, 1)[, -1],
    times = seq(0, min(nodistparms[["start_dist"]], timewindow), 1),
    func = covidmodel,
    parms = nodistparms
  )
  
  if(timewindow <= nodistparms[["start_dist"]]) return(simres)
  distparms <- c(popsize = popsize, parameters)
  time <- distparms[["start_dist"]]
  beta1 <- distparms[["beta1"]]
  beta2 <- distparms[["beta2"]]
  beta3 <- distparms[["beta3"]]

  while(time < distparms[["end_dist"]])
  {
  	distparms[["beta1"]] <- beta1 * lockdown_factor[time - 83 +1]
  	distparms[["beta2"]] <- beta2 * lockdown_factor[time - 83 +1]
  	distparms[["beta3"]] <- beta3 * lockdown_factor[time - 83 +1]

  simres <- rbind(
    head(simres, -1),
    ode(
      y = tail(simres, 1)[, -1],
      times = seq(time, time+1, 1),
      func = covidmodel,
      parms = distparms
    )
  )
  time <- time + 1
  }  

  if(timewindow <= nodistparms[["end_dist"]]) return(simres)
  nodistparms <- c(popsize = popsize, parameters)
  #nodistparms[["beta1"]] <- nodistparms[["beta1"]] * 0.6
  #nodistparms[["beta2"]] <- nodistparms[["beta2"]] * 0.6
  #nodistparms[["beta3"]] <- nodistparms[["beta3"]] * 0.6

  time <- nodistparms[["end_dist"]]
  beta1 <- nodistparms[["beta1"]]
  beta2 <- nodistparms[["beta2"]]
  beta3 <- nodistparms[["beta3"]]

  while(time < timewindow)
  {
	nodistparms[["beta1"]] <- beta1 * after_lockdown_factor[time - 186 +1]
  	nodistparms[["beta2"]] <- beta2 * after_lockdown_factor[time - 186 +1]
  	nodistparms[["beta3"]] <- beta3 * after_lockdown_factor[time - 186 +1]

  	simres <- rbind(
    	head(simres, -1),
    	ode(
      	y = tail(simres, 1)[, -1],
      	times = seq(time, time+1, 1),
      	func = covidmodel,
      	parms = nodistparms
    		)
  	)
  time <- time + 1

  }
  
  return(simres)
}


# READ IN PARAMS FROM FILE ##########################################################

library(XLConnect)

# first drawn an incubation period and then derive
# gamma_1 and sigma under the assumption that
# the exposed time and the time in the infected compartment I1
# are each half the incubation period
# needs to be called before calculating beta
get_gamma1_sigma <- function(params, n)
{
	IP <- runif(n, params$incuPer[2], params$incuPer[3])
	gamma_1 <- 1/(0.5*IP)
	sigma <- 1/(0.5*IP)

	return(list(sigma=sigma, gamma_1=gamma_1))
}

# calculate parameters not defined in file and
# draw n samples from parameter distribution of
# i) p - proportion symptomatic
# ii) sigma + gamma1 - time to infectiousness
# iii) q2 and q3 - relative infectiousness of asymptomatic individuals
getParams <- function(params, n)
{
	R0 <- rep(params$R_0[1], n)
	gamma2 <- rep(params$gamma_2[1], n)
	gamma3 <- rep(params$gamma_3[1], n)

	p <- runif(n, params$p[2], params$p[3])

	alpha <- rep(params$alpha[1], n)
	b <- rep(params$b[1], n)
	m <- rep(params$m[1], n)

	q1 <- rep(params$q_1[1], n)
	q2 <- runif(n, params$q_2[2], params$q_2[3])
	q3 <- runif(n, params$q_3[2], params$q_3[3])

	temp <- get_gamma1_sigma(params, n)
	sigma <- temp[[1]]
	gamma1 <- temp[[2]]

	beta <- R0 / (1/gamma1 + (1 - p) * (1/gamma2 + 1/gamma3) + p/gamma2 * (q2 + gamma2*q3/gamma3))
	beta1 <- beta * q1
	beta2 <- beta * q2
	beta3 <- beta * q3

	parameters <- data.frame(R0=R0, beta1=beta1, beta2=beta2, beta3=beta3, 
					sigma=sigma, gamma1=gamma1, gamma2=gamma2, gamma3=gamma3,
						p_par=p, alpha=alpha, b_par=b, m_par=m)

	return(parameters)
}

getParamsCentral <- function(params, n)
{
	R0 <- params$R_0[1]
	gamma2 <- params$gamma_2[1]
	gamma3 <- params$gamma_3[1]
	p <- params$p[1]
	alpha <- params$alpha[1]
	b <- params$b[1]
	m <- params$m[1]

	q1 <- params$q_1[1]
	q2 <- params$q_2[1]
	q3 <- params$q_3[1]

	IP <- params$incuPer[1]
	gamma1 <- 1/(0.5*IP)
	sigma <- 1/(0.5*IP)

	beta <- R0 / (1/gamma1 + (1 - p) * (1/gamma2 + 1/gamma3) + p/gamma2 * (q2 + gamma2*q3/gamma3))
	beta1 <- beta * q1
	beta2 <- beta * q2
	beta3 <- beta * q3

	parameters <- data.frame(R0=R0, beta1=beta1, beta2=beta2, beta3=beta3, 
					sigma=sigma, gamma1=gamma1, gamma2=gamma2, gamma3=gamma3,
						p_par=p, alpha=alpha, b_par=b, m_par=m)

	return(parameters)
}




# WITH LOCKDOWN AND ISOLATION ########################################################################


set.seed(123)

filename <- "params3.xlsx"

wb <- loadWorkbook(paste0(path, filename))
params.lit <- readWorksheet(wb, sheet="Sheet1")
paramnames <- params.lit$Param
params <- as.data.frame(t(params.lit[,2:4]))
names(params) <- paramnames

n <- 100
parameters <- getParams(params, n)
parameters$start_dist <- rep(83, n)
parameters$end_dist <- rep(186, n)
#parameters$effect_dist <- rep(0.01, n)

#parameters$m_par <- rep(1, n)
#parameters$b_par <- rep(1, n)


get_lockdown_R <- function(params, t)
{
	a <- params$a
	b <- params$b
	C <- params$C

	R <- a * exp(-b*t) + C
	return(R)
}

lockdown_R <- get_lockdown_R(params=list(a=1.8, b=0.06491, C=0.7), t=seq(0,103))
lockdown_factor <- lockdown_R / parameters$R0[[1]]


get_after_lockdown_R <- function(params, t)
{
	a <- params$a
	S <- params$S
	k <- params$k

	R <- a*S / (a + (S-a)*exp(-S*k*t))
	return(R)
}



after_lockdown_R <- get_after_lockdown_R(a=0.71, S=2.5, k=0.0077258), t=seq(0,180))
after_lockdown_factor <- after_lockdown_R / parameters$R0[[1]]


# RUN SIMULATIONS ########################################################################


sims <- list()
for(i in 1:nrow(parameters))
{
	sims[[i]] <- as.data.frame(covidsim(parameters[i,], popsize = 60e+06, seedsize = 100, burnintime = 30, timewindow = 365, 
			lockdown_factor=lockdown_factor, after_lockdown_factor=after_lockdown_factor))
	sims[[i]]$rep <- rep(i, nrow(sims[[i]]))
}

df <- do.call(rbind, sims)
df$rep <- as.factor(df$rep)

df$Incidence <- df$Infected - lag(df$Infected)
df$Inc.reported <- df$Reported - lag(df$Reported)
df$Incidence[is.na(df$Incidence)] <- 0
df$Incidence[df$Incidence<0] <- 0
df$Inc.reported[is.na(df$Inc.reported)] <- 0
df$Inc.reported[df$Inc.reported<0] <- 0



#############################################################################################
# PLOT INDIVIDUAL-LEVEL SIMULATIONS
#############################################################################################

plot.path <- "C:\\Users\\Carolin\\Documents\\Rscripts\\Covid-19\\SmoothFigures\\"

p1 <- ggplot(df, aes(x=time, y=E, colour=rep)) + theme_classic() + geom_line() + 
	theme(legend.position="none")  
print(p1)

p2 <- ggplot(df, aes(x=time, y=Infected, colour=rep)) + theme_classic() + geom_line() +
	theme(legend.position="none")
print(p2)


q1 <- ggplot(df, aes(x=time, y=Inc.reported, colour=rep)) + theme_classic() + geom_line() + 
	theme(legend.position="none")  
print(q1)

q2 <- ggplot(df, aes(x=time, y=Reported, colour=rep)) + theme_classic() + geom_line() +
	theme(legend.position="none")
print(q2)

# output results from all iterations
write.table(file=paste0(plot.path, "lhc.lockdown.csv"), df, row.names=FALSE, sep=",")



#############################################################################################
# PLOT ALL SIMULATIONS MEAN AND 95% CRI OF INCIDENCE
#############################################################################################


inc.mean <- aggregate(df[, c("E")], by=list(df$time), mean)$x
inc.95CrI <- aggregate(df[, c("E")], by=list(df$time), quantile, c(0.025, 0.975))[,2]
df.plot <- as.data.frame(cbind(inc.mean, inc.95CrI))
names(df.plot) <- c("mean", "pc_2.5", "pc_97.5")
df.plot$time <- df$time[df$rep==1]


df.plot.cut <- df.plot[1:274, ]


q <- ggplot(df.plot.cut) + theme_classic() + theme(axis.title=element_text(size=16), axis.text=element_text(size=16)) +
     geom_ribbon(aes(x=time, ymin=pc_2.5, ymax=pc_97.5), fill="cadetblue", alpha=0.5) +
     geom_line(aes(x=time, y=mean), colour="darkslategrey") +
     xlab("Time in days") + ylab("Daily incidence") +
     coord_cartesian(xlim=c(0, 285), ylim=c(0, (max(df.plot.cut$pc_97.5)+0.05*df.plot.cut$pc_97.5)), expand=FALSE) +
     scale_x_continuous(breaks=c(1, 31, 60, 91, 121, 152, 182, 213, 244, 274), 
				labels=c("January", "February", "March", "April", "May", "June", "July", "August",
					"September", "October")) +
     theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=0.5), axis.title.x=element_blank(),
	axis.text.y=element_blank(), axis.ticks.y=element_blank())
print(q)
ggsave(paste0(plot.path, "mean.95CI.lock.down.logistic.growth.pdf"), plot=q, dpi=300)	







