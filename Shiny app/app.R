# working directory
# setwd("~/Documents/R/Shiny")

# set seed for replication
set.seed(1234)

# load required packages
library(kableExtra)
library(knitr)
library(magrittr)

######################################
# Function
######################################

OLE <- function(N, hE, hC, t1, t2, t_star) { 
    
    ### Function input ---
    
    # N = patient horizon 
    # hE = hazard rate experimental treatment 
    # hC = hazard rate control treatment 
    # t1 = starting point of trial (in years)
    # t2 = endpoint of trial (in years)
    # t_star = time until N new onset  
    
    ### Computations ---
    
    # 1. Calculate the sample size of the different groups
    n1 <- 0.5*N*(t1/t_star)           # Patients in the trial receiving the standard treatment  
    n2 <- 0.5*N*(t1/t_star)           # Patients in the trial receiving the experimental treatment
    n3 <- N*((t2-t1)/t_star)          # Patients not in the trial, with disease onset between t1 and t2
    n4 <- N*((t_star-t2)/t_star)      # Patients not in the trial, with disease onset after t2. 
    
    # 2. Calculating the average life expectancy when E is not superior 
    C.LE1 <- 1/hC                                                                   # group only receives C 
    C.LE2 <- (1/hE^2) * (hE + (((exp(t1*hE)-1)*exp(-t2*hE)*(hE-hC))/(t1*hC)))       # group receives E first, after t' C
    C.LE3 <- 1/hC                                                                   # group only receives C
    C.LE4 <- 1/hC                                                                   # group only receives C
    
    # 3. Calculating the average life expectancy when E is superior
    E.LE1 <- (1/hC^2) * (hC + (((exp(t1*hC)-1)*exp(-t2*hC)*(hC-hE))/(t1*hE)))       # group receives C first, after t' E                                                       
    E.LE2 <- 1/hE                                                                   # group only receives E
    E.LE3 <- (1/hC^2) * (hC + (((exp(t1*hC)-1)*exp(-t2*hC)*(hC-hE))/(t1*hE)))       # group receives C first, after t' E                                                      
    E.LE4 <- 1/hE                                                                   # group only receives E
    
    # 4. Calculating the overall life expectancy for the entire population N for situation E is superior and E is not          superior
    
    #### By multiplying group size with LE for individuals 
    C.OLE <- (n1 * C.LE1) + (n2 * C.LE2) + (n3 * C.LE3) + (n4 * C.LE4) 
    E.OLE <- (n1 * E.LE1) + (n2 * E.LE2) + (n3 * E.LE3) + (n4 * E.LE4) 
    
    # 5. Calculating the power 
    
    ### First calculate the expected number of events for each group 
    D.E <- n2 * (1 + ((exp(-1*(t2-t1)*hE) - exp(-1*t2*hE))/((t2-t1)*hE - t2*hE)))
    D.C <- n1 * (1 + ((exp(-1*(t2-t1)*hC) - exp(-1*t2*hC))/((t2-t1)*hC - t2*hC)))
    
    ### Calculate the total expected number of events by adding  both groups 
    D <- D.E + D.C
    
    #### PropE = proportion in E group, in the case of 1:1 allocation this is 0.5
    PropE <- 0.5
    
    ### Calculate the mean under the alternative hypothesis 
    ### Log-rank test statistic T has distribution (Schoenfeld, 1981 Biometrika)
    muHa <- log(hE/hC) * sqrt(PropE * (1-PropE) * D)
    
    ### Determine the critical value for one sided alpha of 0.05. 
    Cr <- qnorm(0.025, 0, 1) 
    
    #### Determine the power by taking the area left of the critical value 
    pwr <- pnorm(Cr, muHa, 1)
    
    # 6. Weigh the two possible outcomes by the power 
    U <- (pwr*E.OLE) + ((1-pwr)*C.OLE) 
    
    
    ### Output --- 
    
    # 7. Take the results and put in a list 
    result <- list(
        'Overall life expectancy E not superior' = C.OLE,
        'Overall life expectancy E superior' = E.OLE,
        'Utility' = U
    )
    
    # 6. Return a list with the results 
    return(result)
    
    # End of function 
}


######################################
# Generating data 
# E treatment better 
######################################
# Data input
hE <-0.2
hC <- 0.4
t_star <- 10
N <- 3000

t1.seq <- seq(0, t_star, 0.1)
t2.seq <- seq(0, t_star, 0.1)

# Create matrix with possible combinations of t1 and t2
m <- expand.grid(hE=hE, hC=hC, t_star=t_star, N=N, t1=t1.seq, t2=t2.seq)    # t2 = column name, t2seq = content
m <- m[m$t2 >= m$t1 & m$t1>0,]                                              # add conditions: end of follow up (t2) needs                                                              to be after start of trial (t1), t1 must be larger than 0. 

# Add extra columns to matrix for C.OLE, E.OLE and Utility (U)
m$uc <- NA # 'Utility' C
m$ue <- NA # 'Utility' E
m$U <- NA # Overall utility weighted by the power 

# For loop
for(i in 1:nrow(m)) {
    m$uc[i] <- OLE(N = m$N[i], hE = m$hE[i], hC = m$hC[i], t1 = m$t1[i], t2 = m$t2[i], t_star = m$t_star[i])$'Overall life expectancy E not superior'
    m$ue[i] <- OLE(N = m$N[i], hE = m$hE[i], hC = m$hC[i], t1 = m$t1[i], t2 = m$t2[i], t_star = m$t_star[i])$'Overall life expectancy E superior'
    m$U[i] <- OLE(N = m$N[i], hE = m$hE[i], hC = m$hC[i], t1 = m$t1[i], t2 = m$t2[i], t_star = m$t_star[i])$'Utility'
}

# save as data
as.data.frame(m)

######################################
# Generating data 
# C treatment better 
######################################
t1.seq <- seq(0, t_star, 0.1)
t2.seq <- seq(0, t_star, 0.1)

# Data input
hE <-0.4
hC <- 0.2

# Create matrix with possible combinations of t1 and t2
c <- expand.grid(hE=hE, hC=hC, t_star=t_star, N=N, t1=t1.seq, t2=t2.seq)    # t2 = column name, t2seq = content
c <- c[c$t2 >= c$t1 & c$t1>0,]                                              # add conditions: end of follow up (t2) needs                                                              to be after start of trial (t1), t1 must be larger than 0. 

# Add extra columns to matrix for C.OLE, E.OLE and Utility (U)
c$uc <- NA # 'Utility' C
c$ue <- NA # 'Utility' E
c$U <- NA # Overall utility weighted by the power 

# For loop
for(i in 1:nrow(c)) {
    c$uc[i] <- OLE(N = c$N[i], hE = c$hE[i], hC = c$hC[i], t1 = c$t1[i], t2 = c$t2[i], t_star = c$t_star[i])$'Overall life expectancy E not superior'
    c$ue[i] <- OLE(N = c$N[i], hE = c$hE[i], hC = c$hC[i], t1 = c$t1[i], t2 = c$t2[i], t_star = c$t_star[i])$'Overall life expectancy E superior'
    c$U[i] <- OLE(N = c$N[i], hE = c$hE[i], hC = c$hC[i], t1 = c$t1[i], t2 = c$t2[i], t_star = c$t_star[i])$'Utility'
}

# save as data
as.data.frame(c)

######################################
# Generating data 
# Treatments equal
######################################
t1.seq <- seq(0, t_star, 0.1)
t2.seq <- seq(0, t_star, 0.1)

# Data input
hE <-0.2
hC <- 0.2

# Create matrix with possible combinations of t1 and t2
x <- expand.grid(hE=hE, hC=hC, t_star=t_star, N=N, t1=t1.seq, t2=t2.seq)    # t2 = column name, t2seq = content
x <- x[x$t2 >= x$t1 & x$t1>0,]                                              # add conditions: end of follow up (t2) needs                                                              to be after start of trial (t1), t1 must be larger than 0. 

# Add extra columns to matrix for C.OLE, E.OLE and Utility (U)
x$uc <- NA # 'Utility' C
x$ue <- NA # 'Utility' E
x$U <- NA # Overall utility weighted by the power 

# For loop
for(i in 1:nrow(c)) {
    x$uc[i] <- OLE(N = x$N[i], hE = x$hE[i], hC = x$hC[i], t1 = x$t1[i], t2 = x$t2[i], t_star = x$t_star[i])$'Overall life expectancy E not superior'
    x$ue[i] <- OLE(N = x$N[i], hE = x$hE[i], hC = x$hC[i], t1 = x$t1[i], t2 = x$t2[i], t_star = x$t_star[i])$'Overall life expectancy E superior'
    x$U[i] <- OLE(N = x$N[i], hE = x$hE[i], hC = x$hC[i], t1 = x$t1[i], t2 = x$t2[i], t_star = x$t_star[i])$'Utility'
}

# save as data
as.data.frame(x)

######################################
# Shiny app
######################################

# Define UI 
ui <- fluidPage(
    title = "Life expectancy",
    sidebarLayout(
        sidebarPanel(
            conditionalPanel(
                'input.dataset === "m"',
                checkboxGroupInput("show_vars", "Columns in m to show:",
                                   names(m), selected = names(m))
            )
            ,
            conditionalPanel(
                'input.dataset === "c"',
                helpText("Click the column header to sort a column.")
            ),
            conditionalPanel(
                'input.dataset === "x"',
                helpText("Display 5 records by default.")
            )
        ),
        mainPanel(
            tabsetPanel(
                id = 'Outcome of clinical trial',
                tabPanel("Experimental treatment superior", DT::dataTableOutput("mytable1")),
                tabPanel("Control treatment treatment superior", DT::dataTableOutput("mytable2")),
                tabPanel("No difference between two treatments", DT::dataTableOutput("mytable3"))
            )
        )
    )
)


# Define server 
server <- function(input, output) {
    
    # choose columns to display
    diamonds2 = m[sample(nrow(m), 1000), ]
    output$mytable1 <- DT::renderDataTable({
        DT::datatable(m[, input$show_vars, drop = FALSE]) 
    })
    
    # sorted columns are colored now because CSS are attached to them
    output$mytable2 <- DT::renderDataTable({
        DT::datatable(c[, input$show_vars, drop = FALSE])
    })
    
    # customize the length drop-down menu; display 5 rows per page by default
    output$mytable3 <- DT::renderDataTable({
        DT::datatable(x[, input$show_vars, drop = FALSE])
    })
    
}

# Run the application 
shinyApp(ui = ui, server = server)


