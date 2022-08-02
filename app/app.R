# App to fit Random Forest models to predict 14-3-3 Binding Sites

library(randomForest)
library(tidyverse)
library(httr)
library(jsonlite)
library(shiny)
library(DT)
library(gtools)

# read in data to predict for and data to train on
pred.dat <- read.csv('Normalization Database including CST.csv',header=TRUE)
train.dat <- read.csv('withID.csv',header=TRUE)
# only want S and T because those are the only ones they can bind to
train.dat <- train.dat[train.dat$Amino.Acid %in% c("S","T"),]
train.dat$Response <- as.factor(train.dat$Response)
train.dat$Amino.Acid <- as.factor(train.dat$Amino.Acid)
train.dat$Site <- train.dat$Position1
# train data should have 0s in where no normalization score??
train.dat$normalization <- as.numeric(train.dat$normalization)
train.dat$normalization[is.na(train.dat$normalization)]<-0
# train various models (for multiple columns of tool)
# they asked for the model to be trained without the protein they are writing
# mainly about, which is TNK1 (Q13470), so taking it out here
train.leaveout <- train.dat#[-which(train.dat$uniprot.ID=="Q15418"),]
full <- randomForest(formula=Response~Amino.Acid+Iupred.score+ANN+PSSM+SVM+
                         normalization,data=train.leaveout,mtry=2,importance=TRUE,
                     na.action=na.roughfix)
red <- randomForest(formula=Response~Amino.Acid+Iupred.score+
                        normalization,data=train.leaveout,mtry=2,importance=TRUE,
                    na.action=na.roughfix)
con <- randomForest(formula=Response~ANN+PSSM+SVM,data=train.leaveout,mtry=2,importance=TRUE,
                    na.action=na.roughfix)
# algorithm mimicking their site
# cutoffs
#ANN - Artificial Neural Network (cut-off = 0.55)
#PSSM - Position-Specific Scoring Matrix (cut-off = 0.80)
#SVM - Support Vector Machine (cut-off = 0.25)
#Consensus - Average of the scores provided by the three methods (cut-off = 0.50)
site.copy <- function(new.dat) {
    annind <- new.dat$ANN > 0.55
    pssmind <- new.dat$PSSM > 0.8
    svmind <- new.dat$SVM > 0.25
    no.high <- annind+pssmind+svmind
    # indices for the sites at each level
    threes <- which(no.high==3)
    twos <- which(no.high==2)
    ones <- which(no.high==1)
    zeros <- which(no.high==0)
    # sort each level by Consensus within it
    s.3 <- new.dat[threes[order(new.dat$Consensus[threes],decreasing=TRUE)],]
    s.2 <- new.dat[twos[order(new.dat$Consensus[twos],decreasing=TRUE)],]
    s.1 <- new.dat[ones[order(new.dat$Consensus[ones],decreasing=TRUE)],]
    s.0 <- new.dat[zeros[order(new.dat$Consensus[zeros],decreasing=TRUE)],]
    rbind(s.3,s.2,s.1,s.0)
}

# Define UI for application
ui <- fluidPage(fluidRow(
  column(
    6,
    fluidRow(
      img(
        src = 'lablogo.jpg',
        align = "left",
        height = 160,
        width = 160
      ),
      img(
        src = 'Final logo 1.png',
        align = "right",
        height = 100,
        width = 160
      ),
      align = "center"
    ),
    strong("Enter a list of proteins for output:"),
    verbatimTextOutput("protexs"),
    textInput(
      inputId = "uniprot.ID.list",
      label = "Protein ID",
      value = "Q96PX8,Q16695,P27986"
    ),
    sliderInput(
      inputId = "percent2",
      label = "Top Percentage to Show",
      min = 5,
      max = 50,
      value = 10
    ),
    strong("Experimentally determined binding sites on these proteins:"),
    verbatimTextOutput("truesites"),
    strong("Column Explanations:"),
    verbatimTextOutput("text"),
    tags$style(type = "text/css", "#text {white-space: pre-wrap;
                        word-break: keep-all;}"),
    strong("Input Descriptions:"),
    verbatimTextOutput("text2"),
    tags$style(type = "text/css", "#text2 {white-space: pre-wrap;
                        word-break: keep-all;}"),
    strong("Input Sources:"),
    tags$br(),
    tags$a(href = "https://www.phosphosite.org/homeAction.action", "PhosphoSitePlus", target =
             "_blank"),
    tags$br(),
    tags$a(href = "https://iupred2a.elte.hu/", "Iupred2A", target =
             "_blank"),
    tags$br(),
    tags$a(href = "https://www.compbio.dundee.ac.uk/1433pred", "14-3-3 Pred", target =
             "_blank")
  ),
  column(
    6,
    strong("Sites in Order of Probability (amino acid numbering)"),
    DT::dataTableOutput("result.comb2"),
    offset = 0
  )
))
    
# Define server logic
server <- function(input, output) {
    # Pull in data with APIs
    pull.data <- function() {
        # if in our data set, pull from there. Otherwise, get from APIs
        if(input$uniprot.ID %in% train.dat$uniprot.ID) {
            new.dat <- train.dat[train.dat$uniprot.ID==input$uniprot.ID,]
        } else {
            consensus <- GET(paste0("http://www.compbio.dundee.ac.uk/1433pred/1433pred_s",
                                    toupper(input$uniprot.ID),".json"))
            cdat <- fromJSON(rawToChar(consensus$content))
            sites <- cdat$Site
            # this has the iupred score
            iupred <- GET(paste0("http://iupred2a.elte.hu/iupred2a/",
                                   toupper(input$uniprot.ID),".json"))
            iupreddat <- fromJSON(rawToChar(iupred$content))
            seqs <- unlist(strsplit(iupreddat$sequence,''))
            amino <- seqs[sites]
            iusc <- iupreddat$iupred2[sites]
            
            # this has site, SVM, PSSM, ANN, Consensus
            # what is pSER/Thr
            
            #get normalization from new data, if no normalization score, set it to zero
            norm <- pred.dat$Percentages[pred.dat$ACC_ID==input$uniprot.ID]
            norm.site <- unlist(str_split(
                pred.dat$MOD_RSD[pred.dat$ACC_ID==input$uniprot.ID],pattern='-'))
            norm.site2 <- substr(norm.site[!(norm.site %in% c("","p"))],start=2,
                                 stop=nchar(norm.site[!(norm.site %in% c("","p"))]))
            grab.norm <- function(s) {if(s%in%norm.site2){norm[which(norm.site2==s)]}else{0}}
            norms <- unlist(sapply(sites,grab.norm))
            new.dat <- cdat
            new.dat$Amino.Acid <- as.factor(amino)
            new.dat$Iupred.score <- iusc
            new.dat$normalization <- norms
            new.dat$SVM <- as.numeric(new.dat$SVM)
            new.dat$PSSM <- as.numeric(new.dat$PSSM)
            new.dat$ANN <- as.numeric(new.dat$ANN)
            new.dat$Consensus <- as.numeric(new.dat$Consensus)
        }
        new.dat[new.dat$Amino.Acid %in% c("S","T"),]
    }
    ## only want S and Ts
    title <- function() {paste0("Top ",input$percent,"% most likely sites")}
    truesites1 <- function() {"True sites on this protein:"}
    truesites2 <- function() {
      prots <- unlist(strsplit(input$uniprot.ID.list,split=','))
      ret <- NA
      for(i in 1:length(prots)) {
        if(prots[i] %in% train.dat$uniprot.ID){
          ret <- c(ret,train.dat$Position1[train.dat$Response==1 & 
                                       train.dat$uniprot.ID==prots[i]],",")
        } else{ret <- c(ret,"Unknown",",")}
      }
      paste(ret[-c(1,length(ret))],collapse=' ')
      }
    output$logo <- renderImage
    output$tabletitle <- renderText({title()})
    output$truesites <- renderText(truesites2())
    output$result.comb2 <- DT::renderDataTable({ #dom='t' in list for just table
        datatable(result.comb2(),options = list(dom  = '<"top">t<"bottom">ip',pageLength=20,ordering=F))%>%formatStyle("Total Score",backgroundColor="lightblue")
    })
    output$text <- renderText({
        paste("Total Score: Random Forest model trained on inputs from phosphorylation, disorder, and 14-3-3-Pred scores.",
              "PTM & Disorder Score: Random Forest model trained on inputs from phosphorylation and disorder data.",
              "14-3-3-Pred Score: Random Forest model trained on 14-3-3-Pred inputs.",
              "Adapted 14-3-3-Pred Score: A rank order list of most probable binding sites found by applying cutoffs for consensus scores used by 14-3-3-Pred.",
              "*Note that sites shown are limited to only those amino acids with previously observed phosphorylations, as catalogued on phosphosite.org",
              "",sep="\n")
    })
    output$text2 <- renderText({
        paste("PTM values are derived from","HTP phosphorylation data gathered from phosphosite.org",
              "Disorder values are obtained from iupred2a.elte.hu",
              "14-3-3-Pred inputs used include the ANN, PSSM, and SVM output score values from compbio.dundee.ac.uk/1433pred",
              sep='\n')
    })
    output$protexs <- renderText({
        paste("Enter one protein ID such as:",
              "Q96PX8",
              "Or multiple with commas and no spaces, e.g.:",
              "Q96PX8,Q16695,P27986",sep='\n')
    })
    pull.data2 <- function() {
        prot.list <- unlist(strsplit(input$uniprot.ID.list,split=','))
        # if in our data set, pull from there. Otherwise, get from APIs
        for(i in 1:length(prot.list)) {
            if(prot.list[i] %in% train.dat$uniprot.ID) {
                new.dat <- train.dat[train.dat$uniprot.ID==prot.list[i],]
            } else {
                consensus <- GET(paste0("http://www.compbio.dundee.ac.uk/1433pred/1433pred_s",
                                        toupper(prot.list[i]),".json"))
                cdat <- fromJSON(rawToChar(consensus$content))
                sites <- cdat$Site
                # this has the iupred score
                iupred <- GET(paste0("http://iupred2a.elte.hu/iupred2a/",
                                     toupper(prot.list[i]),".json"))
                iupreddat <- fromJSON(rawToChar(iupred$content))
                seqs <- unlist(strsplit(iupreddat$sequence,''))
                amino <- seqs[sites]
                iusc <- iupreddat$iupred2[sites]
                
                # this has site, SVM, PSSM, ANN, Consensus
                # what is pSER/Thr
                
                #get normalization from new data, if no normalization score, set it to zero
                norm <- pred.dat$Percentages[pred.dat$ACC_ID==prot.list[i]]
                norm.site <- unlist(str_split(
                    pred.dat$MOD_RSD[pred.dat$ACC_ID==prot.list[i]],pattern='-'))
                norm.site2 <- substr(norm.site[!(norm.site %in% c("","p"))],start=2,
                                     stop=nchar(norm.site[!(norm.site %in% c("","p"))]))
                grab.norm <- function(s) {if(s%in%norm.site2){norm[which(norm.site2==s)]}else{0}}
                norms <- unlist(sapply(sites,grab.norm))
                new.dat <- cdat
                new.dat$Amino.Acid <- as.factor(amino)
                new.dat$Iupred.score <- iusc
                new.dat$normalization <- norms
                new.dat$SVM <- as.numeric(new.dat$SVM)
                new.dat$PSSM <- as.numeric(new.dat$PSSM)
                new.dat$ANN <- as.numeric(new.dat$ANN)
                new.dat$Consensus <- as.numeric(new.dat$Consensus)
                new.dat$uniprot.ID <- rep(prot.list[i],nrow(new.dat))
            }
          new.dat2 <- new.dat[new.dat$Amino.Acid %in% c("S","T"),]
          if(i==1){
            all.new.dat <- new.dat2
          } else {
            #all.new.dat <- rbind(all.new.dat,new.dat2)
            all.new.dat <- smartbind(all.new.dat,new.dat2)
          }
        } 
        all.new.dat
    }
    result.comb2 <- function() {
        new.dat <- pull.data2()
        new.dat$comb.id <- paste(new.dat$uniprot.ID,new.dat$Site,sep=' : ')
        these.sites <- new.dat$Site
        #number of total candidate sites on that protein
        total.sites <- length(these.sites)
        top.no <- ceiling(input$percent2/100*total.sites)
        # just predict directly instead of through other functions
        pred.full <- predict(full,new.dat,"prob")
        pred.red <- predict(red,new.dat,"prob")
        pred.con <- predict(con,new.dat,"prob")
        copy <- site.copy(new.dat) # already sorted
        # sort each table by their probability of 1, output the combined ID
        fullsort <- new.dat[,'comb.id'][order(pred.full[,2],decreasing=TRUE)][1:top.no]
        redsort <- new.dat[,'comb.id'][order(pred.red[,2],decreasing=TRUE)][1:top.no]
        consort <- new.dat[,'comb.id'][order(pred.con[,2],decreasing=TRUE)][1:top.no]
        copysort <- copy$comb.id[1:top.no]
        tab1 <- cbind(fullsort,redsort,consort,copysort)
        # output
        colnames(tab1) <- c("Total Score","PTM & Disorder Score",
                            "14-3-3-Pred Score","Adapted 14-3-3-Pred Score")
        tab1
    }
        
}

# Run the application 
shinyApp(ui = ui, server = server)
