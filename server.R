rm(list=ls())
library("readxl")
library("tidyverse")
library("visNetwork")
library("shiny")
library("shinydashboard")
library("DT")
library("dplyr")

# define a shiny server
server <- function(input, output, session) {
  options(shiny.maxRequestSize=30*1024^2) 
  
  source("functions.R",local = TRUE)
  ##########################################################################
  # Read Excel file and displays its sheets
  ## list all excel sheets 
  sheets <- excel_sheets(path = "NETWORKS.xlsx")
  sheets_ok <- c()
  sheets_notok <- c()
  sheets_flag <- c()
  # loop through all the sheets for four columns
  must_have_columns <- c("TRAITID1","TRAITID2","PVALUE","BETA", "COR")
  annotation_columns <- c("TRAITID","SHORTNAME","PLAT")
  # i=1
  n=0
  m=0
  # or pre-allocate for slightly more efficiency
  datalist = vector("list", length = length(sheets))
  i=1
  for (i in 1:length(sheets)) {
    sheetContent = read_excel(path = "NETWORKS.xlsx", sheet = sheets[i], .name_repair = ~make.unique(.x, sep = "_"))
    musthave_cols <- unlist(lapply(must_have_columns, function(header) {
      grep(paste0("^",header,"$"), names(sheetContent), ignore.case = TRUE)
    }))
    annotation_cols <- unlist(lapply(annotation_columns, function(header) {
      grep(paste("^",header,"$",sep = ""), names(sheetContent), ignore.case = TRUE)
    }))
    if(length(musthave_cols) == 4){
      m=m+1
      idx <- grep("^COR$|^BETA$", names(sheetContent))
      dat <- data.frame(sheetContent %>% dplyr::select(musthave_cols, idx[1]))
      
      dat$type <- rep(sheets[i])
      colnames(dat) <- c("to","from","pvalue","weight","type")
      dat$i <- i  # maybe you want to keep track of which iteration produced it?
      dat$id <- paste0(sheets[i],"_0",rownames(dat))
      datalist[[i]] <- dat # add it to your list
      sheets_ok[m] <- sheets[i]
    } else if(length(annotation_cols) == 3) {
      idx <- grep("^PLAT$", names(sheetContent))
      anno <- sheetContent %>% dplyr::select(all_of(annotation_cols))
      
    } else {
      n=n+1
      columns_notfound <- toString(must_have_columns[!(must_have_columns %in% names(sheetContent))])
      sheets_flag[n] <- paste(sheets[i],": missing columns - ", columns_notfound)
      sheets_notok[n] <- sheets[i]
    }
  }
  
  datalist = datalist[-which(sapply(datalist, is.null))]
  
  all_edges = do.call(rbind, datalist)
  all_edges$sign = sign(all_edges$weight)
  all_edges$weight = abs(all_edges$weight)
  
  all_nodes = data.frame(TRAITID = unique(sort(c(all_edges$from, all_edges$to))))
  # get node annotations
  all_nodes = left_join(all_nodes, anno)
  #all_nodes$SHORTNAME <- paste0(all_nodes$PLAT,":",all_nodes$TRAITNAME)
  
  # replace all_edges$to and $from ids with all_nodes$SHORTNAME
  TRAITID2SHORTNAME = all_nodes$SHORTNAME
  names(TRAITID2SHORTNAME) = all_nodes$TRAITID
  
  SHORTNAME2TRAITID = all_nodes$TRAITID
  names(SHORTNAME2TRAITID) = all_nodes$SHORTNAME
  
  all_edges$from = TRAITID2SHORTNAME[all_edges$from]
  all_edges$to = TRAITID2SHORTNAME[all_edges$to]
  all_nodes$id = all_nodes$SHORTNAME
  names(all_nodes)[which(names(all_nodes) == "PLAT")] = "plat"
  
  # round the p-values and weights
  all_edges$pvalue = signif(all_edges$pvalue, digits = 2)
  all_edges$weight = signif(all_edges$weight, digits = 2)
  ######################################################################################################################
  #########################################################
  # a data structure for a network
  #########################################################
  
  fullnet = list(
    edges = all_edges,
    nodes = all_nodes
  )
  
  ##########################################################################
  # add colors for platforms to network
  # to be used as: PLATcols[fullnet$nodes$plat]
  # 
  
  PLATlist = unique(fullnet$nodes$plat)
  
  PLATcols = rep("#999999",length(PLATlist))
  names(PLATcols) = PLATlist
  
  
  # num_plats <- length(PLATlist)
  # colors <- rainbow(num_plats)
  # PLATcols = colors
  # names(PLATcols) = PLATlist
  
  
  
  PLATshapes = rep("star",length(PLATlist))
  # PLATshapes = rep("icon",length(PLATlist))
  names(PLATshapes) = PLATlist
  
  # # define PLAT colors manually
  PLATcols["DNA"] = "#23bbee"
  PLATcols["SOMA"] = "#a62281"
  PLATcols["BRAIN"] = "#f2921f"
  PLATcols["BM"] = "#ffc815"
  PLATcols["LD"] = "#ffc815"
  PLATcols["CPG"] = "#145da9"
  PLATcols["CLIN"] = "#a0b6a8"
  PLATcols["CM"] = "#57ba47"
  PLATcols["IgA"] = "#e41d30"
  PLATcols["RNA"] = "#5c2d83"
  PLATcols["HD4"] = "#57ba47"
  PLATcols["IgG"] = "#e41d30"
  PLATcols["miRNA"] = "#5c2d83"
  PLATcols["OLINK"] = "#a62281"
  PLATcols["PGP"] = "#e41d30"
  PLATcols["PM"] = "#57ba47"
  PLATcols["SM"] = "#57ba47"
  PLATcols["UM"] = "#57ba47"
  PLATcols["STAT"] = "#EEEEEE"
  PLATcols["GWAS"] = "#EEEEEE"
  
  # define PLAT colors manually
  PLATshapes["UM"] = "square"
  PLATshapes["CM"] = "square"
  PLATshapes["SM"] = "triangle"
  PLATshapes["DNA"] = "diamond"
  PLATshapes["RNA"] = "diamond"
  PLATshapes["CPG"] = "diamond"
  
  PLATshapes["STAT"] = "circle"
  PLATshapes["GWAS"] = "dot"
  
  fullnet$nodes$color = PLATcols[fullnet$nodes$plat]
  
  fullnet$nodes$group = fullnet$nodes$plat
  
  fullnet$nodes$shape = PLATshapes[fullnet$nodes$plat]
  
  fullnet$edges$title = paste(fullnet$edges$type, ": ", fullnet$edges$id,
                              ", p=", fullnet$edges$pvalue,
                              ", beta=",
                              ifelse(fullnet$edges$sign >0, "", "-"), fullnet$edges$weight,
                              sep = "")
  fullnet$edges$color = ifelse(fullnet$edges$sign > 0, "blue", "red")
  
  #convert all pvalues == 0 to a very small number so -log10 won't throw an error
  fullnet$edges$pvalue[fullnet$edges$pvalue == 0] <- 1E-100
  
  edge_thickness = -log10(fullnet$edges$pvalue)
  minVal = min(edge_thickness)-1
  maxVal = max(edge_thickness)-1
  numItv = 10
  seq_range <- round(seq(minVal, maxVal, length.out=numItv))
  
  fullnet$edges$width = cut(edge_thickness, seq_range, labels = F)
  
  # list of ids to choose from
  plat_list = c("ALL", fullnet$nodes$plat %>% unique %>% sort())
  
  # limit the number of nodes here
  max_nodes_list = c(1,20,40,60,80,100,150,200)
  
  ##########################################################################
  ##########################################################################
  
  # storage for reactiveValues
  storage <- reactiveValues(
    focus = "CLIN: HbA1c (%)",
    plat = "CLIN"
  )

  # storage <- reactiveValues(
  #   focus = all_nodes$SHORTNAME[3],
  #   plat = PLATlist[2]
  # )
  
  
  # process uploaded excel file
  observeEvent(input$file_in,{
    req(input$file_in)
    name <- input$file_in
    #file.copy( from = input$file_in$datapath, to = "NETWORKS.xlsx", overwrite=TRUE)
    inFile <- input$file1
    print(inFile)
  })
  
  # set focus if pull-down menu changes
  observeEvent(input$trait,
               {
                 # below is a quick-fix to avoid overwriting of the query string by the first menu value .. not very elegant
                 if ((input$trait != "")&(input$trait != "BM: Arg")) {
                   storage$focus <- input$trait
                   # message("focus set by dropdown to ", storage$focus)
                 }
                 
               })
  
  # # read the URL parameter from session$clientData$url_search
  observe({
    query <- parseQueryString(session$clientData$url_search)
    if (!is.null(query[['focus']])) {
      storage$focus = query[['focus']]
    }
    if (!is.null(query[['maxnodes']])) {
      updateSelectInput(session, "maxnodes", selected = query[['maxnodes']])
    }
  })
  
  #this section is moved to line 631  
  # select trait pair from association tables  
  # output[[paste0(sheets_ok[i],"table")]]
  #observeEvent(input$GWAStable_rows_selected,{
  # observeEvent(input[[paste0(sheets_ok[1],"table_rows_selected")]],{
  #   print("line433")
  #   print(input[[paste0(sheets_ok[1],"table_rows_selected")]])
  #   #storage$asso_selected = GWAS$ID[as.numeric(input$GWAStable_rows_selected)]
  #   #storage$asso_selected = paste0(sheets_ok[1])$ID[as.numeric(input$GWAStable_rows_selected)]
  #   print(head(datalist[1]))
  #   #storage$asso_selected = paste0(sheets_ok[1])$id[as.numeric(input[[paste0(sheets_ok[1],"table_rows_selected")]])]
  #   })
  # observeEvent(input$EWAStable_rows_selected,{storage$asso_selected = EWAS$ID[as.numeric(input$EWAStable_rows_selected)]})
  # observeEvent(input$RWAStable_rows_selected,{storage$asso_selected = RWAS$ID[as.numeric(input$RWAStable_rows_selected)]})
  # observeEvent(input$MBHtable_rows_selected,{storage$asso_selected = MBH$ID[as.numeric(input$MBHtable_rows_selected)]})
  # observeEvent(input$GGMtable_rows_selected,{storage$asso_selected = GGM$ID[as.numeric(input$GGMtable_rows_selected)]})
  # observeEvent(input$GENOtable_rows_selected,{storage$asso_selected = GENO$ID[as.numeric(input$GENOtable_rows_selected)]})
  # observeEvent(input$STATtable_rows_selected,{storage$asso_selected = STAT$ID[as.numeric(input$STATtable_rows_selected)]})
  # observeEvent(input$CATAtable_rows_selected,{storage$asso_selected = CATA$ID[as.numeric(input$CATAtable_rows_selected)]})
  
  # set focus if pull-down menu changes
  observeEvent(storage$asso_selected,
               {
                 message("storage$asso_selected:", storage$asso_selected)
                 
                 if (substr(storage$asso_selected,1,4) == "STAT") {
                   # special case for STAT
                   ix = which(STAT$ID == storage$asso_selected)   
                   storage$selected_TRAIT1 = TRAITID2SHORTNAME[STAT$TRAITID[ix]]
                   storage$selected_TRAIT2 = TRAITID2SHORTNAME[STAT$TRAITID[ix]]
                   if (is.na(storage$selected_TRAIT1)){
                     storage$selected_TRAIT1 = ""
                     storage$selected_TRAIT2 = ""
                   }
                 } else {
                   
                   ix = which(all_edges$id == "MBH_02140")
                   ix = which(all_edges$id == storage$asso_selected)
                   storage$selected_TRAIT1 = all_edges$from[ix]
                   storage$selected_TRAIT2 = all_edges$to[ix]
                 }
                 message("storage$selected_TRAIT1:", storage$selected_TRAIT1)
                 message("storage$selected_TRAIT2:", storage$selected_TRAIT2)
               }
  )             
  
  
  # set focus if network is clicked
  observeEvent(input$network_selected,
               {
                 if (input$network_selected != "") {
                   storage$focus <- input$network_selected 
                   message("focus set by network click to ", storage$focus)
                 }
               })
  
  # observe click onto SubmitPair
  observeEvent(input$submitInfo, {
    updateTabItems(session = session, inputId = "tabs", selected = "Tab1")
    storage$focus = "special:pair"
    message("focus set by table click to ", storage$focus)
  })
  
  # text output of focus
  output$focus <- renderText({
    paste("Click here to focus on", storage$focus)
  })
  
  # text output of asso_selected
  output$asso_selected <- renderText({
    if (length(storage$selected_TRAIT1) == 0) {
      paste("Select an association, then click here to view the network context")
    } else if (storage$selected_TRAIT1 != "") {
      paste("Click here to focus network on association ", storage$asso_selected)
    } else {
      paste("Sorry - no data for association ", storage$asso_selected)
    }   
  })
  
  output$SelectTrait = renderUI({
    
    if (input$plat == "ALL") {
      traitlist = c(fullnet$nodes$id %>% sort())
    } else {
      ix = which(fullnet$nodes$plat == input$plat)
      traitlist = c(fullnet$nodes$id[ix] %>% sort())
    }
    
    tagList(
      selectInput(inputId = "trait", label = "select trait", 
                  choices = traitlist, selected = traitlist[1])
    )
  })
  
  output$network <- renderVisNetwork({
    
    maxnodes = input$maxnodes
    if (length(maxnodes) == 0) {maxnodes = max_nodes_list[2]}
    maxnodes = as.numeric(maxnodes)
    message("limiting to maxnodes = ", maxnodes)
    message("GoButton pressed: ", input$GoButton)
    message("submitInfo pressed: ", input$submitInfo)
    message("storage$focus:", isolate(storage$focus))
    
    # to fix some weird bug
    zwi = isolate(storage$focus)
    zwi = zwi[1]
    
    if (length(zwi) == 0) {
      message("storage$focus empty, returning")
      act_trait = ""
      return()
    } else if (zwi == "") {
      message("storage$focus empty, returning")
      act_trait = ""
      return()
    } else if (zwi == "special:pair") {
      message("selected an edge")
      act_trait = c(isolate(storage$selected_TRAIT1), isolate(storage$selected_TRAIT2))
      # storage$focus = storage$selected_TRAIT1 # set focus to the first trait of the edge
      storage$focus = c(isolate(storage$selected_TRAIT1), isolate(storage$selected_TRAIT2))
    } else {
      act_trait = isolate(storage$focus)
    }
    
    message("getting subnet for act_trait: '", act_trait, "'")
    
    # get a network that is tractable
    # subnet = neighbors(act_trait, fullnet) %>% maxneighbors(fullnet, limit = maxnodes) %>% nodes2network(fullnet)
    subnet = neighbors(act_trait, fullnet) %>% maxneighbors_noSTAT(fullnet, limit = maxnodes) %>% nodes2network(fullnet)
    # summary_network(subnet)
    
    # store in a reactive variable
    storage$subnet = subnet
    
    edges = data.frame(from = subnet$edges$from, 
                       to = subnet$edges$to,
                       title = subnet$edges$title,
                       color = subnet$edges$color,
                       width = subnet$edges$width)
    nodes = data.frame(id = subnet$nodes$id, 
                       label = subnet$nodes$id,
                       title = subnet$nodes$id,
                       color = subnet$nodes$color,
                       shape = subnet$nodes$shape)
    #group = subnet$nodes$group
    
    message("rendering subnet with ", length(nodes$id), " nodes and ", length(edges$from), " edges for act_trait: ", act_trait)
    
    visNetwork(nodes, edges, height = "500px", width = "100%") %>% 
      visNodes(shadow = list(enabled = TRUE, size = 10), 
               scaling = list(min=10,max=30)) %>%
      visLayout(randomSeed = 4711) %>% 
      visOptions(nodesIdSelection = list(enabled = TRUE, style = 'width: 1px; height: 1px;')) %>% 
      visPhysics(stabilization = FALSE) %>% 
      visEdges(smooth = TRUE)
    
  })
  
  # export the current edges as a table
  output$EXPORTtable = DT::renderDataTable({
    
    message("UPDATING EXPORTtable")
    zwi = storage$subnet$edges
    
    out = data.frame(ID = zwi$id,
                     TYPE = zwi$type,
                     TRAITID1 = SHORTNAME2TRAITID[zwi$from],
                     TRAIT1 = zwi$from,
                     BETA = zwi$sign*zwi$weight,
                     PVALUE = zwi$pvalue,
                     TRAITID2 = SHORTNAME2TRAITID[zwi$to],
                     TRAIT2 = zwi$to
    )
    out
  }, filter = 'top', rownames = FALSE, selection = 'single', server = FALSE, extensions = 'Buttons', options = list(pageLength = 10, scrollX = TRUE, dom = 'Blfrtip', 
                                                                                                                    buttons = 
                                                                                                                      list("copy", list(
                                                                                                                        extend = "collection"
                                                                                                                        , buttons = c("csv", "excel", "pdf")
                                                                                                                        , text = "Download"
                                                                                                                      ) ) # end of buttons customization
  ))
  
  
  ############# 
  ############# 
  output$tabs = renderUI({
    nTabs = length(datalist)
    myTabs = lapply(seq_len(nTabs), function(i) {
      tabPanel(paste0(sheets_ok[i]),width =32,
               #DT::dataTableOutput(paste0("datatable_",i))
               DT::dataTableOutput(paste0(sheets_ok[i],"table"))
      )
    })
    do.call(tabBox, c(myTabs, width=32))
  })
  
  lapply(seq_len(length(datalist)), function(i){
    output[[paste0(sheets_ok[i],"table")]] <- DT::renderDataTable({
      data.frame(datalist[i]
      )}, filter = 'top', rownames = FALSE, selection = 'single', 
      options = list(pageLength = 10, scrollX = TRUE))
  })
  
  lapply(seq_len(length(sheets_ok)), function(tables){
    observeEvent(input[[paste0(sheets_ok[tables],"table_rows_selected")]],{
      #print(input[[paste0(sheets_ok[tables],"table_rows_selected")]])
      #print(head(datalist[[tables]]))
      storage$asso_selected = datalist[[tables]]$id[as.numeric(input[[paste0(sheets_ok[tables],"table_rows_selected")]])]
      #print (datalist[[tables]]$ID[as.numeric(input[[paste0(sheets_ok[tables],"table_rows_selected")]])])
    })
  })
  output$SAMPLEtable = DT::renderDataTable(({
    data.frame(sheets_flag)
  }))
  
}
