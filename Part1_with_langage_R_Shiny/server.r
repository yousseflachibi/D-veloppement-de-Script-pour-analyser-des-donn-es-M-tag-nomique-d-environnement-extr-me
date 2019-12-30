options(shiny.maxRequestSize=1024*1024^2) 
library(shiny)
library(shinydashboard)
if (!requireNamespace("rentrez", quietly=TRUE)){
  install.packages("rentrez")
  install.packages("stringr")
  install.packages("filesstrings")
  install.packages("CHNOSZ")
  install.packages("tidyverse")
  install.packages("Biostrings")
  install.packages("seqinr")
  install.packages("stringr")
  install.packages("rex")
  install.packages("phylotools")
  install.packages("qdapRegex")
  install.packages("openxlsx")}
library("rentrez")
library("stringr")
library("filesstrings")
library("CHNOSZ")
library("tidyverse")
library("Biostrings")
library("seqinr")
library("stringr")
library("rex")
library("rex")
library("phylotools")
library("qdapRegex") 
library("openxlsx")

shinyServer(function(input, output, session){
  
  
  
  destDir <-paste(dirname(rstudioapi::getSourceEditorContext()$path),"www", sep="/")
  test_50 <<-data.frame()
  test_75 <<-data.frame()
  observeEvent(input$go_1, {
    shinyjs::disable("go_1")
    shinyjs::show("text_1")
    withProgress(message = 'Traitement Searching', value = 0, {
      
      keywords <- input$input_1
      parametre <- input$input_2
      r_search <- entrez_search(db="genome",
                                term=paste(keywords,"[ORGN] AND ",parametre),
                                retmax=15)
      r_search$ids
      n <- length(r_search$ids)
      cc<-""
      for (i in 1:n){
        tryCatch({
          linked_seq_ids <- entrez_link(dbfrom="genome", id=r_search$ids[i], db="nuccore")
          linked_transripts <- linked_seq_ids$links$genome_nuccore
          all_recs <- entrez_fetch(db="nuccore", id=linked_transripts, rettype="fasta_cds_na")
          temp<-paste(r_search$ids[i],".fna",sep = "")
          temp <- paste(dirname(rstudioapi::getSourceEditorContext()$path),temp,sep = "/")
          write(all_recs, temp)
          file.move(temp,destDir, overwrite = TRUE)
        }, error=function(e){})
        incProgress(1/n, detail = paste("Doing part", i))
      }
      
      #file_names<-list.files(path = "C:/Users/dell/Documents/R/Projet/MedBio/www" , pattern="*.fna")
      #for (d in 1:length(file_names)){
      #  print("***")
      #  print(d)
      #  del_file <- readDNAStringSet(paste("C:/Users/dell/Documents/R/Projet/MedBio/www",file_names[d],sep="/"))
      #  if(length(del_file)==0){
      #    print("**********")
      #    #unlink(paste("C:/Users/dell/Documents/R/Projet/MedBio/www",file_names[d],sep="/"), recursive = TRUE)
      #    file.remove(paste("C:/Users/dell/Documents/R/Projet/MedBio/www",file_names[d],sep="/"))
      #    print("**********")
      #    }
      #}
      shell.exec(destDir)
      shinyjs::enable("go_1")
      shinyjs::hide("text_1")
      updateTextInput(session, inputId="input_1", label = "Organisme : ", value ="")
      updateTextInput(session, inputId="input_2", label = "Parametre : ", value ="")
    })
    
  })
    observeEvent(input$go_2,{
      shinyjs::disable("go_2")
      shinyjs::show("text_2")
      
      
      ###################################################
      #### Extract protein_id from all fasta files   ####
      ###################################################
      prot <- ""
      rm(df)
      df<-data.frame()
      ddf<-data.frame()
      withProgress(message = 'Traitement !!', value = 0, {
        n <- length(input$input_3$datapath)+2
        #for (j in 1:length(input$input_3$datapath)){
        
        for (j in 1:length(input$input_3$datapath)){
          tryCatch({
            print(j)
            fichier_num<-j
            #r_file <- readDNAStringSet(paste("C:\\Users\\dell\\Documents\\R\\Projet\\MedBio\\www\\my_transcripts_",j,".fna",sep = ""))
            r_file <- readDNAStringSet(input$input_3$datapath[fichier_num])
            seq_name = names(r_file)
            prot<-re_matches(seq_name,rex(capture("protein=",anything,"[p")))
            prot<-qdapRegex::ex_between(prot, "protein=", "]")[[1]]
            f<-unique(data.frame(prot))
            prot_2<-rep(input$input_3$name[fichier_num],nrow(f))
            ddf<- rbind(ddf,data.frame(prot_2))
            df <-rbind(df,unique(data.frame(prot,input$input_3$name[fichier_num])))
          }, error=function(e){})
          incProgress(1/n, detail = paste("Doing part", j))
        }
        #Calculating frequences
        n_occur2 <- data.frame(table(Nom_Protein = c(df[,1])))#on aura pas le column num fichier
        output$distPlot <- renderPlot({
          counts <- table(df[,1])
          barplot(counts, main="Test ...", horiz=FALSE)
        })
        incProgress(1/n, detail = paste("Doing part", length(input$input_3$datapath)+1))
        #Extraire la liste des protéines qui apparait dans +50% des fichiers.
        test_50<<- n_occur2[(which((n_occur2[,2]*100 /length(input$input_3$datapath))>50)),]
        #Extraire la liste des protéines qui apparait dans +75% des fichiers.
        test_75<<- n_occur2[(which((n_occur2[,2]*100 /length(input$input_3$datapath))>75)),]
        #Extraire la liste des protéines qui apparait dans +75% des fichiers.
        test_90<<- n_occur2[(which((n_occur2[,2]*100 /length(input$input_3$datapath))>90)),]
        #Extraire la liste des protéines qui apparait dans -50% des fichiers.
        test_moins_50<<- n_occur2[(which((n_occur2[,2]*100 /length(input$input_3$datapath))<=50)),]
        incProgress(1/n, detail = paste("Doing part", length(input$input_3$datapath)+2))
        
      })
      
      
      
      #######################################################
      ##### supprimer seq proteins no commun entre file #####
      #######################################################
      
      for (ss in 1:length(input$input_3$datapath)){
        tryCatch({
          print(ss)
          updateProgressBar(session = session,id = "pb1",value = ((ss*100)/length(input$input_3$datapath)))
          fichier_num<-ss
          r_file <- readDNAStringSet(input$input_3$datapath[fichier_num])
          seq_name = data.frame(r_file,names(r_file))
          id_no_sup<- paste(test_moins_50[,2],sep="\t")
          rm.sequence.fasta(paste(input$input_3$datapath[fichier_num])
                            , outfile = paste(destDir,"/sequence_removed_Pro_commun/sequence.removed_",fichier_num,".fna",sep = "")
                            , to.rm = id_no_sup)
        }, error=function(e){})
      }
      
      shell.exec(paste(destDir,"/sequence_removed_Pro_commun",sep = ""))
      
      
      
      #######################################################
      ######     create file csv of protein commun     ######
      #######################################################
      
      for (file_N in 1:3){
        updateProgressBar(session = session,id = "pb2",value = ((file_N*100)/3))
        if(file_N==1){
          write.xlsx(test_50,paste(destDir,"/protein_commun/protein_commun_plus_50%.xlsx",sep = ""), sheetName="Sheet1",
                     col.names=TRUE, row.names=TRUE, append=FALSE, showNA=TRUE, password=NULL)
        }else if(file_N==2){
          write.xlsx(test_75,paste(destDir,"/protein_commun/protein_commun_plus_75%.xlsx",sep = ""), sheetName="Sheet1",
                     col.names=TRUE, row.names=TRUE, append=FALSE, showNA=TRUE, password=NULL)
        }else{
          write.xlsx(test_90,paste(destDir,"/protein_commun/protein_commun_plus_90%.xlsx",sep = ""), sheetName="Sheet1",
                     col.names=TRUE, row.names=TRUE, append=FALSE, showNA=TRUE, password=NULL)
        }
      }
      shell.exec(paste(destDir,"/protein_commun",sep = ""))
      
      shinyjs::enable("go_2")
        shinyjs::hide("text_2")
        shinyjs::enable("go_3")
      })
      
      
      
    
  })
  
  
  
