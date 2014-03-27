# v.5. allow different test instead of ANOVA, remove ecology plot, 
# remove stats based on cover - this can be discriminated by post-hoc tests when zone differences found.
# v.6. 

library("pgirmess")

sig_code <- function(x_in){
  sig_code_vector <- c()
  for(x in x_in) {
    sig_code = ""
    if(x<=0.001) {sig_code="***"} else {
      if(x<=0.01) {sig_code="**"} else {
        if(x<=0.05) {sig_code="*"} else {
          if(x<=0.1) {sig_code="."}
        }
      }
    }
    sig_code_vector <- c(sig_code_vector,sig_code)
  }
  return(sig_code_vector)
}

sig_code2 <- function(x_in,sym){
  sig_code_list <- array()
  i=0
  for(x in x_in) {
    i <- i+1
    sig_code = ""
    if(x<=0.001) {sig_code=paste(rep(sym[i],3),collapse="")} else {
      if(x<=0.01) {sig_code=paste(rep(sym[i],2),collapse="")} else {
        if(x<=0.05) {sig_code=sym[i]} 
      }
    }
    sig_code_list[i] <- sig_code
  }
  return(sig_code_list)
}
#t <-  sig_code2(c(0.0005),c("+"))
#t

if(FALSE){   # for testing

load("Q:/sync/lab_book/kalahari/files/pyrosequencing/R/R_objects/expt_ord.rda")
my_phyloseq <- expt #.filtered # for testing
topN <- 9 # for testing
verbose <- FALSE # for testing
plot_title <- "Test" # for testing
rename_facets <- "FALSE" # for testing
}

my_OTU_plot <- function(my_phyloseq,topN=9,verbose=FALSE,plot_title="OTU_plot",rename_facets="no") {
  stats_report_file <- paste(STATS_dir,"/stat_tables_for_",plot_title,".txt",sep="")
  writefile <- paste(tables_dir,"/OTUs_",plot_title,".csv",sep="")
  sink(stats_report_file,split=TRUE)
    cat("report for:",plot_title,"\n")
    n_input <- ntaxa(my_phyloseq)
    cat(plot_title,":","Using the top",topN,"of",n_input,"taxa in group",plot_title,"\n\n",sep=" ")
    #cat(plot_title,":","Sum of observations (for whole group, not just the plotted ones), average per sample:", sum(taxa_sums(my_phyloseq))/nsamples(expt.filtered),"\n")
    print(my_phyloseq)
    cat("\nA list of p-values corrected for multiple comparisons using the FDR method can be found in file ",writefile,"\n",sep="")
   # cat("A value of >0.1 is desirable for the Shapiro test of normality, however this in itself does not prove a normal distribution. In many cases the ANOVA results can not be used because the residuals are not normal, therefore we prefer the Kruskal-Wallis test for comparing OTU abundance amongst categories.\n\n")
  
  carbon_plot_file <- paste(extra_figures_dir,"/",plot_title,"_carbon.png",sep="")
  zone_plot_file <- paste(extra_figures_dir,"/",plot_title,"_zone.png",sep="")
  OTU_tree_file <- paste(extra_figures_dir,"/",plot_title,"_tree.png",sep="")
  
  # keep only top N taxa
  TopNOTUs <- names(sort(taxa_sums(my_phyloseq), TRUE)[1:topN])
  my_phyloseq <- prune_taxa(TopNOTUs, my_phyloseq)
  
  if(rename_facets=="Phylum"){
    taxa_names(my_phyloseq) <- tax_table(my_phyloseq)[,"Phylum"]
    taxa_names(my_phyloseq) <- gsub(taxa_names(my_phyloseq),pattern="-",replacement=".")
   } 
  if(rename_facets=="Class"){
    taxa_names(my_phyloseq) <- tax_table(my_phyloseq)[,"Class"]
    taxa_names(my_phyloseq) <- gsub(taxa_names(my_phyloseq),pattern="-",replacement=".")
  } 
  
  # get the OTU tables
  OTUs.table <- data.frame(otu_table(my_phyloseq))
  
  OTUs.t <- t(OTUs.table)
  
  sample_data(my_phyloseq)$zone_type <- paste(sample_data(my_phyloseq)$zone_code,sample_data(my_phyloseq)$type,sep="_")
  
  sample_data(my_phyloseq)$sample <- row.names(sample_data(my_phyloseq))
  sample_data <- sample_data(my_phyloseq)[,c("sample","zone","zone_code","zone2line","cover","type","depth2","month")]
  sample_data.extended <- sample_data(my_phyloseq)[,c("sample","zone","zone_code","zone2line","cover","type","month","Carbon","Q10_soil")]
  #,"Carbon","Q10_soil"
  # combine OTU and sample data 
  OTUs_and_env_data <- data.frame(sample_data,OTUs.t)
  OTUs_and_env_data.extended <- data.frame(sample_data.extended,OTUs.t)
  if(verbose) {
    cat("\nTable of OTU abundance per sample, with sample metadata:\n")
    print(OTUs_and_env_data)
  }
  # Put dataframes in order of abundance
  # Already in order. Leave alone. #OTUs_and_env_data <- OTUs_and_env_data[order(rowSums(OTUs.t),decreasing=TRUE) , ]
  # colSums(OTUs_and_env_data[7:15])
  # melt the data for faceting
  OTUs_and_env_data.melt <- melt(OTUs_and_env_data)
  
  # get some extra variables into the melted dataframe
  # Weird way of doing it but it's because the OTU data needs to be the only numerical data during melting
  index <- match(OTUs_and_env_data.melt$sample, OTUs_and_env_data.extended$sample)
  #OTUs_and_env_data.melt$sample_check <- OTUs_and_env_data.extended[index,"sample"]
  OTUs_and_env_data.melt$Carbon <- OTUs_and_env_data.extended[index,"Carbon"]
  #OTUs_and_env_data.melt$Q10_soil <- OTUs_and_env_data.extended[index,"Q10_soil"]
  
  #identical(OTUs_and_env_data.melt$sample,OTUs_and_env_data.melt$sample_check)
  
  # re-name the OTU's that had names changed by melting
  OTUs_and_env_data.melt$variable <- gsub(pattern="X",replacement="",x=OTUs_and_env_data.melt$variable)
  # plot the OTUs
  
  # OTU_names <- row.names(tax_table(my_phyloseq))
  OTU_names <- row.names(OTUs.table)

  my_phyloseq.merged <- merge_samples(my_phyloseq,"zone_type") 
  
  # These counts need to be divided by 6 manually because Using the phyloseq method doesn't work as expected. The division is required because when merging samples they are "summed", not "averaged"
  OTU_prob_table <- as.data.frame(t(otu_table(my_phyloseq.merged)))/6
  colnames2 <- names(OTU_prob_table)
  colnames <- names(as.data.frame(tax_table(my_phyloseq)))
  # prepare columns to hold the taxonomy data
  OTU_prob_table[,colnames] <- "x"
  # re-ordering the columns for easier reading
  OTU_prob_table <- OTU_prob_table[,c(colnames,colnames2)]

 # post_hoc_list <- list()
  OTU_name <- OTU_names[1] # for testing
  for (OTU_name in OTU_names) { 
    cat("\n",rep("-",30),"\n")
    cat(rep(":",10),"OTU ",OTU_name,rep(":",10),"\n")
    cat(rep("-",30),"\n\n")

    OTU_name <- gsub(OTU_name,pattern="-",replacement=".")
    this_OTU <- OTUs_and_env_data.melt[OTUs_and_env_data.melt$variable==OTU_name,]

    # Stats to compare abundance by factors
    # First set defaults
    stat_result_type <- NA
    shapiro.test <- NA
    stat_result_zone <- NA
    stat_result_month <- NA
    
    # Now check that we have factors to compare and then do the tests
    # The factors of interest in the paper are zone and type
    test_zone <-  length(unique(sample_data(my_phyloseq)$zone)) >1
    test_type <- length(unique(sample_data(my_phyloseq)$type)) >1

    # If we have both factors in the data then it is best to compare them together in the same model:
    if(test_zone && test_type) {
      # cat ("### ANOVA ###\n")
      aov_zt <- aov(this_OTU$value ~ factor(this_OTU$type) + factor(this_OTU$zone) )
       anova_zt <- anova(aov_zt) 
      anova_p_zone <- anova_zt$'Pr(>F)'[1]
      anova_p_type <- anova_zt$'Pr(>F)'[2]

      # check if the residuals on the anova are normal 
      # (other checks were also performed separately)
      # A p-value above 0.1 is desirable but not in itself proof of normality
       shapiro.test <- shapiro.test(aov_zt$residuals) 
      # qqnorm(aov_zt$residuals)
      
      # Perform a post-hoc test if significant difference found:
      if(anova_p_type < 0.05 | anova_p_zone <0.05 ) {
        # cat ( "ANOVA significant so doing post-hoc test:\n")
        tukey <- TukeyHSD(x=aov_zt,ordered=TRUE) 
      }
      
      # same again using Kruskal Wallis non parametric alternative
      cat ("### Kruskal Wallis ###\n")
      # kruskal_zt <- kruskal.test(this_OTU$value ~ factor(this_OTU$type) + factor(this_OTU$zone))
      # the above gives only 1 p-value, so we have to test separately for zone and type:
      print ( kruskal_z <- kruskal.test(this_OTU$value ~ factor(this_OTU$zone)) )
      print ( kruskal_t <- kruskal.test(this_OTU$value ~ factor(this_OTU$type)) )
      print ( kruskal_m <- kruskal.test(this_OTU$value ~ factor(this_OTU$month)) )
      kruskal_p_zone <- kruskal_z$p.value
      kruskal_p_type <- kruskal_t$p.value
      kruskal_p_month <- kruskal_m$p.value
      
      # Perform a post-hoc test if significant difference found:
      if (kruskal_p_zone < 0.05 | kruskal_p_type < 0.05 | kruskal_p_month < 0.05) {
        cat ( "kruskal wallis significant so doing post-hoc test:\n\n")
      }
      if (kruskal_p_zone < 0.05) {
        print ( kruskal_z <- kruskalmc(this_OTU$value ~ factor(this_OTU$zone)) )
      }
      if (kruskal_p_type < 0.05) {
        print ( kruskal_t <- kruskalmc(this_OTU$value ~ factor(this_OTU$type)) )
      }
      if (kruskal_p_month < 0.05) {
        print ( kruskal_t <- kruskalmc(this_OTU$value ~ factor(this_OTU$month)) )
      }
      # Keep a record of the p-values - these variables are stored in the OTU_prob_table which is the main output of this script.
      # In this case we are keeping the Kruskal Wallis results because for some of our data the residuals are not normally distributed
      stat_result_type <- kruskal_p_type
      stat_result_zone <- kruskal_p_zone
      stat_result_month <- kruskal_p_month
    }
    # No tests are done if we don't have both factors,
    # but the script can easily be changed to allow that
    
    OTU_prob_table[OTU_name,colnames] <- as.character(tax_table(my_phyloseq)[ taxa_names(my_phyloseq)==OTU_name])
    
    OTU_prob_table[OTU_name,"prob_type"] <- stat_result_type
    if(length(unique(sample_data(my_phyloseq)$zone)) >1 ) {
      OTU_prob_table[OTU_name,"prob_zone"] <- stat_result_zone
    } else {
      OTU_prob_table[OTU_name,"prob_zone"] <- NULL
    }
    OTU_prob_table[OTU_name,"prob_month"] <- stat_result_month
    
    OTUs_and_env_data.melt.OTU <- OTUs_and_env_data.melt[OTUs_and_env_data.melt$variable==OTU_name,]
    #lm.carbon <- lm( OTUs_and_env_data.melt.OTU$Carbon ~ OTUs_and_env_data.melt.OTU$value )
    #p.carbon <- summary(lm.carbon)$coefficients["OTUs_and_env_data.melt.OTU$value","Pr(>|t|)"]
    #m.carbon <- summary(lm.carbon)$coefficients["OTUs_and_env_data.melt.OTU$value","Estimate"]
    # In v.4 changed the method to Spearman:
    result <- cor.test(OTUs_and_env_data.melt.OTU$Carbon, OTUs_and_env_data.melt.OTU$value,method="spearman")
    p.carbon <- result$p.value
    m.carbon <- result$estimate
    
    OTU_prob_table[OTU_name,"carbon_rho"] <- m.carbon
    OTU_prob_table[OTU_name,"prob_carbon"] <- p.carbon
    
    #put some numbers in the OTU prob table for the purpose positioning markers on plots later:
    OTU_prob_table[OTU_name,"OTU_min"] <- min(this_OTU$value)
    OTU_prob_table[OTU_name,"OTU_max"] <- max(this_OTU$value)
    OTU_prob_table[OTU_name,"OTU_mid"] <- min(this_OTU$value) + (0.5*max(this_OTU$value) )
  }# end stats
  
  OTU_prob_table$prob_type_FDR <- p.adjust(OTU_prob_table$prob_type,method="fdr")
  OTU_prob_table$prob_zone_FDR <- p.adjust(OTU_prob_table$prob_zone,method="fdr")
  OTU_prob_table$prob_month_FDR <- p.adjust(OTU_prob_table$prob_month,method="fdr")
  OTU_prob_table$prob_carbon_FDR <- p.adjust(OTU_prob_table$prob_carbon,method="fdr")

  prob <- OTU_prob_table$prob_carbon_FDR < 0.05
  slope.copiotroph <- OTU_prob_table$carbon_rho > 0
  slope.oligotroph <- OTU_prob_table$carbon_rho < 0
  
  OTU_prob_table[prob&slope.copiotroph,"classification"] <- "copiotroph"
  OTU_prob_table[prob&slope.oligotroph,"classification"] <- "oligotroph"
  OTU_prob_table[OTU_prob_table$Phylum=="Cyanobacteria","classification"] <- "Cyanobacteria"
  OTU_prob_table[OTU_prob_table$Phylum=="Chloroflexi","classification"] <- "Chloroflexi"
  OTU_prob_table[is.na(OTU_prob_table$classification),"classification"] <- "N D"

  OTU_prob_table$c_cor <- ""
  OTU_prob_table[prob&slope.copiotroph,"c_cor"] <- "+"
  OTU_prob_table[prob&slope.oligotroph,"c_cor"] <- "-"
  
  OTU_prob_table$c_cor <- sig_code2(x_in=OTU_prob_table$prob_carbon_FDR,sym=OTU_prob_table$c_cor)
  
  OTU_prob_table$type_sig <- sig_code(OTU_prob_table$prob_type_FDR)
  OTU_prob_table$zone_sig <- sig_code(OTU_prob_table$prob_zone_FDR)
  OTU_prob_table$month_sig <- sig_code(OTU_prob_table$prob_month_FDR)
  
  candidate_copiotroph <- OTU_prob_table[OTU_prob_table$classification=="copiotroph",c("Phylum","Class", "Order",  "Family",  "Genus",	"Species","carbon_rho","prob_carbon_FDR")]
  candidate_oligotroph <- OTU_prob_table[OTU_prob_table$classification=="oligotroph",c("Phylum","Class", "Order",  "Family",	"Genus",	"Species","carbon_rho","prob_carbon_FDR")]
 
  # get ecological classification into the melted dataframe
  index <- match(OTUs_and_env_data.melt$variable, row.names(OTU_prob_table))
  OTUs_and_env_data.melt$ecology <- OTU_prob_table[index,"classification"]
  OTUs_and_env_data.melt$c_cor <- OTU_prob_table[index,"c_cor"]
  
  ################################
  
  if(topN < 20) { # suppress the plot if there are a lot of comparisons
    # TREE
    my_phyloseq.merged <- merge_samples(my_phyloseq, "zone_type")
    # my_phyloseq.merged <- merge_samples(my_phyloseq, "type") # alternative
    #sample_data(my_phyloseq.merged)$zone_type <- factor(sample_names(my_phyloseq.merged))
    #sample_data(my_phyloseq.merged)$zone <- gsub(pattern="_.*",replacement="",x=sample_data(my_phyloseq.merged)$zone_type)
    #sample_data(my_phyloseq.merged)$type <- gsub(pattern=".*_",replacement="",x=sample_data(my_phyloseq.merged)$zone_type)
    
    # NOTE - phyloseq doesn't offer the option of fill by factor, so I have just set all to white
    #fill_pal_bodge <- fill_pal
    #fill_pal_bodge <- rep("white",9)
    #names(fill_pal_bodge) <- names(fill_pal)
    #OTU_tree <- plot_tree(my_phyloseq.merged,color="zone",shape="type",size="abundance",label.tips = "taxa_names") # ,base.spacing = 0.04
    #OTU_tree = OTU_tree + scale_shape_manual(values=shape_scale,name="Depth") + scale_fill_manual(values=fill_pal_bodge)+ scale_colour_manual(values=fill_pal)
    #OTU_tree = OTU_tree + labs(title=plot_title) 
    
    # ZONE
    #fill2 <- as.character(fill_pal[1:2])
    zone_plot <- (ggplot(data=OTUs_and_env_data.melt,aes(x=zone_code,y=value, fill=depth2)) + geom_point(size=2,position=position_dodge(width=0.5)) 
                  + facet_wrap(~variable) 
                  + theme_bw() # base_size = 12, base_family = ""
                  + theme(aspect.ratio=1)
                  + labs(title=plot_title) 
                  + scale_x_discrete(name="vegetation") 
                  + scale_y_continuous(name="% abundance",trans = 'log10',breaks=c(0.1,1,10,100))  
                  + geom_boxplot()  
                  + theme(axis.text.x = element_text(angle = -90, vjust = 0.4))
                  + scale_fill_manual(values=fill_pal,name="Depth") )
                  #+ scale_shape_manual(values=zone_shape_scale,name="test"))
    #zone_plot
    
    # CARBON
    ypos <- (OTU_prob_table$OTU_max - 0.1*OTU_prob_table$OTU_max) 
    ann_text <- data.frame(Carbon = 2.2,value = ypos,lab = OTU_prob_table$c_cor,variable=row.names(OTU_prob_table))
    
    carbon_plot <- (ggplot(data=OTUs_and_env_data.melt,aes(x=Carbon,y=value)) + geom_point(size=3,aes(colour=depth2,shape=depth2)) 
                    + facet_wrap(~variable,scales="free_y") 
                    + theme_bw() # base_size = 12, base_family = "" 
                    + theme(aspect.ratio=1)
                    + labs(title=plot_title)   
                    + scale_y_continuous(name="% abundance") 
                    + scale_shape_manual(values=shape_scale, name="Depth")
                    + scale_colour_manual(values=color_pal, name="Depth")
                    + theme(axis.text.x = element_text(angle = -90, vjust = 0.4))
                    + geom_text(data = ann_text,aes(label = lab),size=12,alpha=0.6, hjust=1) )

    #png(filename=OTU_tree_file,width=800,height=800)  
    #print(OTU_tree)
    #dev.off()
    
    png(filename=carbon_plot_file,width=800,height=800)  
    print(carbon_plot)
    dev.off()
    
    png(filename=zone_plot_file,width=800,height=800)
    print(zone_plot)
    dev.off()
    
  }  else {
    cat(plot_title,":","Plots suppressed due to large number of comparisons")
    OTU_tree <- NULL
    zone_plot <- NULL
    carbon_plot <- NULL
  }
  OTU_tree <- NULL
  OTU_prob_table$OTU <- row.names(OTU_prob_table)
  write.csv(OTU_prob_table,writefile)

  write.csv(OTU_prob_table,writefile)
  
  if(verbose){ 
    cat(plot_title,":",nrow(candidate_oligotroph),"possible oligotrophs found\n")
    print(candidate_oligotroph)
    cat(plot_title,":",nrow(candidate_copiotroph),"possible copiotrophs found\n")
    print(candidate_copiotroph)
    cat("Note that copiotroph/oligotroph detection will not work effectively if zones are omitted from the phyloseq object")
  
  cat("\nOTU taxonomy:\n")
  print(OTU_prob_table[c("Phylum", "Class",  "Order", "Family", "Genus" )])
  cat("Full results saved to:",writefile,"\n")
  cat("stat tables were saved in folder",STATS_dir,"\n")
  #return(paste("<a href='",writefile,"'>",plot_title,"</a>\n(top ",topN," of ",n_input," taxa in group ",plot_title,")",sep=""))
  }
  #names(OTU_prob_table)
  
  t <- OTU_prob_table[c("OTU","Phylum", "Class",  "Order", "Family", "Genus" )]
  t.string <- ""
  for(n in 1:length(row.names(t))) {
   t.string = paste(t.string , (paste(t[n,],collapse=";")), sep="<br>\n")
  }

  plot_list <- list(
    "zone" = list(plot = zone_plot, file = zone_plot_file, name="zone plot"),
    "carbon" = list(plot = carbon_plot, file = carbon_plot_file, name="carbon plot"),
    "tree" = list(plot = OTU_tree , file = OTU_tree_file,name="tree") )
  
  links <- ""
  for(this_plot in plot_list) {
    links = paste(links, paste("<a href='",this_plot$file,"'>",plot_title," [top ",topN,"] ",this_plot$name,"</a>",sep=""),sep="<br>\n")
  }
#  for(this_plot in plot_list) {
#    links = paste(links, paste("[",plot_title," top ",topN,this_plot$name,"](",this_plot$file,"]",sep=""),sep=" <br>\n")
#  }
  
  results_link <- paste("<a href='",writefile,"'>",plot_title,": OTU taxonomy and statistical results</a><br>",sep="")
  stat_link <- paste("<a href='",STATS_dir,"'>Directory containing stat tables for each test</a>",sep="")
  # post_hoc_link <- paste("<a href='",post_hoc_file,"'>Post-hoc tests for significant results</a><br>",sep="")
  #post_hoc_link = post_hoc_link,
  
  result_list <- list(
    results_link = results_link,
    stat_link = stat_link,
    plot_links = links)

  return_list <- list(
     "results" = result_list,
     "plots" = plot_list,
     "taxonomy" = t.string,
     "normality" = shapiro.test,
     "results_table" =  OTU_prob_table,
     "OTU_table" = OTUs_and_env_data.extended
     )
  save(return_list,file=paste(objects_dir,"/",plot_title,"_OTU_analysis.rda",sep=""))
  
  sink()
  # close(outfile)
  return(return_list)
}

# Example of how to use and how to get the results out after:
if(FALSE){
test <- my_OTU_plot(expt,9,FALSE,"test")

test.plots <- test$plots
# print all the plot filenames
for (item in test.plots) {
 print(item$file) 
}

with(test, normality$p.value)
with(test, post_hoc)
tab <- with(test, results_table[c("OTU","Phylum", "Class",  "Order", "Family", "Genus" )])
tab <- with(test, results_table)
tab[c("OTU","Phylum", "Class",  "Order", "Family", "Genus" )]
tab$type_sig <- ""
tab[tab$prob_type_FDR < 0.05,"type_sig"] <- "*"

carbon_plot + annotate("text", label = "+", x = 2, y = 15, size = 8, colour = "red")
carbon_plot + geom_text(data = NULL, x = 2, label = "+",aes(y = max(value)))

max(4:9)
library(xtable)
print(xtable(aov), type="html")

# print one of the plots
test.plots[["zone"]]$plot

# print the results list
print(test$results)
#dev.off()
}
