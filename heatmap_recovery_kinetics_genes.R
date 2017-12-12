#Creating a heatmap for recovery kinetics genes (core upregulated, core downregulated, important non-canonical immune genes)
#May 8th, 2017
#Joo Hyun Im (ji72)

####################################################
#1. Plot a heatmap of genes
####################################################

rm(list=ls(all=TRUE))
library(gplots); library(RColorBrewer); library(cluster)

#1) Filter out the EE-EE-EE genes
#2) Split the genes based on situation (12h>0, 36h>0, 132h>0)
#3) Remove any genes that have 12hr>0 and 36hr<0 OR 12hr<0 and 36hr>0 (decided after talking to Nicolas and Brian)

get_Kinetics_with_directions = function(core_genes_address, direction){
    
    path_data_NS_removed = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.based.on.edgeR.results/recovery_genes/recovery_pattern_of_genes.txt", header=T)
    cat("The initial number of genes: ", dim(path_data_NS_removed)[1])
    expression_data = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_basic_comparison_FDR_converted_to_Y-N_at_least_one_sig.txt", header=T) #fc
    colnames(expression_data)[1:2] = c("gene_id", "gene_name")
    core_regulated = read.table(core_genes_address, header=T)
    
    path_data_NS_removed_core_dir_only = path_data_NS_removed[which(path_data_NS_removed$gene.id %in% core_regulated$gene_id),]
    expression_data = expression_data[which(expression_data$gene_id %in% core_regulated$gene_id),]
    expression_data_fc_only = expression_data[,c(1:2,seq(3,43,2), seq(53,64,2))] #for fc
    
    heatmap_number = matrix(NA, ncol = (dim(path_data_NS_removed_core_dir_only)[2]-2), nrow = dim(path_data_NS_removed_core_dir_only)[1])
    
    counting_cases_of_opposite_directions = 0
    for (i in 1:dim(path_data_NS_removed_core_dir_only)[1]){ #gene
        for (j in 1:(dim(path_data_NS_removed_core_dir_only)[2]-2)){ #condition
            pattern = as.character(path_data_NS_removed_core_dir_only[i,(j+2)])
            if (pattern == "EE-EE-EE"){ #mark the gene/condition that did not significantly change
                #heatmap_number[i,j] = NA
            }
            
            else{
                #cat("i: ", i, " j: ", j, " j*3: ", j*3, " (j*3)+1: ",(j*3)+1, "\n")
                twelve = expression_data_fc_only[i,j*3]
                thirtysix = expression_data_fc_only[i,(j*3)+1]
                onethirtytwo = expression_data_fc_only[i,(j*3)+2]
                
                ratio=0; onethirtytwo_added =0
                
                if (direction == "Up"){
                    if (twelve > 0 & thirtysix > 0 & onethirtytwo > 0){
                        if (onethirtytwo > twelve & onethirtytwo > thirtysix){
                            pickmax = max(twelve, thirtysix)
                            ratio = round( -abs(onethirtytwo - pickmax)/abs(pickmax), 3)
                            heatmap_number[i,j] = ratio
                        }
                        else{
                            pickmax = max(twelve, thirtysix)
                            ratio = round( abs(onethirtytwo-pickmax)/abs(pickmax), 3)
                            heatmap_number[i,j] = ratio
                        }
                    }
                    else if (twelve > 0 & thirtysix > 0 & onethirtytwo < 0){
                        pickmax = max(twelve, thirtysix)
                        onethirtytwo_added = pickmax + abs(onethirtytwo)
                        ratio = round(( onethirtytwo_added/pickmax ), 3)
                        heatmap_number[i,j] = ratio
                    }
                    else if (twelve < 0 & thirtysix > 0 & onethirtytwo > 0){
                        counting_cases_of_opposite_directions = counting_cases_of_opposite_directions + 1
                        #heatmap_number[i,j] = NA
                        #                     if ( abs(twelve) < abs(onethirtytwo) & abs(onethirtytwo) < abs(thirtysix) ){
                        #                         ratio = round(( abs(onethirtytwo)/abs(thirtysix) ), 3)
                        #                         heatmap_number[i,j] = ratio
                        #                     }
                        #                     else if ( abs(onethirtytwo) < abs(twelve) & abs(twelve) < abs(thirtysix) ){
                        #                         ratio = round(( abs(onethirtytwo)/abs(thirtysix) ), 3)
                        #                         heatmap_number[i,j] = ratio
                        #                     }
                        #                     else {
                        #                         onethirtytwo_added = abs(twelve) + abs(onethirtytwo)
                        #                         ratio = round(( abs(onethirtytwo_added)/abs(twelve) ), 3)
                        #                             heatmap_number[i,j] = ratio
                        #                     }
                    } #NA as of 5/1/17
                    else if (twelve < 0 & thirtysix < 0 & onethirtytwo > 0){
                        pickmax = max( abs(twelve), abs(thirtysix) )
                        onethirtytwo_added = abs(pickmax) + onethirtytwo
                        ratio = round(( onethirtytwo_added/pickmax ), 3)
                        heatmap_number[i,j] = ratio
                    }
                    else if (twelve > 0 & thirtysix < 0 & onethirtytwo > 0){
                        counting_cases_of_opposite_directions = counting_cases_of_opposite_directions + 1
                        #heatmap_number[i,j] = NA
                        #                     if ( abs(twelve) > abs(thirtysix)){
                        #                         ratio = round(( onethirtytwo/twelve ), 3)
                        #                         heatmap_number[i,j] = ratio
                        #                     }
                        #                     else{
                        #                         onethirtytwo_added = onethirtytwo + abs(thirtysix)
                        #                         ratio = round(( onethirtytwo_added/ abs(thirtysix) ), 3)
                        #                         heatmap_number[i,j] = ratio
                        #                      }
                    } #NA as of 5/1/17
                    else if (twelve > 0 & thirtysix < 0 & onethirtytwo < 0){
                        counting_cases_of_opposite_directions = counting_cases_of_opposite_directions + 1
                        #heatmap_number[i,j] = NA
                        #                     if ( abs(twelve) < abs(onethirtytwo) & abs(onethirtytwo) < abs(thirtysix) ){
                        #                         ratio = round(( abs(onethirtytwo)/abs(thirtysix) ), 3)
                        #                         heatmap_number[i,j] = ratio
                        #                     }
                        #                     else if ( abs(onethirtytwo) < abs(twelve) & abs(twelve) < abs(thirtysix) ){
                        #                         ratio = round(( abs(onethirtytwo)/abs(thirtysix) ), 3)
                        #                         heatmap_number[i,j] = ratio
                        #                     }
                        #                     else {
                        #                         onethirtytwo_added = abs(twelve) + abs(onethirtytwo)
                        #                         ratio = round(( abs(onethirtytwo_added)/abs(twelve) ), 3)
                        #                         heatmap_number[i,j] = ratio
                        #                     }
                    } #NA as of 5/1/17
                    else if (twelve < 0 & thirtysix > 0 & onethirtytwo < 0){
                        counting_cases_of_opposite_directions = counting_cases_of_opposite_directions + 1
                        #heatmap_number[i,j] = NA
                        #                     if ( abs(twelve) > abs(thirtysix)){
                        #                         ratio = round(( abs(onethirtytwo)/ abs(twelve) ), 3)
                        #                         heatmap_number[i,j] = ratio
                        #                     }
                        #                     else{
                        #                         onethirtytwo_added = abs(onethirtytwo) + thirtysix
                        #                         ratio = round(( onethirtytwo_added / thirtysix ), 3)
                        #                         heatmap_number[i,j] = ratio
                        #                     }
                    } #NA as of 5/1/17
                    else if (twelve < 0 & thirtysix < 0 & onethirtytwo < 0){
                        cat("exclude this: ", as.character(path_data_NS_removed_core_dir_only[i,2]), "\n")
                        #                         if ( abs(onethirtytwo) > abs(twelve) & abs(onethirtytwo) > abs(thirtysix)){
                        #                             pickmax = max( abs(twelve), abs(thirtysix) )
                        #                             ratio = round( -abs(onethirtytwo - (-pickmax))/abs(pickmax), 3)
                        #                             heatmap_number[i,j] = ratio
                        #                         }
                        #                         else{
                        #                             pickmax = max( abs(twelve), abs(thirtysix) )
                        #                             ratio = round(( abs(onethirtytwo-(-pickmax)) / pickmax ), 3)
                        #                             heatmap_number[i,j] = ratio
                        #                         }
                    } #remove genes that are all negatively induced
                }
                else{
                    if (twelve > 0 & thirtysix > 0 & onethirtytwo > 0){
                        cat("exclude this: ", as.character(path_data_NS_removed_core_dir_only[i,2]), "\n")
                        #                     if (onethirtytwo > twelve & onethirtytwo > thirtysix){
                        #                         pickmax = max(twelve, thirtysix)
                        #                         ratio = round( -abs(onethirtytwo - pickmax)/abs(pickmax), 3)
                        #                         heatmap_number[i,j] = ratio
                        #                     }
                        #                     else{
                        #                         pickmax = max(twelve, thirtysix)
                        #                         ratio = round( abs(onethirtytwo-pickmax)/abs(pickmax), 3)
                        #                         heatmap_number[i,j] = ratio
                        #                     }
                    }#remove genes that are all positively induced
                    else if (twelve > 0 & thirtysix > 0 & onethirtytwo < 0){
                        pickmax = max(twelve, thirtysix)
                        onethirtytwo_added = pickmax + abs(onethirtytwo)
                        ratio = round(( onethirtytwo_added/pickmax ), 3)
                        heatmap_number[i,j] = ratio
                    }
                    else if (twelve < 0 & thirtysix > 0 & onethirtytwo > 0){
                        counting_cases_of_opposite_directions = counting_cases_of_opposite_directions + 1
                        #heatmap_number[i,j] = NA
                        #                     if ( abs(twelve) < abs(onethirtytwo) & abs(onethirtytwo) < abs(thirtysix) ){
                        #                         ratio = round(( abs(onethirtytwo)/abs(thirtysix) ), 3)
                        #                         heatmap_number[i,j] = ratio
                        #                     }
                        #                     else if ( abs(onethirtytwo) < abs(twelve) & abs(twelve) < abs(thirtysix) ){
                        #                         ratio = round(( abs(onethirtytwo)/abs(thirtysix) ), 3)
                        #                         heatmap_number[i,j] = ratio
                        #                     }
                        #                     else {
                        #                         onethirtytwo_added = abs(twelve) + abs(onethirtytwo)
                        #                         ratio = round(( abs(onethirtytwo_added)/abs(twelve) ), 3)
                        #                             heatmap_number[i,j] = ratio
                        #                     }
                    } #NA as of 5/1/17
                    else if (twelve < 0 & thirtysix < 0 & onethirtytwo > 0){
                        pickmax = max( abs(twelve), abs(thirtysix) )
                        onethirtytwo_added = abs(pickmax) + onethirtytwo
                        ratio = round(( onethirtytwo_added/pickmax ), 3)
                        heatmap_number[i,j] = ratio
                    }
                    else if (twelve > 0 & thirtysix < 0 & onethirtytwo > 0){
                        counting_cases_of_opposite_directions = counting_cases_of_opposite_directions + 1
                        #heatmap_number[i,j] = NA
                        #                     if ( abs(twelve) > abs(thirtysix)){
                        #                         ratio = round(( onethirtytwo/twelve ), 3)
                        #                         heatmap_number[i,j] = ratio
                        #                     }
                        #                     else{
                        #                         onethirtytwo_added = onethirtytwo + abs(thirtysix)
                        #                         ratio = round(( onethirtytwo_added/ abs(thirtysix) ), 3)
                        #                         heatmap_number[i,j] = ratio
                        #                      }
                    } #NA as of 5/1/17
                    else if (twelve > 0 & thirtysix < 0 & onethirtytwo < 0){
                        counting_cases_of_opposite_directions = counting_cases_of_opposite_directions + 1
                        #heatmap_number[i,j] = NA
                        #                     if ( abs(twelve) < abs(onethirtytwo) & abs(onethirtytwo) < abs(thirtysix) ){
                        #                         ratio = round(( abs(onethirtytwo)/abs(thirtysix) ), 3)
                        #                         heatmap_number[i,j] = ratio
                        #                     }
                        #                     else if ( abs(onethirtytwo) < abs(twelve) & abs(twelve) < abs(thirtysix) ){
                        #                         ratio = round(( abs(onethirtytwo)/abs(thirtysix) ), 3)
                        #                         heatmap_number[i,j] = ratio
                        #                     }
                        #                     else {
                        #                         onethirtytwo_added = abs(twelve) + abs(onethirtytwo)
                        #                         ratio = round(( abs(onethirtytwo_added)/abs(twelve) ), 3)
                        #                         heatmap_number[i,j] = ratio
                        #                     }
                    } #NA as of 5/1/17
                    else if (twelve < 0 & thirtysix > 0 & onethirtytwo < 0){
                        counting_cases_of_opposite_directions = counting_cases_of_opposite_directions + 1
                        #heatmap_number[i,j] = NA
                        #                     if ( abs(twelve) > abs(thirtysix)){
                        #                         ratio = round(( abs(onethirtytwo)/ abs(twelve) ), 3)
                        #                         heatmap_number[i,j] = ratio
                        #                     }
                        #                     else{
                        #                         onethirtytwo_added = abs(onethirtytwo) + thirtysix
                        #                         ratio = round(( onethirtytwo_added / thirtysix ), 3)
                        #                         heatmap_number[i,j] = ratio
                        #                     }
                    } #NA as of 5/1/17
                    else if (twelve < 0 & thirtysix < 0 & onethirtytwo < 0){
                        if ( abs(onethirtytwo) > abs(twelve) & abs(onethirtytwo) > abs(thirtysix)){
                            pickmax = max( abs(twelve), abs(thirtysix) )
                            ratio = round( -abs(onethirtytwo - (-pickmax))/abs(pickmax), 3)
                            heatmap_number[i,j] = ratio
                        }
                        else{
                            pickmax = max( abs(twelve), abs(thirtysix) )
                            ratio = round(( abs(onethirtytwo-(-pickmax)) / pickmax ), 3)
                            heatmap_number[i,j] = ratio
                        }
                    } 
                    
                    
                }
                
            }
        }
    }
    cat("counting_cases_of_opposite_directions: ", counting_cases_of_opposite_directions) #responding to reviewers (12/1/2017)
    
    heatmap_number_with_name = cbind(path_data_NS_removed_core_dir_only[,c(1:2)], heatmap_number)
    colnames(heatmap_number_with_name) = c("gene_id", "gene_name", "clean.prick","M.luteus","E.coli","S.marcescens","E.faecalis","P.rettgeri","Ecc15","E.fae.heatkilled", "P.rett.heatkilled")
    
    heatmap_number_with_name = heatmap_number_with_name[match(core_regulated$gene_id, heatmap_number_with_name$gene_id),]
    return(heatmap_number_with_name)
    
}


#Core upregulated genes:
#core_upregulated_kinetics = get_Kinetics_with_directions("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_(7_or_more_conditions)/core_genes_7_plus_live_bacteria_upregulated_gene_list.txt", "Up") #all
core_upregulated_kinetics = get_Kinetics_with_directions("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.based.on.edgeR.results/recovery_genes/quantification_of_recovery/spectrum/core_upgenes_selected.txt", "Up") #selected

rownames(core_upregulated_kinetics) = core_upregulated_kinetics$gene_name
core_upregulated_kinetics_value_only = as.matrix(core_upregulated_kinetics[,3:11])
#core_upregulated_kinetics_value_only = core_upregulated_kinetics_value_only[-which(rownames(core_upregulated_kinetics_value_only) == c("CG32284")),]  #remove the case that has NA in every single condition (165 genes), I need to remove this to get the clustering work. This is not needed if the input is selected core genes
core_upregulated_kinetics_value_only = core_upregulated_kinetics_value_only[,c(1,9,8,7,2,3,4,6,5)] #change the order to: sw, heatkilleds, ecc15, rest


#my_palette = colorRampPalette(c("#762a83", "#9970ab", "#c2a5cf", "#e7d4e8", "gray93", "#d9f0d3", "#a6dba0", "#5aae61", "#1b7837"))(n = 299) #tri-color, green, purle, light gray
#my_palette = colorRampPalette(c("#f7fcb9","#d9f0a3","#addd8e","#78c679","#41ab5d","#238443","#006837","#004529"))(n = 299) #uni-directional sequential green
#my_palette = colorRampPalette(c("#edf8b1","#c7e9b4","#7fcdbb","#41b6c4","#1d91c0","#225ea8","#253494","#081d58"))(n = 299) #uni-directional sequential blue
my_palette = colorRampPalette(c("#762a83", "#9970ab", "#c2a5cf", "#e7d4e8", "gray87", "#d9f0d3", "#a6dba0", "#5aae61", "#1b7837"))(n = 299) #tri-color, green, purple, gray

#Get the heatmap
heatmap.2(core_upregulated_kinetics_value_only, Rowv=TRUE, Colv=FALSE, dendrogram="none", symm=F, scale="none", na.rm=TRUE, col=my_palette, na.color="white", trace="none", key=F, keysize=1, density.info="none", margin=c(6,16), cexRow=1) #works
heatmap.2(core_upregulated_kinetics_value_only, Rowv=TRUE, Colv=FALSE, dendrogram="row", symm=F, scale="none", na.rm=TRUE, col=my_palette, na.color="white", trace="none", key=F, keysize=1, density.info="none", margin=c(6,16), cexRow=1) #works, add a dendrogram

#Get the color key (square)
# heatmap.2(core_upregulated_kinetics_value_only, Rowv=TRUE, Colv=FALSE, dendrogram="row", symm=F, scale="none", na.rm=TRUE, col=my_palette, na.color="white", trace="none", key=T, keysize=5, density.info="none", key.xlab="Magnitude of expression changes", srtCol=45, lmat = rbind(c(4,0),c(2,1),c(3,0)), lwid = c(1.5,4), lhei = c(1.5,4,1))
#Get the color key (rectangle)
heatmap.2(core_upregulated_kinetics_value_only, Rowv=TRUE, Colv=FALSE, dendrogram="row", symm=F, scale="none", na.rm=TRUE, col=my_palette, na.color="white", trace="none", key=T, keysize=1, density.info="none", key.xlab="Spectrum of expression patterns", srtCol=45, lmat = rbind(c(0,4),c(2,1),c(3,0)), lwid = c(1.5,4), lhei = c(1,4,1)) 


#plotting the rows based on clustering of the following samples: Ml, Ec, Sm, Ef, Pr
# core_upregulated_kinetics_value_only = as.matrix(core_upregulated_kinetics[,c(4:8)])
# test = heatmap.2(core_upregulated_kinetics_value_only, Rowv=TRUE, Colv=FALSE, dendrogram="none", symm=F, scale="column", na.rm=TRUE, col=my_palette, na.color="white", trace="none", key=F, density.info="none", margin=c(6,16), cexRow=0.8) #works
# clustered_rows = core_upregulated_kinetics_value_only[rev(test$rowInd),]
# core_upregulated_kinetics_value_only = core_upregulated_kinetics_value_only[match(rownames(clustered_rows), rownames(core_upregulated_kinetics_value_only)),]
# 
# #Get the heatmap
# heatmap.2(core_upregulated_kinetics_value_only, Rowv=TRUE, Colv=FALSE, dendrogram="none", symm=F, scale="column", na.rm=TRUE, col=my_palette, na.color="white", trace="none", key=F, density.info="none", margin=c(6,16), cexRow=0.9) #works

#Change colors








###########################################################################
#Core downregulated genes:
#core_downregulated_kinetics = get_Kinetics_with_directions("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_(7_or_more_conditions)/core_genes_7_plus_live_bacteria_downregulated_gene_list.txt", "Down") #all
core_downregulated_kinetics = get_Kinetics_with_directions("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.based.on.edgeR.results/recovery_genes/quantification_of_recovery/spectrum/core_downgenes_selected.txt", "Down") #selected

rownames(core_downregulated_kinetics) = core_downregulated_kinetics$gene_name
core_downregulated_kinetics_value_only = as.matrix(core_downregulated_kinetics[,3:11])
core_downregulated_kinetics_value_only = core_downregulated_kinetics_value_only[,c(1,9,8,7,2,3,4,6,5)] #change the condition order
# core_downregulated_kinetics_value_only = core_downregulated_kinetics_value_only[c("CG11200","PKD","CG9034","CR43051", "CG16926", "CG14629", "CG12091", "CG12338", "CG17124", "CG9664", "Non1", "Tps1", "Lectin-galC1", "alpha-Est2", "Try29F", "Cyp6a8", "Jheh1", "ninaD","CG34166","CG3739","CG3940","CG9510", "CG3609","Jon99Ci", "tobi", "CG42369", "CG34331",  "GstE8", "CG18003", "CG31104", "CG31075", "CG3699", "Cyp6a18", "fit", "CG3301", "CG16986", "CG16985", "CG7997", "Mur18B", "CG15263", "glob1", "Cyp6g1", "Actbeta", "alpha-Est7", "Cyp4e2", "Spat", "Lsp1beta", "Lsp2", "fat-spondin", "Ugt37b1", "CG13365", "CG3603", "CG2233", "CG15369", "regucalcin", "CG9512", "dob", "Iris",  "dsb", "Cyp4d20", "CG4950", "CG18135", "CG5618", "Est-Q", "CG16749",  "CG10824", "CG10208", "to", "Obp99b", "Jon99Fii", "lectin-28C", "Nplp2", "CG31548", "CG33301", "26-29-p","CG42788", "Acsl","Ugt36Bc","Npc2h", "lectin-33A","Cyp4ac2", "CG17108", "Mal-B1", "CG2070", "CG16898","Cyp6t1"),] #manually change the gene order to group blues, reds, and the mix


#my_palette = colorRampPalette(c("#762a83", "#9970ab", "#c2a5cf", "#e7d4e8", "gray93", "#d9f0d3", "#a6dba0", "#5aae61", "#1b7837"))(n = 299) #tri-color
#my_palette = colorRampPalette(c("#f7fcb9","#d9f0a3","#addd8e","#78c679","#41ab5d","#238443","#006837","#004529"))(n = 299) #uni-directional sequential green
#my_palette = colorRampPalette(c("#edf8b1","#c7e9b4","#7fcdbb","#41b6c4","#1d91c0","#225ea8","#253494","#081d58"))(n = 299) #uni-directional sequential blue
my_palette = colorRampPalette(c("#762a83", "#9970ab", "#c2a5cf", "#e7d4e8", "gray87", "#d9f0d3", "#a6dba0", "#5aae61", "#1b7837"))(n = 299) #tri-color, green, purple, gray



#Get the heatmap
heatmap.2(core_downregulated_kinetics_value_only, Rowv=TRUE, Colv=FALSE, dendrogram="none", symm=F, scale="none", na.rm=TRUE, col=my_palette, na.color="white", trace="none", key=F, density.info="none", margin=c(6,16), cexRow=0.9) #works
heatmap.2(core_downregulated_kinetics_value_only, Rowv=TRUE, Colv=FALSE, dendrogram="row", symm=F, scale="none", na.rm=TRUE, col=my_palette, na.color="white", trace="none", key=F, keysize=1, density.info="none", margin=c(6,16), cexRow=0.9) #works, add a dendrogram


#Get the color key (square)
# heatmap.2(core_downregulated_kinetics_value_only, Rowv=FALSE, Colv=FALSE, dendrogram="none", symm=F, scale="none", na.rm=TRUE, col=my_palette, na.color="white", trace="none", key=T, keysize=5, density.info="none", key.xlab="Magnitude of expression changes", srtCol=45, lmat = rbind(c(4,0),c(2,1),c(3,0)), lwid = c(1.5,4), lhei = c(1.5,4,1))
# #Get the color key (rectangle)
heatmap.2(core_downregulated_kinetics_value_only, Rowv=TRUE, Colv=FALSE, dendrogram="row", symm=F, scale="none", na.rm=TRUE, col=my_palette, na.color="white", trace="none", key=T, keysize=1, density.info="none", key.xlab="Spectrum of expression patterns", srtCol=45, lmat = rbind(c(0,4),c(2,1),c(3,0)), lwid = c(1.5,4), lhei = c(1,4,1)) 


#Extra genes that Nicolas picked from PCs that are not core-genes:
other_genes_kinetics = get_Kinetics_with_directions("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.based.on.edgeR.results/recovery_genes/quantification_of_recovery/spectrum/list_of_non_imm_genes_from_PCs_upreg_row_sized_051217.txt", "Up") #the first 6 genes
#other_genes_kinetics = get_Kinetics_with_directions("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.based.on.edgeR.results/recovery_genes/quantification_of_recovery/spectrum/list_of_non_imm_genes_from_PCs_downreg_051217.txt", "Down") #13 genes

rownames(other_genes_kinetics) = other_genes_kinetics$gene_name
#rownames(other_genes_kinetics) = c("spirit", "Sp7", "cact", "SPE", "MP1", "Gs1", "sample2","sample3","sample4","sample5","sample6","sample7") #for upreg only
other_genes_kinetics_value_only = as.matrix(other_genes_kinetics[,3:11])
other_genes_kinetics_value_only = other_genes_kinetics_value_only[,c(1,9,8,7,2,3,4,6,5)] #change the condition order
my_palette = colorRampPalette(c("#762a83", "#9970ab", "#c2a5cf", "#e7d4e8", "gray93", "#d9f0d3", "#a6dba0", "#5aae61", "#1b7837"))(n = 299) #tri-color

#Get the heatmap with the color key -- for upregulated ones
heatmap.2(other_genes_kinetics_value_only, Rowv=FALSE, Colv=FALSE, dendrogram="none", symm=F, scale="none", na.rm=TRUE, col=my_palette, na.color="white", trace="none", key=T, keysize=1, density.info="none", key.xlab="Spectrum of expression patterns", margin=c(26,13), cexRow=1.15, lmat = rbind(c(0,4),c(2,1),c(3,0)), lwid = c(1.5,4), lhei = c(1,4,1)) 

#legend
heatmap.2(other_genes_kinetics_value_only, Rowv=FALSE, Colv=FALSE, dendrogram="none", symm=F, scale="none", na.rm=TRUE, col=my_palette, na.color="white", trace="none", key=T, keysize=1, density.info="none", key.xlab="Spectrum of expression patterns", lmat = rbind(c(0,4),c(2,1),c(3,0)), lwid = c(1.5,4), lhei = c(1,4,1)) 
1
#otherwise (downregulated)
my_palette_extra_down_only = colorRampPalette(c("gray93", "#d9f0d3", "#a6dba0", "#5aae61", "#1b7837"))(n = 299) #tri-color
heatmap.2(other_genes_kinetics_value_only, Rowv=FALSE, Colv=FALSE, dendrogram="none", symm=F, scale="none", na.rm=TRUE, col=my_palette_extra_down_only, na.color="white", trace="none", key=T, keysize=1, density.info="none", key.xlab="Spectrum of expression patterns", margin=c(26,13), cexRow=1.15, lmat = rbind(c(0,4),c(2,1),c(3,0)), lwid = c(1.5,4), lhei = c(1,4,1)) 

#legend - All of them are downregulated at first and then upregulated at the end
heatmap.2(other_genes_kinetics_value_only, Rowv=FALSE, Colv=FALSE, dendrogram="none", symm=F, scale="none", na.rm=TRUE, col=my_palette_extra_down_only, na.color="white", trace="none", key=T, keysize=1, density.info="none", key.xlab="Spectrum of expression patterns", lmat = rbind(c(0,4),c(2,1),c(3,0)), lwid = c(1.5,4), lhei = c(1,4,1)) 


##Get the heatmap with the color key -- for upregulated ones
#heatmap.2(other_genes_kinetics_value_only, Rowv=FALSE, Colv=FALSE, dendrogram="none", symm=F, scale="none", na.rm=TRUE, col=my_palette, na.color="white", trace="none", key=T, keysize=1, density.info="none", key.xlab="Spectrum of expression patterns", margin=c(30,13), cexRow=1.15, lmat = rbind(c(0,4),c(2,1),c(3,0)), lwid = c(1.5,4), lhei = c(1,4,1)) 


