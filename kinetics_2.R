#Kinetics II (4/27/2017)
rm(list=ls(all=TRUE))

#Only look at the genes that changed in at least one time point in a given condition (excluding EE-EE-EE).
#Only look at the genes that are consistently upregulated (in the case of core upregulated) or downregulated (in the case of core downregulated). For instance, a core upregulated genes in 7 conditions that is downregulated in another condition has been excluded.

get_Kinetics = function(core_genes_address, direction){
    
    path_data_NS_removed = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.based.on.edgeR.results/recovery_genes/recovery_pattern_of_genes.txt", header=T)
    expression_data = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_basic_comparison_FDR_converted_to_Y-N_at_least_one_sig.txt", header=T) #fc
    #expression_data = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_normalized_counts_from_all_genes_averaged_with_all_UCs.txt", header=T) #count
    colnames(expression_data)[1:2] = c("gene_id", "gene_name")
    core_regulated = read.table(core_genes_address, header=T)
    
    path_data_NS_removed_core_dir_only = path_data_NS_removed[which(path_data_NS_removed$gene.id %in% core_regulated$gene_id),]
    
    expression_data = expression_data[which(expression_data$gene_id %in% core_regulated$gene_id),]
    expression_data_fc_only = expression_data[,c(1:2,seq(3,43,2), seq(53,64,2))] #for fc
    #expression_data_count_only = expression_data[,c(1:2,4:24,29:34)] #for count
    
    heatmap_number = matrix(NA, ncol = 9, nrow = dim(path_data_NS_removed_core_dir_only)[1])
    late_start = c()
    for (i in 1:dim(path_data_NS_removed_core_dir_only)[1]){ #gene
        for (j in 1:(dim(path_data_NS_removed_core_dir_only)[2]-2)){ #condition
            pattern = as.character(path_data_NS_removed_core_dir_only[i,(j+2)])
            late_start_row = c()
            #print(pattern)
            
            twelve = expression_data_fc_only[i,j*3]
            thirtysix = expression_data_fc_only[i,(j*3)+1]
            onethirtytwo = expression_data_fc_only[i,(j*3)+2]
            
            if (direction == "Up") {
                if (twelve > 0 & thirtysix > 0 & onethirtytwo >0){
                    #if 132hr is maximum: either continuously (Up-Up-Up) or exponentially (EE-EE/Up-Up) increasing or Up-EE-Up
                    if ( max(abs(twelve), abs(thirtysix), abs(onethirtytwo)) == abs(onethirtytwo)) { 
                        heatmap_number[i,j] = NaN #"late"
                        #late_start_row = matrix(c(as.character(path_data_NS_removed_core_dir_only[i,2]),as.character(colnames(path_data_NS_removed_core_dir_only)[j+2])), 1,2)
                        #late_start = rbind(late_start, late_start_row)
                    }
                    #otherwise
                    else{
                        pickmax = max(twelve, thirtysix)
                        ratio = round((onethirtytwo/pickmax)*100, 2)
                        heatmap_number[i,j] = ratio
                    }
                }
                else{
                    heatmap_number[i,j] = NaN #"up-down"
                }
            }
            
            else if (direction == "Down") {
                if (twelve < 0 & thirtysix < 0 & onethirtytwo < 0){
                    #if 132hr is maximum: either continuously (Down-Down-Down) or exponentially (EE-EE/Down-Down) decreasing or Down-EE-Down
                    if ( max(abs(twelve), abs(thirtysix), abs(onethirtytwo)) == abs(onethirtytwo)) { 
                        heatmap_number[i,j] = NaN #"late"
                        #late_start_row = matrix(c(as.character(path_data_NS_removed_core_dir_only[i,2]),as.character(colnames(path_data_NS_removed_core_dir_only)[j+2])), 1,2)
                        #late_start = rbind(late_start, late_start_row)
                    }
                    #otherwise
                    else{ 
                        pickmax = max(abs(twelve), abs(thirtysix)) #already change the minus to plus if
                        ratio = round(( abs(onethirtytwo) /pickmax)*100, 2)
                        heatmap_number[i,j] = ratio
                    }
                }
                else{
                    heatmap_number[i,j] = NaN #"up-down"
                }

                
                

            }
        }
    }
    
    #cat("The ones excluded because the expression level at 132h was at its highest: ", late_start)
    
    heatmap_number_with_name = cbind(path_data_NS_removed_core_dir_only[,c(1:2)], heatmap_number)
    colnames(heatmap_number_with_name) = c("gene_id", "gene_name", "clean.prick","M.luteus","E.coli","S.marcescens","E.faecalis","P.rettgeri","Ecc15","E.fae.heatkilled", "P.rett.heatkilled")
    return(heatmap_number_with_name)
}

core_upregulated_kinetics = get_Kinetics("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_(7_or_more_conditions)/core_genes_7_plus_live_bacteria_upregulated_gene_list.txt", "Up")
core_downregulated_kinetics = get_Kinetics("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_(7_or_more_conditions)/core_genes_7_plus_live_bacteria_downregulated_gene_list.txt", "Down")

#core_downregulated_kinetics_up_down_removed = core_downregulated_kinetics
#core_downregulated_kinetics_up_down_removed$up_down_count = rowSums(core_downregulated_kinetics_up_down_removed == "up-down")

library(gplots); library(RColorBrewer)
my_palette = colorRampPalette(c("#cad3f9","#6b83fb","#041890"))(n = 299) 
heatmap.2(as.matrix(core_upregulated_kinetics[,3:11]), Rowv=FALSE, density.info="none", dendrogram="none",trace="none", symm=F, scale="none", col=my_palette, srtCol=45, key=T, keysize = 1, cexCol=1.5, cexRow=1,key.xlab="Expression")

my_palette2 = colorRampPalette(c("pink","red","purple"))(n = 299) 
heatmap.2(as.matrix(core_downregulated_kinetics[,3:11]), Rowv=FALSE, density.info="none", dendrogram="none",trace="none", symm=F, scale="none", col=my_palette2, srtCol=45, key=T, keysize = 1, cexCol=1.5, cexRow=1,key.xlab="Expression")



#Trying something different (4/27/2017 at 4pm and 5/1/2017 at 10am, 5/3/2017) 
rm(list=ls(all=TRUE))
library(gplots); library(RColorBrewer); library(cluster)

#1) Filter out the EE-EE-EE genes
#2) Split the genes based on situation (12h>0, 36h>0, 132h>0)
#3) Remove any genes that have 12hr>0 and 36hr<0 OR 12hr<0 and 36hr>0 (decided after talking to Nicolas and Brian)

get_Kinetics_with_directions = function(core_genes_address){
    
    path_data_NS_removed = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.based.on.edgeR.results/recovery_genes/recovery_pattern_of_genes.txt", header=T)
    expression_data = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_basic_comparison_FDR_converted_to_Y-N_at_least_one_sig.txt", header=T) #fc
    colnames(expression_data)[1:2] = c("gene_id", "gene_name")
    core_regulated = read.table(core_genes_address, header=T)
    
    path_data_NS_removed_core_dir_only = path_data_NS_removed[which(path_data_NS_removed$gene.id %in% core_regulated$gene_id),]
    expression_data = expression_data[which(expression_data$gene_id %in% core_regulated$gene_id),]
    expression_data_fc_only = expression_data[,c(1:2,seq(3,43,2), seq(53,64,2))] #for fc
    
    #testset
    #expression_data_fc_only = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.based.on.edgeR.results/recovery_genes/kinetics_example.txt", header=T)
    
    heatmap_number = matrix(NA, ncol = (dim(path_data_NS_removed_core_dir_only)[2]-2), nrow = dim(path_data_NS_removed_core_dir_only)[1])
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
                
                #cat("twelve: ", twelve, " thirtysix: ", thirtysix, " onethirtytwo: ", onethirtytwo, "\n")
                
                ratio=0; onethirtytwo_added =0
                if (twelve > 0 & thirtysix > 0 & onethirtytwo > 0){
                    pickmax = max(twelve, thirtysix)
                    ratio = round((onethirtytwo/pickmax), 3)
                    heatmap_number[i,j] = ratio
                }
                else if (twelve > 0 & thirtysix > 0 & onethirtytwo < 0){
                    pickmax = max(twelve, thirtysix)
                    onethirtytwo_added = pickmax + abs(onethirtytwo)
                    ratio = round(( onethirtytwo_added/pickmax ), 3)
                    heatmap_number[i,j] = ratio
                }
                else if (twelve < 0 & thirtysix > 0 & onethirtytwo > 0){
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
                    pickmax = max( abs(twelve), abs(thirtysix) )
                    ratio = round(( abs(onethirtytwo) / pickmax ), 3)
                    heatmap_number[i,j] = ratio
                }
                
                }
            }
    }
    #a = cbind(heatmap_number, expression_data_fc_only$Answer_Key_100)
    #cat("The ones excluded because the expression level at 132h was at its highest: ", late_start)
    
    heatmap_number_with_name = cbind(path_data_NS_removed_core_dir_only[,c(1:2)], heatmap_number)
    colnames(heatmap_number_with_name) = c("gene_id", "gene_name", "clean.prick","M.luteus","E.coli","S.marcescens","E.faecalis","P.rettgeri","Ecc15","E.fae.heatkilled", "P.rett.heatkilled")
    return(heatmap_number_with_name)
}
core_upregulated_kinetics = get_Kinetics_with_directions("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_(7_or_more_conditions)/core_genes_7_plus_live_bacteria_upregulated_gene_list.txt")
rownames(core_upregulated_kinetics) = core_upregulated_kinetics$gene_name

# my_palette = colorRampPalette(c("lightgrey","#80cdc1", "#018571"))(n = 299) 
# heatmap.2(as.matrix(core_upregulated_kinetics[,3:11]), density.info="none", dendrogram="row",trace="none", symm=F, scale="none", col=my_palette, srtCol=45, key=T, keysize = 0.5, cexCol=1, cexRow=0.8, key.xlab="Magnitude of expression changes", lmat = rbind(c(4,0),c(2,1),c(3,0)), lwid = c(1.5,4), lhei = c(1.5,4,1) ) 

#if scale by column (0 in the scale is raw number 1 = not recovered/similart to highest peak)
core_upregulated_kinetics_value_only = as.matrix(core_upregulated_kinetics[,3:11])
core_upregulated_kinetics_value_only = core_upregulated_kinetics_value_only[-which(rownames(core_upregulated_kinetics_value_only) == c("CG32284")),] #remove the case that has NA in every single condition (165 genes)
core_upregulated_kinetics_value_only = core_upregulated_kinetics_value_only[,c(1,9,8,2,3,4,7,6,5)] #change the order

my_palette = colorRampPalette(c("#2166ac", "lightgrey", "#b2182b"))(n = 299) #colorblind-safe blue, grey, red
heatmap.2(core_upregulated_kinetics_value_only, Rowv=TRUE, Colv=FALSE, dendrogram="row", symm=F, scale="column", na.rm=TRUE, col=my_palette, na.color="white", trace="none", key=F, density.info="none", key.xlab="Magnitude of expression changes", srtCol=45, margin=c(6,6)) 
heatmap.2(core_upregulated_kinetics_value_only, Rowv=TRUE, Colv=FALSE, dendrogram="none", symm=F, scale="column", na.rm=TRUE, col=my_palette, na.color="white", trace="none", key=F, density.info="none", key.xlab="Magnitude of expression changes", srtCol=45, margin=c(6,13), cexRow=1.3) #works

#just to get the key
heatmap.2(core_upregulated_kinetics_value_only, Rowv=TRUE, Colv=FALSE, dendrogram="row", symm=F, scale="column", na.rm=TRUE, col=my_palette, na.color="white", trace="none", key=T, keysize=5, density.info="none", key.xlab="Magnitude of expression changes", srtCol=45,  lmat = rbind(c(4,0),c(2,1),c(3,0)), lwid = c(1.5,4), lhei = c(1.5,4,1)) 
#If I want the horizontal key
heatmap.2(core_upregulated_kinetics_value_only, Rowv=TRUE, Colv=FALSE, dendrogram="row", symm=F, scale="column", na.rm=TRUE, col=my_palette, na.color="white", trace="none", key=T, keysize=1, density.info="none", key.xlab="Magnitude of expression changes", srtCol=45,  lmat = rbind(c(0,4),c(2,1),c(3,0)), lwid = c(1.5,4), lhei = c(1.5,4,1)) 


#####
core_downregulated_kinetics = get_Kinetics_with_directions("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_(7_or_more_conditions)/core_genes_7_plus_live_bacteria_downregulated_gene_list.txt")
rownames(core_downregulated_kinetics) = core_downregulated_kinetics$gene_name
core_downregulated_kinetics_value_only = as.matrix(core_downregulated_kinetics[,3:11])
core_downregulated_kinetics_value_only = core_downregulated_kinetics_value_only[,c(1,9,8,2,3,4,7,6,5)] #change the condition order
core_downregulated_kinetics_value_only = core_downregulated_kinetics_value_only[c("CR43051", "CG16926", "CG14629", "CG12091", "CG11200", "CG12338", "CG17124", "CG9664", "Non1", "Tps1", "Lectin-galC1", "alpha-Est2", "Try29F", "Cyp6a8", "Jheh1", "ninaD", "tobi", "CG42369", "CG34331", "CG34166", "GstE8", "CG18003", "CG31104", "CG31075", "CG3699", "Cyp6a18", "fit", "CG3301", "CG3739", "CG3940", "CG16986", "CG16985", "CG7997", "CG9510", "CG3609", "Mur18B", "CG15263", "glob1", "Cyp6g1", "Actbeta", "alpha-Est7", "Cyp4e2", "Spat", "Jon99Ci", "Lsp1beta", "Lsp2", "fat-spondin", "Ugt37b1", "CG13365", "CG3603", "CG2233", "CG15369", "regucalcin", "CG9512", "dob","Cyp6t1", "Iris", "Cyp4ac2", "CG17108", "Mal-B1", "CG2070", "CG16898", "dsb", "Cyp4d20", "CG4950", "CG18135", "CG5618", "Est-Q", "CG16749", "PKD", "CG10824", "CG10208", "to", "Obp99b", "Jon99Fii", "Npc2h", "lectin-33A", "lectin-28C","Ugt36Bc", "Nplp2", "CG9034", "CG31548", "CG33301", "26-29-p","CG42788", "Acsl"),] #manually change the gene order to group blues, reds, and the mix


#no scale
# my_palette = colorRampPalette(c("lightgrey","#dfc27d", "#a6611a"))(n = 299) 
# heatmap.2(as.matrix(core_downregulated_kinetics[,3:11]), density.info="none", dendrogram="row",trace="none", symm=F, scale="none", col=my_palette, srtCol=45, key=T, keysize = 0.5, cexCol=1, cexRow=0.8, key.xlab="Magnitude of expression changes", lmat = rbind(c(4,0),c(2,1),c(3,0)), lwid = c(1.5,4), lhei = c(1.5,4,1) ) 

#if scale by column (0 in the scale is raw number 1 = not recovered/similart to highest peak)
my_palette = colorRampPalette(c("#2166ac", "lightgrey", "#b2182b"))(n = 299) #colorblind-safe blue, grey, red
heatmap.2(core_downregulated_kinetics_value_only, Rowv=TRUE, Colv=FALSE, dendrogram="row", symm=F, scale="column", na.rm=TRUE, col=my_palette, na.color="white", trace="none", key=F, density.info="none", key.xlab="Magnitude of expression changes", srtCol=45, margin=c(6,6)) #doesn't work
heatmap.2(core_downregulated_kinetics_value_only, Rowv=FALSE, Colv=FALSE, dendrogram="none", symm=F, scale="column", na.rm=TRUE, col=my_palette, na.color="white", trace="none", key=F, density.info="none", key.xlab="Magnitude of expression changes", srtCol=45, margin=c(6,13), cexRow=1.3) #works
3
#just to get the key
heatmap.2(core_downregulated_kinetics_value_only, Rowv=FALSE, Colv=FALSE, dendrogram="none", symm=F, scale="column", na.rm=TRUE, col=my_palette, na.color="white", trace="none", key=T, keysize=2, density.info="none", key.xlab="Magnitude of expression changes", srtCol=45,  lmat = rbind(c(4,0),c(2,1),c(3,0)), lwid = c(1.5,4), lhei = c(1.5,4,1)) 
#getting a horizontal key
heatmap.2(core_downregulated_kinetics_value_only, Rowv=FALSE, Colv=FALSE, dendrogram="none", symm=F, scale="column", na.rm=TRUE, col=my_palette, na.color="white", trace="none", key=T, keysize=1, density.info="none", key.xlab="Magnitude of expression changes", srtCol=45,  lmat = rbind(c(0,4),c(2,1),c(3,0)), lwid = c(1.5,4), lhei = c(1.5,4,1)) 



#Updated the Trying something different (4/27/2017 at 4pm and 5/1/2017 at 10am, 5/3/2017, 5/4/2017)
#Treating the samples that have 132h > 12hr or 36h differently
#Making the ratio based off of the maximum value, not based off of 0

rm(list=ls(all=TRUE))
library(gplots); library(RColorBrewer); library(cluster)

#1) Filter out the EE-EE-EE genes
#2) Split the genes based on situation (12h>0, 36h>0, 132h>0)
#3) Remove any genes that have 12hr>0 and 36hr<0 OR 12hr<0 and 36hr>0 (decided after talking to Nicolas and Brian)

get_Kinetics_with_directions = function(core_genes_address){
    
    path_data_NS_removed = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.based.on.edgeR.results/recovery_genes/recovery_pattern_of_genes.txt", header=T)
    expression_data = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_basic_comparison_FDR_converted_to_Y-N_at_least_one_sig.txt", header=T) #fc
    colnames(expression_data)[1:2] = c("gene_id", "gene_name")
    core_regulated = read.table(core_genes_address, header=T)
    
    path_data_NS_removed_core_dir_only = path_data_NS_removed[which(path_data_NS_removed$gene.id %in% core_regulated$gene_id),]
    #path_data_NS_removed_core_dir_only= path_data_NS_removed[1:24,1:3] #for test run
    #path_data_NS_removed_core_dir_only[,3]= rep("Up-Up-Up",24) #for test run
    expression_data = expression_data[which(expression_data$gene_id %in% core_regulated$gene_id),]
    expression_data_fc_only = expression_data[,c(1:2,seq(3,43,2), seq(53,64,2))] #for fc
    
    #testset
    #expression_data_fc_only = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/clustering/clustering.based.on.edgeR.results/recovery_genes/kinetics_example_050417.txt", header=T)
    
    heatmap_number = matrix(NA, ncol = (dim(path_data_NS_removed_core_dir_only)[2]-2), nrow = dim(path_data_NS_removed_core_dir_only)[1])
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
                
                #cat("twelve: ", twelve, " thirtysix: ", thirtysix, " onethirtytwo: ", onethirtytwo, "\n")
                
                ratio=0; onethirtytwo_added =0
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
    heatmap_number_with_name = cbind(path_data_NS_removed_core_dir_only[,c(1:2)], heatmap_number)
    colnames(heatmap_number_with_name) = c("gene_id", "gene_name", "clean.prick","M.luteus","E.coli","S.marcescens","E.faecalis","P.rettgeri","Ecc15","E.fae.heatkilled", "P.rett.heatkilled")
    return(heatmap_number_with_name)
}

#Core upregulated genes:
core_upregulated_kinetics = get_Kinetics_with_directions("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_(7_or_more_conditions)/core_genes_7_plus_live_bacteria_upregulated_gene_list.txt")
rownames(core_upregulated_kinetics) = core_upregulated_kinetics$gene_name

core_upregulated_kinetics_value_only = as.matrix(core_upregulated_kinetics[,3:11])
core_upregulated_kinetics_value_only = core_upregulated_kinetics_value_only[-which(rownames(core_upregulated_kinetics_value_only) == c("CG32284")),] #remove the case that has NA in every single condition (165 genes), I need to remove this to get the clustering work
core_upregulated_kinetics_value_only = core_upregulated_kinetics_value_only[,c(1,9,8,2,3,4,7,6,5)] #change the order to reflect virulence

my_palette = colorRampPalette(c("#2166ac", "lightgrey", "#b2182b"))(n = 299) #colorblind-safe blue, grey, red
#my_palette = colorRampPalette(c("#efedf5", "#3f007d"))(n = 299) #colorblind-safe uni-directional purple
#Get the heatmap
heatmap.2(core_upregulated_kinetics_value_only, Rowv=TRUE, Colv=FALSE, dendrogram="none", symm=F, scale="column", na.rm=TRUE, col=my_palette, na.color="white", trace="none", key=F, density.info="none", key.xlab="Magnitude of expression changes", srtCol=45, margin=c(6,13))#, cexRow=0.8) #works

#Get the color key (square)
# heatmap.2(core_upregulated_kinetics_value_only, Rowv=TRUE, Colv=FALSE, dendrogram="row", symm=F, scale="column", na.rm=TRUE, col=my_palette, na.color="white", trace="none", key=T, keysize=5, density.info="none", key.xlab="Magnitude of expression changes", srtCol=45, lmat = rbind(c(4,0),c(2,1),c(3,0)), lwid = c(1.5,4), lhei = c(1.5,4,1))
# #Get the color key (rectangle)
# heatmap.2(core_upregulated_kinetics_value_only, Rowv=TRUE, Colv=FALSE, dendrogram="row", symm=F, scale="column", na.rm=TRUE, col=my_palette, na.color="white", trace="none", key=T, keysize=1, density.info="none", key.title="Spectrum", key.xlab="Spectrum of expression patterns", srtCol=45, lmat = rbind(c(0,4),c(2,1),c(3,0)), lwid = c(1.5,4), lhei = c(1,4,1)) 


#Core downregulated genes:
core_downregulated_kinetics = get_Kinetics_with_directions("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/finding.core.genes/true_core_genes_(7_or_more_conditions)/core_genes_7_plus_live_bacteria_downregulated_gene_list.txt")
rownames(core_downregulated_kinetics) = core_downregulated_kinetics$gene_name
core_downregulated_kinetics_value_only = as.matrix(core_downregulated_kinetics[,3:11])
core_downregulated_kinetics_value_only = core_downregulated_kinetics_value_only[,c(1,9,8,2,3,4,7,6,5)] #change the condition order
core_downregulated_kinetics_value_only = core_downregulated_kinetics_value_only[c("CG11200","PKD","CG9034","CR43051", "CG16926", "CG14629", "CG12091", "CG12338", "CG17124", "CG9664", "Non1", "Tps1", "Lectin-galC1", "alpha-Est2", "Try29F", "Cyp6a8", "Jheh1", "ninaD","CG34166","CG3739","CG3940","CG9510", "CG3609","Jon99Ci", "tobi", "CG42369", "CG34331",  "GstE8", "CG18003", "CG31104", "CG31075", "CG3699", "Cyp6a18", "fit", "CG3301", "CG16986", "CG16985", "CG7997", "Mur18B", "CG15263", "glob1", "Cyp6g1", "Actbeta", "alpha-Est7", "Cyp4e2", "Spat", "Lsp1beta", "Lsp2", "fat-spondin", "Ugt37b1", "CG13365", "CG3603", "CG2233", "CG15369", "regucalcin", "CG9512", "dob", "Iris",  "dsb", "Cyp4d20", "CG4950", "CG18135", "CG5618", "Est-Q", "CG16749",  "CG10824", "CG10208", "to", "Obp99b", "Jon99Fii", "lectin-28C", "Nplp2", "CG31548", "CG33301", "26-29-p","CG42788", "Acsl","Ugt36Bc","Npc2h", "lectin-33A","Cyp4ac2", "CG17108", "Mal-B1", "CG2070", "CG16898","Cyp6t1"),] #manually change the gene order to group blues, reds, and the mix

my_palette = colorRampPalette(c("#2166ac", "lightgrey", "#b2182b"))(n = 299) #colorblind-safe blue, grey, red
#my_palette = colorRampPalette(c("#efedf5", "#3f007d"))(n = 299) #colorblind-safe uni-directional purple
#Get the heatmap
heatmap.2(core_downregulated_kinetics_value_only, Rowv=FALSE, Colv=FALSE, dendrogram="none", symm=F, scale="column", na.rm=TRUE, col=my_palette, na.color="white", trace="none", key=F, density.info="none", srtCol=45, margin=c(6,13))#, cexRow=1.1) #works

#Get the color key (square)
# heatmap.2(core_downregulated_kinetics_value_only, Rowv=FALSE, Colv=FALSE, dendrogram="none", symm=F, scale="column", na.rm=TRUE, col=my_palette, na.color="white", trace="none", key=T, keysize=5, density.info="none", key.xlab="Magnitude of expression changes", srtCol=45, lmat = rbind(c(4,0),c(2,1),c(3,0)), lwid = c(1.5,4), lhei = c(1.5,4,1))
# #Get the color key (rectangle)
# heatmap.2(core_downregulated_kinetics_value_only, Rowv=FALSE, Colv=FALSE, dendrogram="none", symm=F, scale="column", na.rm=TRUE, col=my_palette, na.color="white", trace="none", key=T, keysize=1, density.info="none", key.title="Spectrum", key.xlab="Spectrum of expression patterns", srtCol=45, lmat = rbind(c(0,4),c(2,1),c(3,0)), lwid = c(1.5,4), lhei = c(1,4,1)) 


#checking colors for the heatmap
par(mfrow = c(2,2))
my_palette = colorRampPalette(c("#2166ac", "white", "#b2182b"))(n = 299) #colorblind-safe blue, grey, red
heatmap.2(core_downregulated_kinetics_value_only[c(3:5,7:9), c(4:6,8,9)], Rowv=FALSE, Colv=FALSE, dendrogram="none", symm=F, scale="column", na.rm=TRUE, col=my_palette, na.color="white", trace="none", key=T, keysize=1, density.info="none", key.title="Spectrum", key.xlab="#2166ac, white, #b2182b", srtCol=45, lmat = rbind(c(0,4),c(2,1),c(3,0)), lwid = c(1.5,4), lhei = c(1,4,1)) 

my_palette = colorRampPalette(c("#4575b4", "white", "#d73027"))(n = 299) #colorblind-safe blue, grey, red
heatmap.2(core_downregulated_kinetics_value_only[c(3:5,7:9), c(4:6,8,9)], Rowv=FALSE, Colv=FALSE, dendrogram="none", symm=F, scale="column", na.rm=TRUE, col=my_palette, na.color="white", trace="none", key=T, keysize=1, density.info="none", key.title="Spectrum", key.xlab="#4575b4, white, #d73027", srtCol=45, lmat = rbind(c(0,4),c(2,1),c(3,0)), lwid = c(1.5,4), lhei = c(1,4,1)) 

my_palette = colorRampPalette(c("yellow", "white", "blue"))(n = 299) #colorblind-safe blue, grey, red
heatmap.2(core_downregulated_kinetics_value_only[c(3:5,7:9), c(4:6,8,9)], Rowv=FALSE, Colv=FALSE, dendrogram="none", symm=F, scale="column", na.rm=TRUE, col=my_palette, na.color="white", trace="none", key=T, keysize=1, density.info="none", key.title="Spectrum", key.xlab="yellow, white, blue", srtCol=45, lmat = rbind(c(0,4),c(2,1),c(3,0)), lwid = c(1.5,4), lhei = c(1,4,1)) 

my_palette = colorRampPalette(c("yellow", "grey", "blue"))(n = 299) #colorblind-safe blue, grey, red
heatmap.2(core_downregulated_kinetics_value_only[c(3:5,7:9), c(4:6,8,9)], Rowv=FALSE, Colv=FALSE, dendrogram="none", symm=F, scale="column", na.rm=TRUE, col=my_palette, na.color="white", trace="none", key=T, keysize=1, density.info="none", key.title="Spectrum", key.xlab="yellow, grey, blue", srtCol=45, lmat = rbind(c(0,4),c(2,1),c(3,0)), lwid = c(1.5,4), lhei = c(1,4,1)) 

my_palette = colorRampPalette(c("#2166ac", "grey", "#b2182b"))(n = 299) #colorblind-safe blue, grey, red
heatmap.2(core_downregulated_kinetics_value_only[c(3:5,7:9), c(4:6,8,9)], Rowv=FALSE, Colv=FALSE, dendrogram="none", symm=F, scale="column", na.rm=TRUE, col=my_palette, na.color="white", trace="none", key=T, keysize=1, density.info="none", key.title="Spectrum", key.xlab="#2166ac, black, #b2182b", srtCol=45, lmat = rbind(c(0,4),c(2,1),c(3,0)), lwid = c(1.5,4), lhei = c(1,4,1)) 


#===
# #works, but old code:
# heatmap.2(core_upregulated_kinetics_value_only, Rowv=TRUE, Colv="NA", density.info="none", dendrogram="row", trace="none", symm=F, scale="column", col=my_palette, srtCol=45, key=T, keysize = 0.5, cexCol=1, cexRow=0.8, key.xlab="Magnitude of expression changes", lmat = rbind(c(4,0),c(2,1),c(3,0)), lwid = c(1.5,4), lhei = c(1.5,4,1) ) 