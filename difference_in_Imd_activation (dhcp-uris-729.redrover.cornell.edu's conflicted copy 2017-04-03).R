#Difference in the IMD activation
#February 18th, 2017
#Joo Hyun Im (ji72)

rm(list=ls(all=TRUE)) #delete any previous entry

###Option 1 (qualitative): Create a list of Imd genes and gather information on whether a gene in a particular condition had its expression level return to the base. 

#script from finding_recovery_genes.R

#Get the genes with excluding those that were EE-EE-EE in every condition.
#Read in data
original_data = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015/edgeR_basic_comparison_FDR_converted_to_Y-N_at_least_one_sig.txt", header=T)

#Assign the direction
assign_direction = function(data){
    length.of.table <- as.numeric(dim(data)[1])
    width.of.table <- as.numeric(dim(data)[2])
    indicator <- NULL
    
    for (i in 1:length.of.table){
        for (s in seq(3, width.of.table, 2)){
            
            if (data[i,s+1] == "Y"){ #if significant
                indicator <- isTRUE(data[i,s] > 0) #indicator shows that the direction is positive
                if (indicator == TRUE) { #If the case is Up-DEG
                    data[i,s] = "Up"
                }
                else { #If the caseis Down-DEG
                    data[i,s] = "Down"
                }
            }
            else { #if not significant
                data[i,s] = "EE"
            }
            
        }
    }
    return(data)
}
direction_data = assign_direction(original_data)

#Put the expression path together
put_paths_together = function(data){
    full.expression.path.table = matrix(NA, nrow = dim(data)[1], ncol= 9)
    for (b in 1:dim(data)[1]){ #for each gene
        count=1
        for (c in c(3,9,15,21,27,33,39,53,59)) { #for each condition of cleanprick, 6 live conditions and 2 heatkilled conditions 
            start = c; end = c+5
            subset = data[b,start:end]; subset
            path = paste(subset[1,1],"-",subset[1,3],"-",subset[1,5], sep="")
            full.expression.path.table[b, count] = path
            count=count+1
        }
    }
    return(full.expression.path.table) 
}
path_data = put_paths_together(direction_data)
path_data = cbind(direction_data[,c(1:2)], path_data)

#Remove any genes that havent' changed significantly in any of the conditions.
remove_NS = function(data){
    full.expression.path.table.EEs.removed = c()
    for (k in 1:dim(data)[1]){ #1, 2, 3, ... 11911
        for (m in seq(3, 11, 1)){ #1, 4, ..., 9
            indicator <- isTRUE(data[k,m] != "EE-EE-EE") #When the gene is up or down
            if (indicator == T){
                full.expression.path.table.EEs.removed = rbind(full.expression.path.table.EEs.removed, data[k,])
                break
            }
        }
    }
    return(full.expression.path.table.EEs.removed)
}
path_data_NS_removed = remove_NS(path_data) #1981 entries
colnames(path_data_NS_removed) = c("gene.id","gene.name","clean.prick", "M.luteus", "E.coli", "S.marcescens", "E.faecalis","P.rettgeri","Ecc15", "E.fae.heatkilled","P.rett.heatkilled")

#Only get the Imd genes
imd_genes = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/Imd_genes_curated_by_Early_et_al_2017.txt", header=T, sep="\t")
imd_genes_id = imd_genes[,1]
intersection_id = path_data_NS_removed[path_data_NS_removed$gene.id %in% imd_genes_id,]; dim(intersection_id) #31 -- also tried with name but this worked better
#imd_genes_name = imd_genes[,3]; intersection_name = path_data_NS_removed[path_data_NS_removed$gene.name %in% imd_genes_name,]; dim(intersection_name) #13

#Change the expression paths to 'Y' or 'N' depending on the status of recovery
convert_path_to_direction = function(input_data){
    for (n in 3:11){
        input_data[,n] = as.character(input_data[,n])
    }
    path_converted_to_direction = input_data
    possible_patterns = c("Up-Up-EE","Up-EE-EE","Up-Down-EE","EE-Up-EE","EE-Down-EE","Down-Up-EE","Down-EE-EE", "Down-Down-EE")
    
    for (j in 1:dim(input_data)[1]){
        for (m in 3:11){ #excluding the gene id and gene name
            pattern = as.character(input_data[j,m])
            if (pattern %in% possible_patterns){ #if the pattern is one of the recovery patterns
                path_converted_to_direction[j,m] = toString("Y")
            }
            else{
                path_converted_to_direction[j,m] = toString("N")
            }
        }
    }
    return(path_converted_to_direction)
}
Imd_genes_with_expression_path_direction_table = convert_path_to_direction(intersection_id)
write.table(Imd_genes_with_expression_path_direction_table, file = "/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/PCA_using_DESeq2/difference_in_Imd_activation_based_on_expression_path_recovery_status.txt", quote=F, col.names = T, row.names = F)


###Option 2 (quantitative): Create a list of Imd genes and create a heatmap of fold change (0h-12hr) and see if the responses differ by conditions.

#script from pulling.out.genes.only.activated.in.one.bacteria.R

#Get the relevant data:
rm(list=ls(all=TRUE)) #delete any previous entry
setwd("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/edgeR_results_with_cpm_filter_all_UCs_Nov_2015")
data_table <- read.table("edgeR_basic_comparison_pval_at_least_one_sig.txt",header=T) #FC data, not counts.

odd_index = seq(3,64,2)
data_table_FC_only = data_table[,c(odd_index)]; rownames(data_table_FC_only) = data_table[,1]
data_table_FC_only_12hr = data_table_FC_only[,c(1,4,7,10,13,16,19,22,23,24,25,26,29)] #12hr sample only

#Get Imd genes
imd_genes = read.table("/Users/JooHyun/Dropbox/Cornell/Lab/Projects/Mega_RNA-seq/gene_sets/Imd_genes_curated_by_Early_et_al_2017.txt", header=T, sep="\t")
imd_genes_id = imd_genes[,1]
data_table_FC_only_12hr = data_table_FC_only_12hr[rownames(data_table_FC_only_12hr) %in% imd_genes_id,]; dim(data_table_FC_only_12hr) #33 genes
colnames(data_table_FC_only_12hr) = c("SterileWound", "M.luteus", "E.coli", "S.marcescens Type", "E.faecalis live", "P.rettgeri live","Ecc15","S.aureus","P.sneebia","S.marcescens Db11", "P.entomophila", "E.faecalis heatkilled","P.rettgeri heatkilled")
imd_genes_names = imd_genes[match(rownames(data_table_FC_only_12hr),imd_genes$gene_id),]
data_table_FC_only_12hr = cbind(data_table_FC_only_12hr,imd_genes_names[,3], imd_genes_names[,6])
function_order = c("Recognition", "Signaling", "Signaling: JNK", "AMP", "JAK-STAT")
data_table_FC_only_12hr = data_table_FC_only_12hr[order(match(data_table_FC_only_12hr$`imd_genes_names[, 6]`, function_order)),]
rownames(data_table_FC_only_12hr) = data_table_FC_only_12hr[,14]
data_table_FC_only_12hr = data_table_FC_only_12hr[,1:13]

#Draw a heatmap
library(gplots); library(RColorBrewer)
my_palette = colorRampPalette(c("#67a9cf","#f7f7f7","#ef8a62"))(n = 299)
#scale the log2FC
#heatmap.2(as.matrix(data_table_FC_only_12hr), Rowv=FALSE, density.info="none", dendrogram="column", trace="none",symm=F, scale=c("column"), col=my_palette, margins = c(12,12), srtCol=45)

#do not scale the log2FC
other_colors = brewer.pal(5, "Set2")
table(imd_genes[,4])
Label = c(rep(other_colors[1],5),rep(other_colors[2],9),rep(other_colors[3],5),rep(other_colors[4],11),rep(other_colors[5],3))

#pdf(file="/Users/JooHyun/Desktop/t.pdf", height=10, width=10) 
heatmap.2(as.matrix(data_table_FC_only_12hr), Rowv=FALSE, density.info="none", dendrogram="column", trace="none",symm=F, scale="none", col=my_palette, margins = c(12,12), srtCol=45, RowSideColors = Label, keysize = 1,key.xlab="Expression")
#dev.off()


