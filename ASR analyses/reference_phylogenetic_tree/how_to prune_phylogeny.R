library(ape)
setwd("~/Desktop/Irene_script/reference_phylogenetic_tree")



trees=read.tree("languages_tree.originals.tre")

spatialRate_tree = trees[[1]]
cognateRate_tree = trees[[2]]

to_keep<- c("Danish", "Faroese", "Icelandic_ST", "Riksmal", "Swedish_List", "Old_Norse" )


plot(spatialRate_tree)
spatialRate_tree_sub<-keep.tip(spatialRate_tree, tip = to_keep)
plot(spatialRate_tree_sub)
write.tree(spatialRate_tree_sub, file="spatialRate_subtree.tre")

plot(cognateRate_tree)
cognateRate_tree_sub<-keep.tip(cognateRate_tree, tip = to_keep)
plot(cognateRate_tree_sub)
write.tree(spatialRate_tree_sub, file="cognateRate_subtree.tre")


# After pruning the trees. 
1. change the name of the tips, replace COUNTRY NAMES by NUMBERS, in a txt reader. 
  1 Danish,
  2 Faroese,
  3 Swedish_List,
  4 Icelandic_ST,
  5 Riksmal,
  6 Old_Norse;

2. copier the tree in the python script: 
    change this line: outfile.write("\tTREE pruned_tree_nordic_languages = ((6:176.5572341,(2:354.7826639,4:354.7826639):581.4051686):376.9084392,((5:336.8580069,1:336.8580069):385.5014575,3:722.3594644):590.7368073);;\n")
    change only that part: ((6:176.5572341,(2:354.7826639,4:354.7826639):581.4051686):376.9084392,((5:336.8580069,1:336.8580069):385.5014575,3:722.3594644):590.7368073)
    
Done for the tree. 




