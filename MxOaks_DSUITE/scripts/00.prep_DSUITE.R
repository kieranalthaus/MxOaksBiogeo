library(ape)
library(geoscale)
library(phytools)
library(phyloch)

## 1. Load Data & Reformat -----------------------------------------------------
# Load tree
tree = read.tree(file = "DATING/out/20240321_p04_finaltree.tre")

# Replace tip labels with new labels
tree$tip.label[tree$tip.label == "Quercus_striatula|QUE002051"] = "Quercus_cf._striatula|QUE002051"
tree$tip.label[tree$tip.label == "Quercus_depressipes|QUE002026"] = "Quercus_cf._depressipes|QUE002026"
tree$tip.label[tree$tip.label == "Quercus_viminea|QUE001270"] = "Quercus_bolanyosensis|QUE001270"
tree$tip.label[tree$tip.label == "Quercus_crassifolia|QUE001276"] = "Quercus_crassifolia_s.l|QUE001276"

# Drop Quercus_crassifolia|QUE002021
tree = treeio::drop.tip(object = tree, tip = "Quercus_crassifolia|QUE002021")

# Drop Quercus_microphylla|QUE002858
tree = treeio::drop.tip(object = tree, tip = "Quercus_microphylla|QUE002858")

# Drop Quercus_carmenensis|QUE003185
tree = treeio::drop.tip(object = tree, tip = "Quercus_carmenensis|QUE003185")

# Drop Quercus_carmenensis|QUE003179
tree = treeio::drop.tip(object = tree, tip = "Quercus_carmenensis|QUE003179")

# Drop Quercus_arizonica|QUE002003
tree = treeio::drop.tip(object = tree, tip = "Quercus_arizonica|QUE002003")

# Drop Quercus_laeta|QUE001323
tree = treeio::drop.tip(object = tree, tip = "Quercus_laeta|QUE001323")

# Drop Quercus_potosina|QUE000205
tree = treeio::drop.tip(object = tree, tip = "Quercus_potosina|QUE000205")

# Remove duplicated tips
tree = treeio::drop.tip(object = tree, tip = "Quercus_greggii|QUE002745")

# Drop fussy species
tree = treeio::drop.tip(object = tree, tip = "Quercus_salicifolia|QUE002920") # A randomly chosen of the two Q. salicifolia

tree = treeio::drop.tip(object = tree, tip = "Quercus_cf._miquihuanensis-_autopista|QUE002738")


# List of Old World oaks that are not to be included in the final analysis
bad.oaks = c("Quercus_robur_ssp._imeretina|QUE001600",
             "Quercus_hartwissiana|QUE001505",
             "Quercus_robur_ssp._pedunculiflora|QUE000609", 
             "Quercus_robur|QUE001514", 
             "Quercus_canariensis|QUE000325", 
             "Quercus_petraea_ssp._iberica|QUE000113", 
             "Quercus_dalechampii|QUE000658", 
             "Quercus_petraea|QUE001513", 
             "Quercus_lusitanica|QUE001506", 
             "Quercus_pyrenaica|QUE000579", 
             "Quercus_frainetto|QUE001509", 
             "Quercus_pubescens|QUE001986", 
             "Quercus_faginea|QUE000993", 
             "Quercus_infectoria|QUE001510", 
             "Quercus_infectoria_subsp._veneris|QUE002906",
             "Quercus_boissieri|QUE000995", 
             "Quercus_cedrorum|QUE001563",
             "Quercus_kotschyana|QUE001560", 
             "Quercus_petraea_ssp._iberica|QUE001597", 
             "Quercus_macranthera|QUE001599", 
             "Quercus_vulcanica|QUE001508", 
             "Quercus_griffithii|QUE000829", 
             "Quercus_aliena|QUE000923", 
             "Quercus_fabri|QUE000814",
             "Quercus_aliena_var._acutiserrata_|QUE002903", 
             "Quercus_yunnanensis|QUE000924", 
             "Quercus_dentata|QUE001007", 
             "Quercus_mongolica|QUE000148", 
             "Quercus_serrata|QUE001059", 
             "Quercus_mongolica_var._grosseserrata|QUE002913", 
             "Quercus_liaotungensis|QUE002912",
             'Quercus_senescens|QUE000917',
             'Quercus_rehderiana|QUE000921',
             'Quercus_rehderiana|QUE001495',
             'Quercus_monimotricha|QUE000945',
             'Quercus_pannosa|QUE004304',
             'Quercus_gujavifolia|QUE001515',
             'Quercus_longispica|QUE000919',
             'Quercus_spinosa|QUE004306',
             'Quercus_sp._nov.|QUE001568',
             'Quercus_semecarpifolia|QUE001494',
             'Quercus_floribunda|QUE001493',
             'Quercus_marlipoensis|QUE004303',
             'Quercus_engleriana|QUE000935',
             'Quercus_utilis|QUE000937',
             'Quercus_setulosa|QUE001574',
             'Quercus_bawangliensis|QUE004296',
             'Quercus_tarokoensis|QUE004307',
             'Quercus_yiwuensis|QUE000936',
             'Quercus_cocciferoides|QUE004297',
             'Quercus_kingiana|QUE004299',
             'Quercus_coccifera|QUE001988',
             'Quercus_calliprinos|QUE002667',
             'Quercus_aucheri|QUE001497',
             'Quercus_alnifolia|QUE000999',
             'Quercus_rotundifolia|QUE002755',
             'Quercus_ilex_ssp._rotundifolia|QUE001498',
             'Quercus_ilex|QUE002668',
             'Quercus_baloot|QUE001492',
             'Quercus_dolicholepis|QUE000944',
             'Quercus_baronii|QUE001064',
             'Quercus_phillyreoides|QUE000925',
             'Quercus_acrodonta|QUE000918',
             'Quercus_leucotrichophora|QUE004301',
             'Quercus_lanata|QUE003056',
             'Quercus_franchetii|QUE000922',
             'Quercus_look|QUE000983',
             'Quercus_cerris|QUE001564',
             'Quercus_castaneifolia|QUE001502',
             'Quercus_trojana|QUE001011',
             'Quercus_libani|QUE001000',
             'Quercus_afares|QUE001501',
             'Quercus_macrolepis|QUE000998',
             'Quercus_brantii|QUE000991',
             'Quercus_ithaburensis|QUE001103',
             'Quercus_suber|QUE001601',
             'Quercus_crenata|QUE000954',
             'Quercus_variabilis|QUE000920',
             'Quercus_acutissima|QUE001003',
             'Quercus_chenii|QUE001500',
             'Quercus_blakei|QUE000899',
             'Quercus_bella|QUE000930',
             'Quercus_phanera|QUE000938',
             'Quercus_chapensis|QUE000909',
             'Quercus_litoralis|QUE000943',
             'Quercus_fleuryi|QUE000904',
             'Quercus_langbianensis|QUE000940',
             'Quercus_pachyloma|QUE000905',
             'Quercus_poilanei|QUE000939',
             'Quercus_daimingshanensis|QUE000933',
             'Quercus_patelliformis|QUE000906',
             'Quercus_ciliaris|QUE000896',
             'Quercus_arbutifolia|QUE000901',
             'Quercus_stewardiana|QUE000931',
             'Quercus_kiukiangensis|QUE004300',
             'Quercus_sessilifolia|QUE000912',
             'Quercus_acuta|QUE000932',
             'Quercus_chrysocalyx|QUE000914',
             'Quercus_annulata|QUE000915',
             'Quercus_oxyodon|QUE000900',
             'Quercus_schottkyana|QUE000903',
             'Quercus_lineata|QUE004302',
             'Quercus_salicina|QUE002669',
             'Quercus_multinervis|QUE000902',
             'Quercus_myrsinifolia|QUE000911',
             'Quercus_kouangsiensis|QUE000913',
             'Quercus_lamellosa|QUE000908',
             'Quercus_eumorpha|QUE004298',
             'Quercus_augustini|QUE000941',
             'Quercus_jenseniana|QUE001573',
             'Quercus_delavayi|QUE000907',
             'Quercus_chungii|QUE001572',
             'Quercus_kerrii|QUE000897',
             'Quercus_austrocochinchinensis|QUE000928',
             'Quercus_rex|QUE000898',
             'Quercus_gilva|QUE000934',
             'Quercus_championii|QUE000910',
             'Quercus_sichourensis|QUE004305',
             'Notholithocarpus_densiflorus|QUE001252',
             'Lithocarpus_mairei|QUE000950',
             'Lithocarpus_hancei|QUE000951',
             'Lithocarpus_litseifolius|QUE000949',
             'Lithocarpus_longinux|QUE000926',
             'Chrysolepis_chrysophylla|QUE001253',
             'Castanea_mollissima|QUE002799',
             'Castanea_henryi|QUE002802',
             'Castanea_dentata|QUE000948',
             # 'Castanopsis_fissa|QUE000927',
             "Quercus_sinuata_var._breviloba|QUE000665",
             "Quercus_greggii|QUE002737",
             "Quercus_new_species|QUE000227",
             "Quercus_sp.|QUE002663",
             "Quercus_pontica|QUE000601"
)

# Drop tips
tree = treeio::drop.tip(object = tree, tip = bad.oaks)   

## 2. Create data subsets ------------------------------------------------------
# For red oaks
rednode = findMRCA(tree = tree, 
                   tips = c("Quercus_agrifolia|QUE001242", "Quercus_radiata|QUE001304"))
# For white oaks
whitenode = findMRCA(tree = tree, 
                     tips = c("Quercus_cedrosensis|QUE001187", "Quercus_turbinella|QUE000176"))

# Get taxa subtending each node
white_taxa = descendants(phy = tree, node = whitenode, labels = TRUE)
red_taxa = descendants(phy = tree, node = rednode, labels = TRUE)

# Create subtrees of red and white oaks (+ Castanopsis_fissa|QUE000927)
red.tree = drop.tip(phy = tree, tip = white_taxa)
white.tree = drop.tip(phy = tree, tip = red_taxa)

# Export trees
write.tree(phy = red.tree, file = "DSUITE/data/dsuite_red.tree")
write.tree(phy = white.tree, file = "DSUITE/data/dsuite_white.tree")

## Create dataset to filter VCF file
# For red oaks
write.table(matrix(c(red_taxa, "Castanopsis_fissa|QUE000927", red_taxa, "Outgroup"), 
                   ncol = 2),
            file = "DSUITE/data/dsuite_red_filter.txt", 
            quote = F,
            row.names = F,
            col.names = F,
            sep = "\t")

# For white oaks
write.table(matrix(c(white_taxa, "Castanopsis_fissa|QUE000927", white_taxa, "Outgroup"), 
                   ncol = 2), 
            file = "DSUITE/data/dsuite_white_filter.txt", 
            quote = F,
            row.names = F,
            sep = "\t")
















