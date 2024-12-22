# 4 variants: Delta, Omicron BA_one, Omicron BA_two, Others
colors_lineage <- RColorBrewer::brewer.pal(4, "Set1")
variants = c("Delta", "Omicron BA.1", "Omicron BA.2", "Others")
names(colors_lineage) <- variants

colors_lineage_old <- colors_lineage
names(colors_lineage_old) <- c("Delta", "Omicron BA_one", "Omicron BA_two", "Others")

# Imported/Community: Dark/Light
colors_lineage_type <- c("#9f1214","#ed5e5f",
                        "#25567d","#69a3d2",
                        "#357933","#80c87d",
                        "#68356f","#b87dc1")
names(colors_lineage_type) <- c("Delta (imported)", "Delta (community)",
                      "Omicron BA.1 (imported)", "Omicron BA.1 (community)",
                      "Omicron BA.2 (imported)", "Omicron BA.2 (community)",
                      "Others (imported)", "Others (community)")

# Infected/Detected/Sequenced:
colors_cases_level <- c("#FC8D59","#FFFFBF","#91CF60")
names(colors_cases_level) <- c("Infected","Detected","Sequenced")

# 29 spatial units:
data_fitting <- readRDS(paste0("results/model_data/data_fitting_", "Omicron20", ".rds"))
source("scripts/model_fitting/helper/geo_grid.R")
df_layout_geo_grid_29 <- geo_grid_29()
df_layout_geo_grid_29 <- df_layout_geo_grid_29 %>% left_join(data_fitting$df_new_code %>% select(code=new_code, continent) %>% unique(), by="code") %>% mutate(n_char=nchar(code), continent=factor(continent, levels=c("NA", "EU", "AS", "SA", "AF", "OC"))) %>% arrange(continent, n_char, col, row)

# num_units_continents <- table(df_layout_geo_grid_29$continent)
# colors_af <- MetBrewer::met.brewer("Lakota", num_units_continents[names(num_units_continents)=="AF"], "discrete")
# colors_as <- MetBrewer::met.brewer("Austria", num_units_continents[names(num_units_continents)=="AS"], "discrete")
# colors_eu <- MetBrewer::met.brewer("Klimt", num_units_continents[names(num_units_continents)=="EU"], "discrete")
# colors_na <- MetBrewer::met.brewer("Navajo", num_units_continents[names(num_units_continents)=="NA"], "discrete")
# colors_oc <- MetBrewer::met.brewer("Egypt", num_units_continents[names(num_units_continents)=="OC"], "discrete")
# colors_sa <- MetBrewer::met.brewer("Kandinsky", num_units_continents[names(num_units_continents)=="SA"], "discrete")
# colors_spatial_units <- c(colors_af, colors_as, colors_eu, colors_na, colors_oc, colors_sa)

# ## evenly generate colors for continents
# colors_spatial_units <- MetBrewer::met.brewer("Austria", max(table(df_layout_geo_grid_29$continent))*length(table(df_layout_geo_grid_29$continent)), "continuous", override.order = TRUE)
# ## remove unused colors in continents
# id_to_rm <- sapply(seq_along(table(df_layout_geo_grid_29$continent)), function(i){
#   n_unused_i <- max(table(df_layout_geo_grid_29$continent)) - table(df_layout_geo_grid_29$continent)[i]
#   if(n_unused_i==0) return(NULL)
#   id_unused_start <- (i-1)*max(table(df_layout_geo_grid_29$continent))+1+table(df_layout_geo_grid_29$continent)[i]
#   id_unused_end <- id_unused_start + n_unused_i - 1
#   return(id_unused_start:id_unused_end)
# })
# colors_spatial_units <- colors_spatial_units[-unlist(id_to_rm)]
# names(colors_spatial_units) <- df_layout_geo_grid_29$code
# colors_spatial_units <- colors_spatial_units[order(names(colors_spatial_units))]

# The above code results in the below colors
#        AE AF_others        AR AS_others        AU        BR 
# "#318D29" "#883058" "#EFB52E" "#ECC714" "#642E4A" "#E5A540" 
#        CA        CL        DE        ES EU_others        FR 
# "#A40000" "#F9C51B" "#0B5458" "#0F4964" "#085F4D" "#123E6F" 
#        GB        HK        IN        JP        KR        MA 
# "#15327B" "#C7BB18" "#579825" "#7CA420" "#A1B01C" "#A64F7D" 
#        MX NA_others        NG        NZ OC_others        PG 
# "#7A0E24" "#651536" "#B15989" "#425B69" "#327278" "#53445A" 
# SA_others        TN        TR        US        ZA 
# "#DA9553" "#9C4571" "#0C812D" "#8F0712" "#923A64" 
colors_spatial_units <- c("#318D29", "#883058", "#EFB52E", "#ECC714", "#642E4A", "#E5A540", "#A40000", "#F9C51B", "#0B5458", "#0F4964", "#085F4D", "#123E6F", "#15327B", "#C7BB18", "#579825", "#7CA420", "#A1B01C", "#A64F7D", "#7A0E24", "#651536", "#B15989", "#425B69", "#327278", "#53445A", "#DA9553", "#9C4571", "#0C812D", "#8F0712", "#923A64")

# colors_spatial_units <- scico::scico(nrow(df_layout_geo_grid_29), palette = 'batlow')
names(colors_spatial_units) <- data_fitting$country_under_investigation

# 5 continents: Africa, Asia, Europe, North America, Oceania
library(MetBrewer)
colors_continent <- met.brewer("Wissing", 6)
names(colors_continent) <- c("Africa", "Asia", "Europe", "North America", "South America", "Oceania")

# Fig5abc
colors_fig5abc <- c("#8a3324", "#997a8d", "#b29882", "#666e57")

# IDR, DSR and Travel Level
colors_IDR_DSR_TraLevel <- c("#8b7764", "#65927a", "#f99b95")
names(colors_IDR_DSR_TraLevel) <- c("IDR", "DSR", "Travel control level")
