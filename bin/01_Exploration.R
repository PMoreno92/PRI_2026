source("bin/00_Packages.R")

# Data ----
hb_obs = lapply(list.files(path = "data/", full.names = T), function (x) read.csv(x))
names(hb_obs) = gsub("identified_plants_|\\.csv", "", list.files(path = "data/"))

hb_obs = lapply(1:length(hb_obs), function (x) cbind(hb_obs[[x]], hbSpecies = names(hb_obs)[x]))

hb_numvar = c("gbifIDPlant", "scorePlant", "decimalLatitude", "decimalLongitude", "gbifIDBird", 
              "coordinateUncertaintyInMeters", "year", "month", "day")

hb_df = data.table::rbindlist(hb_obs)|>
  mutate(eventDate = na_if(eventDate, "nan"))|>
  mutate(eventDate = ifelse(grepl("/", eventDate) & !grepl("-", eventDate),
                            format(as.Date(eventDate, format = "%d/%m/%Y"), "%Y-%m-%d"),
                            eventDate))|>
  mutate(year = ifelse(is.na(year), format(as.Date(eventDate),"%Y"), year),
         month = ifelse(is.na(month), format(as.Date(eventDate),"%m"), month),
         day = ifelse(is.na(day), format(as.Date(eventDate),"%d"), day))|>
  mutate_at(hb_numvar, as.numeric) # Problematic date formats: 149163, 149164, 149165, 149166 (I changed the eventDate BUT may need to change also the identificationDateTime)

## Solving problems with species (added at family level in the phylo tree)
hb_df$hbSpecies = gsub("Acestrura_heliodor", "Chaetocercus_heliodor", hb_df$hbSpecies)
hb_df$hbSpecies = gsub("Clytolaema_rubricauda", "Heliodoxa_rubricauda", hb_df$hbSpecies)
hb_df$hbSpecies = gsub("Leucolia_violiceps", "Ramosomyia_violiceps", hb_df$hbSpecies)
hb_df$hbSpecies = gsub("Leucolia_viridifrons", "Ramosomyia_viridifrons", hb_df$hbSpecies)

### Solving problems with species (added at genus level in the phylo tree)
hb_df$hbSpecies = gsub("Amazilia_amabilis", "Polyerata_amabilis", hb_df$hbSpecies)
hb_df$hbSpecies = gsub("Amazilia_boucardi", "Chrysuronia_boucardi", hb_df$hbSpecies)
hb_df$hbSpecies = gsub("Amazilia_chionogaster", "Elliotomyia_chionogaster", hb_df$hbSpecies)
hb_df$hbSpecies = gsub("Amazilia_fimbriata", "Chionomesa_fimbriata", hb_df$hbSpecies)
hb_df$hbSpecies = gsub("Amazilia_franciae", "Uranomitra_franciae", hb_df$hbSpecies)
hb_df$hbSpecies = gsub("Amazilia_saucerottei", "Saucerottia_saucerottei", hb_df$hbSpecies)
hb_df$hbSpecies = gsub("Amazilia_versicolor", "Chrysuronia_versicolor", hb_df$hbSpecies)
hb_df$hbSpecies = gsub("Amazilia_violiceps", "Ramosomyia_violiceps", hb_df$hbSpecies)
hb_df$hbSpecies = gsub("Calliphlox_bryantae", "Philodice_bryantae", hb_df$hbSpecies)
hb_df$hbSpecies = gsub("Calliphlox_mitchellii", "Philodice_mitchellii", hb_df$hbSpecies)
hb_df$hbSpecies = gsub("Chrysuronia_boucardi", "Amazilia_boucardi", hb_df$hbSpecies)
hb_df$hbSpecies = gsub("Hylocharis_cyanus", "Chlorestes_cyanus", hb_df$hbSpecies)
hb_df$hbSpecies = gsub("Hylocharis_eliciae", "Chlorestes_eliciae", hb_df$hbSpecies)
hb_df$hbSpecies = gsub("Lampornis_cinereicauda", "Lampornis_castaneoventris", hb_df$hbSpecies)
hb_df$hbSpecies = gsub("Philodice_evelynae", "Nesophlox_evelynae", hb_df$hbSpecies)

## Selected area
WH_territories = c("Aruba", "Anguilla", "Argentina", "Antigua", "Barbuda",
                   "Bahamas", "Saint Barthelemy", "Belize", "Bermuda", "Bolivia",
                   "Brazil", "Barbados", "Canada", "Chile", "Colombia", "Costa Rica", 
                   "Cuba", "Curacao", "Cayman Islands", "Dominica", "Dominican Republic", 
                   "Ecuador", "French Guiana", "Martinique", "Guadeloupe", 
                   "Grenada", "Guatemala", "Guyana", "Honduras", "Haiti", 
                   "Jamaica", "Nevis", "Saint Kitts", "Saint Lucia", "Saint Martin", 
                   "Mexico" , "Montserrat", "Nicaragua", "Bonaire", "Sint Eustatius", 
                   "Saba", "Panama", "Peru", "Puerto Rico", "Paraguay", "El Salvador", 
                   "Saint Pierre and Miquelon", "Suriname", "Sint Maarten", 
                   "Turks and Caicos Islands", "Trinidad", "Tobago", "Uruguay", 
                   "USA", "Grenadines", "Saint Vincent", "Venezuela", "Virgin Islands")

# Checks
## check if all important data is conserved
sum(is.na(hb_df$imageFileName))
sum(is.na(hb_df$speciesPlant))
sum(is.na(hb_df$scorePlant))
sum(is.na(hb_df$hbSpecies))

## Plot 
world_coordinates = map_data("world")

whem_coordinates = world_coordinates|>
  filter(region %in% WH_territories)

ggplot() +
  geom_map(
    data = whem_coordinates, map = whem_coordinates,
    aes(long, lat, map_id = region), fill = "gray") +
  geom_point(data = hb_df,
             aes(decimalLongitude, decimalLatitude), color = "red", alpha = 0.3) +
  xlim(c(-180, 0)) +
  theme_bw()

# Add the cells and coordinates (HERE I need to allow the grain size to change)
world_sf = maps::map('world', fill = TRUE, col = 1:10, plot = F)
whem_sf = st_as_sf(world_sf)|>
  filter(ID %in% WH_territories)|>
  sf::st_union()|>
  sf::st_crop(xmin = -180, xmax = 0, ymin = -90, ymax = 90)|>
  tidyterra::as_spatvector()

plot(whem_sf)

whem_sf_AEAP = terra::project(whem_sf, crs("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 
                                 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")) # Transform to Albers Equal Area Projection
whem_sf_EEA = terra::project(whem_sf, crs("EPSG:8858")) # Transform to Equal Earth Americas
plot(whem_sf_EEA)


## Elevation
elev = get_wc2.1(var = "elev", res = 2.5, path = "results") # Getting elevation data worldwide (WorldClim) (the date is given by month)
selelev_mean = terra::crop(x = elev, whem_sf) # Crop elevation to area of study
plot(selelev_mean)
###writeRaster(selelev_mean, "results/selelev_mean.tiff")

selelev_mean_EEA = project(selelev_mean, rast(whem_sf_EEA, resolution = 10000)) #Transform to Equal Earth Americas
plot(selelev_mean_EEA)

hb_df = cbind(hb_df,
               get_nth_cell_xy(coord_df = hb_df[,c("decimalLongitude","decimalLatitude")], 
                                   ras_file = selelev_mean_EEA))|>
  as.data.frame()


hb_df = hb_df[, -which(colnames(hb_df) %in% c("decimalLongitude","decimalLatitude"))[c(1,2)]]

# Choosing the plant score
hb_summ_scPl = data.table::rbindlist(lapply(seq(0, 95, 5), function (x)
  cbind(percent = x,
        psych::describe(hb_df|>
                          mutate(scorePlant = as.numeric(scorePlant))|>
                          filter(scorePlant > x)|>
                          count(hbSpecies)|>
                          dplyr::select(n)))))


plot(hb_summ_scPl$median~hb_summ_scPl$percent)

plot_scPl_median = ggplot(hb_summ_scPl, aes(percent, median)) +
  # Individual points and trends per site
  #geom_point(data = hb_summ_scPl, aes(percent, median), color = "#FF6970",
  #           size = 5, alpha = 1) +
  geom_line(data = hb_summ_scPl, aes(percent, median), color = "#FF6970",
             linewidth = 1.5, alpha = 1) +
  
  scale_x_continuous(breaks = seq(0, 90, 10), minor_breaks = seq(5, 95, 10), limits = c(0, 100), expand = c(0, 0)) +
  labs(y = expression(Median~Spp.~Obs.),
       x =  expression(AI~CI)) +
  theme(panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(0.25, 0, 0.25, 0.25, "cm"),
        axis.line.x = element_line(color = "black"), 
        axis.line.y = element_line(color = "black"), 
        text = element_text(size = 45),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        legend.position = "none")

plot_scPl_spp = ggplot(hb_summ_scPl, aes(percent, n)) +
  # Individual points and trends per site
  #geom_point(data = hb_summ_scPl, aes(percent, n), color = "#00C487",
  #           size = 5, alpha = 1) +
  geom_line(data = hb_summ_scPl, aes(percent, n), color = "#00C487",
            linewidth = 1.5, alpha = 1) +
  
  labs(y = expression(Number~Spp.~Obs.),
       x =  expression(AI~CI)) +
  scale_x_continuous(breaks = seq(0, 90, 10), minor_breaks = seq(5, 95, 10), limits = c(0, 100), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 400)) +
  theme(panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(0.25, 0, 0.25, 0.25, "cm"),
        axis.line.x = element_line(color = "black"), 
        axis.line.y = element_line(color = "black"), 
        text = element_text(size = 45),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        #axis.text.y = element_text(color = "#00C487"),
        legend.position = "none")

plots_scPl = cowplot::plot_grid(plot_scPl_spp, plot_scPl_median, nrow = 1, align = "hv")

# Choosing AI CI = 0.8
hb_80_df = hb_df|>
  mutate(scorePlant = as.numeric(scorePlant),
         decimalLongitude = as.numeric(decimalLongitude),
         decimalLatitude = as.numeric(decimalLatitude))|>
  filter(scorePlant > 80)

# Diversity ----
hb_abun <- table(hb_df$nth_cell, hb_df$hbSpecies)
hb_abuninv <- t(hb_abun)
hb_preinv <- as.matrix((hb_abuninv > 0) + 0)

hb_abuninv_5S = hb_abuninv[, colSums(hb_preinv) > 4]

## Taxonomic
hb_TD = parallel::mclapply(1:ncol(hb_abuninv),
                              function (x) try(iNEXT3D(hb_abuninv[, x],
                                                       diversity = 'TD',
                                                       q = c(0, 1, 2),
                                                       datatype = "abundance")),
                              mc.cores = 30)

names(hb_TD) = colnames(hb_abuninv)

min(sapply(hb_TD, function (x) x$TDInfo$`SC(2n)`)) # Min. coverage for n = 2n
sum(sapply(hb_TD, function (x) x$TDInfo$`SC(2n)`) > 0.95)/length(hb_TD)

saveRDS(hb_TD, "results/hb_TD.R")

### Estimate at SC = 95 %
hb_TD_est95 = parallel::mclapply(which(unname(unlist(lapply(hb_TD, function (x) x$TDInfo$`SC(2n)`) > 0.95))),
                                 function (x) try(estimate3D(hb_abuninv [,x], 
                                                             diversity = 'TD',
                                                             q = c(0), 
                                                             level = c(0.95),
                                                             datatype = "abundance")),
                                 mc.cores = 30)


names(hb_TD_est95) = colnames(hb_abuninv)[which(unname(unlist(lapply(hb_TD, function (x) x$TDInfo$`SC(2n)`) > 0.95)))]

hb_TD_est95 <- lapply(1:length(hb_TD_est95), FUN = function(x){
  hb_TD_est95[[x]]|>
    mutate(Assemblage = names(hb_TD_est95)[x])
})

hb_TD_est95 <- data.frame(data.table::rbindlist(hb_TD_est95))

write.csv(hb_TD_est95, "results/hb_TD_est95.csv")

## Phylogenetic
hb_spplist <- data.frame("species" = unique(hb_df$hbSpecies),
                         "genus" = sub( "_.*", "", unique(hb_df$hbSpecies)),
                         "family" = "Trochilidae") |>
  distinct()|>
  arrange(species)|>
  tibble::tibble()

hb_ptree <- get_tree(sp_list = hb_spplist,
                     megatrees::tree_bird_McTavish,
                     taxon = "bird",
                     scenario = "at_basal_node")

plot(hb_ptree)

hb_abuninv_phylo = hb_abuninv[hb_ptree$tip.label,]
setdiff(hb_ptree$tip.label, rownames(hb_abuninv_phylo))
setdiff(rownames(hb_abuninv_phylo), hb_ptree$tip.label)


hb_80_ptree = drop.tip(hb_ptree, hb_ptree$tip.label[-match(rownames(hb_80_abuninv_5S),
                                                           hb_ptree$tip.label)])

hb_abuninv_3S = hb_abuninv[, colSums(hb_preinv) > 2]
hb_abuninv_3S_phylo = hb_abuninv_3S[hb_ptree$tip.label,]
setdiff(hb_ptree$tip.label, rownames(hb_abuninv_3S_phylo))
setdiff(rownames(hb_abuninv_3S_phylo), hb_ptree$tip.label)

### Estimate
hb_PD_3S = parallel::mclapply(1:ncol(hb_abuninv_3S_phylo),
                                 function (x) try(iNEXT3D(hb_abuninv_3S_phylo [,x], 
                                                          diversity = 'PD',
                                                          q = c(0), 
                                                          PDtype = "meanPD",
                                                          datatype = "abundance",
                                                          PDtree = hb_ptree)),
                                 mc.cores = 30)

sapply(hb_PD_3S, class)

### Dealing with errors (do it recursively)
hb_PD_3S_errors = parallel::mclapply(which(sapply(hb_PD_3S, class) == "try-error"),
                                 function (x) try(iNEXT3D(hb_abuninv_3S_phylo [,x], 
                                                          diversity = 'PD',
                                                          q = c(0), 
                                                          PDtype = "meanPD",
                                                          datatype = "abundance",
                                                          PDtree = hb_ptree)),
                                 mc.cores = 30)

sum(sapply(hb_PD_3S_errors, class) == "try-error")

hb_PD_3S[which(sapply(hb_PD_3S, class) == "try-error")] = hb_PD_3S_errors




names(hb_PD_3S) = colnames(hb_abuninv_3S_phylo)

min(sapply(hb_PD_3S, function (x) x$PDInfo$`SC(2n)`)) # Min. coverage for n = 2n
summary(sapply(hb_PD_3S, function (x) x$PDInfo$`SC(2n)`)) # Min. coverage for n = 2n
sum(sapply(hb_PD_3S, function (x) x$PDInfo$`SC(2n)`) > 0.95)/length(hb_PD_3S)

saveRDS(hb_PD_3S, "results/hb_PD_3S.R")

### Estimate at SC = 95 %
hb_PD_3S_est95 = parallel::mclapply(which(unname(unlist(lapply(hb_PD_3S, function (x) x$PDInfo$`SC(2n)`) > 0.95))),
                                       function (x) try(estimate3D(hb_abuninv_3S_phylo [,x], 
                                                                   diversity = 'PD',
                                                                   q = c(0), 
                                                                   PDtype = "meanPD",
                                                                   datatype = "abundance",
                                                                   PDtree = hb_ptree,
                                                                   level = c(0.95))),
                                       mc.cores = 30)

### Dealing with errors (do it recursively)
hb_PD_3S_est95_errors = parallel::mclapply(which(sapply(hb_PD_3S_est95, class) == "try-error"),
                                     function (x) try(estimate3D(hb_abuninv_3S_phylo [,x], 
                                                              diversity = 'PD',
                                                              q = c(0), 
                                                              PDtype = "meanPD",
                                                              datatype = "abundance",
                                                              PDtree = hb_ptree,
                                                              level = c(0.95))),
                                     mc.cores = 30)

sum(sapply(hb_PD_3S_est95_errors, class) == "try-error")

hb_PD_3S_est95[which(sapply(hb_PD_3S_est95, class) == "try-error")] = hb_PD_3S_est95_errors

names(hb_PD_3S_est95) = colnames(hb_abuninv_3S_phylo)[which(sapply(hb_PD_3S, function (x) x$PDInfo$`SC(2n)`) > 0.95)]

hb_PD_3S_est95 <- lapply(1:length(hb_PD_3S_est95), FUN = function(x){
  hb_PD_3S_est95[[x]]|>
    mutate(Assemblage = names(hb_PD_3S_est95)[x])
})

hb_PD_3S_est95 <- data.frame(data.table::rbindlist(hb_PD_3S_est95))

write.csv(hb_PD_3S_est95, "results/hb_PD_3S_est95.csv")

## Functional
hb_spmiss = data.frame(hbSpecies = setdiff(unique(hb_df$hbSpecies), 
                          paste(traitdata::elton_birds$Genus, traitdata::elton_birds$Species, sep = "_")))

lapply(gsub(".*_", "", hb_spmiss$hbSpecies), 
       function (x) paste(traitdata::elton_birds$Genus, traitdata::elton_birds$Species, sep = "_")[which(elton_birds$Species == x)])

hb_spmiss$hbSpecies_elton = c("Aglaiocercus_kingi", "Amazilia_amabilis", 
                                  "Amazilia_boucardi", "Amazilia_chionogaster", 
                                  "Amazilia_fimbriata", "Amazilia_franciae",
                                  "Amazilia_saucerrottei", "Amazilia_versicolor", 
                                  "Amazilia_violiceps", "Amazilia_amazilia",
                                  "Anthocephala_floriceps", "Anthracothorax_dominicus",
                                  "Hylocharis_leucotis", "Hylocharis_xantusii",
                                  "Calliphlox_bryantae", "Calliphlox_mitchellii",
                                  "Campylopterus_largipennis", "Campylopterus_largipennis",
                                  "Amazilia_lactea", "Amazilia_candida",
                                  "Hylocharis_cyanus", "Hylocharis_eliciae",
                                  "Damophila_julie",  "Amazilia_brevirostris",
                                  "Lepidopyga_coeruleogularis", "Lepidopyga_goudoti",
                                  "Hylocharis_grayi", "Hylocharis_humboldtii",
                                  "Amazilia_leucogaster", "Lepidopyga_lilliae",
                                  "Clytolaema_rubricauda", "Coeligena_torquata",
                                  "Coeligena_bonapartei", "Coeligena_bonapartei",
                                  "Coeligena_torquata", "Colibri_thalassinus",
                                  "Chlorostilbon_auriceps", "Chlorostilbon_canivetii",
                                  "Cynanthus_latirostris", "Chlorostilbon_forficatus",
                                  "Cynanthus_latirostris", "Amazilia_viridicauda",
                                  "Eriocnemis_alinae", "Eugenes_fulgens",
                                  "Aphantochroa_cirrochloris", "Thalurania_ridgwayi",
                                  "Goethalsia_bella", "Heliangelus_amethysticollis",
                                  "Heliangelus_amethysticollis", "Amazilia_viridifrons",
                                  "Lophornis_chalybeus", "Elvira_chionura",
                                  "Elvira_cupreiceps", "Calliphlox_evelynae",
                                  "Calliphlox_evelynae", "Ocreatus_underwoodii",
                                  "Ocreatus_underwoodii", "Oreotrochilus_melanogaster",
                                  "Oreotrochilus_estella", "Oxypogon_guerinii",
                                  "Oxypogon_guerinii", "Oxypogon_guerinii",
                                  "Campylopterus_curvipennis", "Campylopterus_rufus",
                                  "Cynanthus_sordidus", NA,
                                  "Phaethornis_longirostris", "Amazilia_amabilis",
                                  "Amazilia_rosenbergi", "Cyanophaia_bicolor",
                                  "Chlorostilbon_maugaeus", "Chlorostilbon_ricordii",
                                  "Chlorostilbon_swainsonii", "Sappho_sparganura",
                                  "Amazilia_beryllina", "Amazilia_castaneiventris",
                                  "Amazilia_cyanifrons", "Amazilia_cyanocephala",
                                  "Amazilia_cyanura", "Amazilia_edward",
                                  "Amazilia_saucerrottei", "Amazilia_tobaci",
                                  "Amazilia_viridigaster", "Schistes_geoffroyi",
                                  "Stellula_calliope", "Atthis_ellioti",
                                  "Atthis_heloisa", "Stephanoxis_lalandi",
                                  "Leucippus_chlorocercus", "Leucippus_baeri", 
                                  "Leucippus_taczanowskii", "Urochroa_bougueri")


hb_spmiss = rbind(hb_spmiss, data.frame(hbSpecies = intersect(unique(hb_df$hbSpecies), 
                                                              paste(traitdata::elton_birds$Genus, traitdata::elton_birds$Species, sep = "_")),
                                        hbSpecies_elton = intersect(unique(hb_df$hbSpecies), 
                                                              paste(traitdata::elton_birds$Genus, traitdata::elton_birds$Species, sep = "_"))))

hb_spmiss$notes = NA
hb_spmiss$notes[hb_spmiss$hbSpecies %in% c("Anthocephala_berlepschi",
                                                       "Anthracothorax_aurulentus",
                                                       "Campylopterus_calcirupicola",
                                                       "Campylopterus_diamantinensis",
                                                       "Coeligena_conradii",
                                                       "Coeligena_consita",
                                                       "Coeligena_eos",
                                                       "Coeligena_inca",
                                                       "Colibri_cyanotus",
                                                       "Cynanthus_doubledayi",
                                                       "Cynanthus_lawrencei",
                                                       "Eugenes_spectabilis",
                                                       "Heliangelus_clarisse",
                                                       "Heliangelus_spencei",
                                                       "Lophornis_verreauxii",
                                                       "Nesophlox_lyrura",
                                                       "Ocreatus_addae",
                                                       "Ocreatus_peruanus",
                                                       "Oreotrochilus_stolzmanni",
                                                       "Oxypogon_cyanolaemus",
                                                       "Oxypogon_lindenii",
                                                       "Oxypogon_stuebelii",
                                                       "Phaethornis_mexicanus",
                                                       "Polyerata_decora",
                                                       "Saucerottia_hoffmanni",
                                                       "Schistes_albogularis",
                                                       "Stephanoxis_loddigesii",
                                                       "Urochroa_leucura")] = "Has since been splitted into a different species"

hb_spmiss$notes[hb_spmiss$hbSpecies %in% c("Oreotrochilus_cyanolaemus")] = "Newly discovered. A related species was used instead"

hb_spmiss$notes[hb_spmiss$hbSpecies %in% c("Phaethornis_aethopygus")] = "Thought to be a hybrid Phaethornis ruber x rupurumii before 2009"


## Getting df of traits
hb_traits_df = traitdata::elton_birds|>
  mutate(hbSpecies_elton = paste(Genus, Species, sep = "_"))|>
  left_join(hb_spmiss, by = "hbSpecies_elton")|>
  mutate(hbSpecies = ifelse(hbSpecies_elton %in% c("Phaethornis_ruber", "Phaethornis_rupurumii"),
                            "Phaethornis_aethopygus", hbSpecies))|>
  filter(!is.na(hbSpecies))|>

  dplyr::select(hbSpecies, contains(c("Diet", "ForStrat",
                                                       "Nocturnal", "BodyMass")))|>
  dplyr::select(!contains(c("Source", "EnteredBy", "Comment")))|>
  mutate(across(-all_of(c("hbSpecies", "Diet.5Cat", "Diet.Certainty")),
    ~ as.numeric(as.character(.x))))|>
  group_by(hbSpecies)|>
  summarise_all(funs(if(is.numeric(.)) mean(., na.rm = TRUE) else first(.)))|>
  as.data.frame()

hb_traits_df = rbind(hb_traits_df, data.frame("hbSpecies" = c("Phaethornis_rupurumii", "Phaethornis_ruber"),
                               traitdata::elton_birds[grep("Phaethornis ruber|Phaethornis rupurumii", traitdata::elton_birds$scientificNameStd),
                                                      colnames(hb_traits_df)[-1]])) # at first I checked adding species to make sure that the order was right
### Functional distances
hb_traits_sel = hb_traits_df[, which(colnames(hb_traits_df) %in% setdiff(c("hbSpecies",
                                                   names(which(colSums(hb_traits_df[which(lapply(1:ncol(hb_traits_df),
                                                                                                 function(x) class(hb_traits_df[, x])) == "numeric")]) > 0))),
                                                   c("ForStrat.SpecLevel", "BodyMass.SpecLevel")))]|>
  remove_rownames()|>
  column_to_rownames(var = 'hbSpecies')


hb_fdist <- gawdis(hb_traits_sel, w.type = "optimized") #Get functional distances
hb_fdist <- as.matrix(hb_fdist)
hb_fdist = hb_fdist[rownames(hb_abuninv), rownames(hb_abuninv)]

hb_fdist_3S = hb_fdist[rownames(hb_abuninv_3S_phylo), rownames(hb_abuninv_3S_phylo)]

hb_fdist_5S = hb_fdist[rownames(hb_80_abuninv_5S), rownames(hb_80_abuninv_5S)]
  
### FD
hb_FD_3S = parallel::mclapply(1:ncol(hb_abuninv_3S_phylo),
                                 function (x) try(iNEXT3D(hb_abuninv_3S_phylo[, x],
                                                          diversity = 'FD',
                                                          q = c(0),
                                                          datatype = "abundance", 
                                                          nboot = 20, 
                                                          FDdistM = hb_fdist_3S, 
                                                          FDtype = "AUC")),
                                 mc.cores = 30)

sapply(hb_FD_3S, class)
names(hb_FD_3S) = colnames(hb_abuninv_3S_phylo)

min(sapply(hb_FD_3S, function (x) x$FDInfo$`SC(2n)`)) # Min. coverage for n = 2n
sum(sapply(hb_FD_3S, function (x) x$FDInfo$`SC(2n)`) > 0.95)/length(hb_FD_3S)

saveRDS(hb_FD_3S, "results/hb_FD_3S.R")

### Estimate at SC = 95 %
hb_FD_3S_est95 = parallel::mclapply(which(unname(unlist(lapply(hb_FD_3S, function (x) x$FDInfo$`SC(2n)`) > 0.95))),
                                       function (x) try(estimate3D(hb_abuninv_3S_phylo [,x], 
                                                                   diversity = 'FD',
                                                                   q = c(0), 
                                                                   level = c(0.95),
                                                                   datatype = "abundance", 
                                                                   nboot = 20, 
                                                                   FDdistM = hb_fdist_3S, 
                                                                   FDtype = "AUC")),
                                       mc.cores = 30)


names(hb_FD_3S_est95) = colnames(hb_abuninv_3S_phylo)[which(unname(unlist(lapply(hb_FD_3S, function (x) x$FDInfo$`SC(2n)`) > 0.95)))]

hb_FD_3S_est95 <- lapply(1:length(hb_FD_3S_est95), FUN = function(x){
  hb_FD_3S_est95[[x]]|>
    mutate(Assemblage = names(hb_FD_3S_est95)[x])
})

hb_FD_3S_est95 <- data.frame(data.table::rbindlist(hb_FD_3S_est95))

write.csv(hb_FD_3S_est95, "results/hb_FD_3S_est95.csv")

# Networks ----
hb_80_webs = lapply(sort(unique(hb_80_df$nth_cell)), 
                    function(x) table(hb_80_df$hbSpecies[hb_80_df$nth_cell == x],
                                      hb_80_df$speciesPlant[hb_80_df$nth_cell == x]))

names(hb_80_webs) = sort(unique(hb_80_df$nth_cell))
  
netst_metrics = c("connectance", "cluster coefficient", "compartment diversity", "NODF", "H2")

hb_80_netst = data.table::rbindlist(lapply(hb_80_webs, 
                                              function (x) data.frame("Assemblage" = names(hb_80_webs)[x], t(bipartite::networklevel(x, index = netst_metrics)))))


hb_80_webs_5S = hb_80_webs[lapply(hb_80_webs, nrow) > 4]
hb_80_netst_5S = data.table::rbindlist(lapply(1:length(hb_80_webs_5S), 
                                           function (x) data.frame("Assemblage" = names(hb_80_webs_5S)[x], t(bipartite::networklevel(hb_80_webs_5S[[x]], index = netst_metrics)))))



bipartite::plotweb(hb_80_webs_5S[[1]])

# Put everything together ----
hb_80_allvars_5S = hb_80_df[, c("nth_cell", "cell_long", "cell_lat")]|>
  group_by(nth_cell, cell_long, cell_lat)|>
  summarise(.groups = "drop")|>
  rename(Assemblage = nth_cell)|>
  mutate(Assemblage = as.character(Assemblage))|>
  right_join(hb_80_TD_5S_est50[, c("Assemblage", "qTD")], by = "Assemblage")|>
  left_join(hb_80_PD_5S_est50[, c("Assemblage", "qPD")], by = "Assemblage")|>
  left_join(hb_80_FD_5S_est50[, c("Assemblage", "qFD")], by = "Assemblage")|>
  left_join(hb_80_netst_5S, by = "Assemblage")
  #pivot_longer(cols = qTD:cluster.coefficient.LL,
  #             names_to = "variable",
  #             values_to = "value")

write.csv(hb_80_allvars_5S, "results/hb_80_allvars_5S.csv")

# Maps ----
## Diversity
hb_raster_qTD_est95 <- terra::setValues(selelev_mean_EEA, NA)
hb_raster_qTD_est95[as.numeric(hb_TD_est95$Assemblage)] = hb_TD_est95$qTD
writeRaster(hb_raster_qTD_est95, "results/hb_raster_qTD_est95.tiff")

hb_raster_qTD_est95_wgs84 = terra::project(hb_raster_qTD_est95, CRS("EPSG:4326"))
writeRaster(hb_raster_qTD_est95_wgs84, "results/hb_raster_qTD_est95_wgs84.tiff")

hb_raster_qPD_3S_est95 <- terra::setValues(selelev_mean_EEA, NA)
hb_raster_qPD_3S_est95[as.numeric(hb_PD_3S_est95$Assemblage)] = hb_PD_3S_est95$qPD
writeRaster(hb_raster_qPD_3S_est95, "results/hb_raster_qPD_3S_est95.tiff")

hb_raster_qPD_3S_est95_wgs84 = terra::project(hb_raster_qPD_3S_est95, CRS("EPSG:4326"))
writeRaster(hb_raster_qPD_3S_est95_wgs84, "results/hb_raster_qPD_3S_est95_wgs84.tiff")

hb_raster_qFD_3S_est95 <- terra::setValues(selelev_mean_EEA, NA)
hb_raster_qFD_3S_est95[as.numeric(hb_FD_3S_est95$Assemblage)] = hb_FD_3S_est95$qFD
writeRaster(hb_raster_qFD_3S_est95, "results/hb_raster_qFD_3S_est95.tiff")

hb_raster_qFD_3S_est95_wgs84 = terra::project(hb_raster_qFD_3S_est95, CRS("EPSG:4326"))
writeRaster(hb_raster_qFD_3S_est95_wgs84, "results/hb_raster_qFD_3S_est95_wgs84.tiff")

## Network metrics
hb_80_raster_connectance_5S <- terra::setValues(selelev_mean_EEA, NA)
hb_80_raster_connectance_5S[as.numeric(hb_80_allvars_5S$Assemblage)] = hb_80_allvars_5S$connectance
writeRaster(hb_80_raster_connectance_5S, "results/hb_80_raster_connectance_5S.tiff")

hb_80_raster_connectance_5S_wgs84 = terra::project(hb_80_raster_connectance_5S, CRS("EPSG:4326"))
writeRaster(hb_80_raster_connectance_5S_wgs84, "results/hb_80_raster_connectance_5S_wgs84.tiff")

hb_80_raster_compdiv_5S <- terra::setValues(seltemp_mean_EEA, NA)
hb_80_raster_compdiv_5S[as.numeric(hb_80_allvars_5S$Assemblage)] = hb_80_allvars_5S$compartment.diversity
writeRaster(hb_80_raster_compdiv_5S, "results/hb_80_raster_compdiv_5S.tiff")

hb_80_raster_compdiv_5S_wgs84 = terra::project(hb_80_raster_compdiv_5S, CRS("EPSG:4326"))
writeRaster(hb_80_raster_compdiv_5S_wgs84, "results/hb_80_raster_compdiv_5S_wgs84.tiff")

hb_80_raster_cc_5S <- terra::setValues(seltemp_mean_EEA, NA)
hb_80_raster_cc_5S[as.numeric(hb_80_allvars_5S$Assemblage)] = hb_80_allvars_5S$cluster.coefficient
writeRaster(hb_80_raster_cc_5S, "results/hb_80_raster_cc_5S.tiff")

hb_80_raster_cc_5S_wgs84 = terra::project(hb_80_raster_cc_5S, CRS("EPSG:4326"))
writeRaster(hb_80_raster_cc_5S_wgs84, "results/hb_80_raster_cc_5S_wgs84.tiff")

hb_80_raster_NODF_5S <- terra::setValues(seltemp_mean_EEA, NA)
hb_80_raster_NODF_5S[as.numeric(hb_80_allvars_5S$Assemblage)] = hb_80_allvars_5S$NODF
writeRaster(hb_80_raster_NODF_5S, "results/hb_80_raster_NODF_5S.tiff")

hb_80_raster_NODF_5S_wgs84 = terra::project(hb_80_raster_NODF_5S, CRS("EPSG:4326"))
writeRaster(hb_80_raster_NODF_5S_wgs84, "results/hb_80_raster_NODF_5S_wgs84.tiff")

hb_80_raster_H2_5S <- terra::setValues(seltemp_mean_EEA, NA)
hb_80_raster_H2_5S[as.numeric(hb_80_allvars_5S$Assemblage)] = hb_80_allvars_5S$H2
writeRaster(hb_80_raster_H2_5S, "results/hb_80_raster_H2_5S.tiff")

hb_80_raster_H2_5S_wgs84 = terra::project(hb_80_raster_H2_5S, CRS("EPSG:4326"))
writeRaster(hb_80_raster_H2_5S_wgs84, "results/hb_80_raster_H2_5S_wgs84.tiff")

hb_80_raster_cchb_5S <- terra::setValues(seltemp_mean_EEA, NA)
hb_80_raster_cchb_5S[as.numeric(hb_80_allvars_5S$Assemblage)] = hb_80_allvars_5S$cluster.coefficient.HL
writeRaster(hb_80_raster_cchb_5S, "results/hb_80_raster_cchb_5S.tiff")

hb_80_raster_cchb_5S_wgs84 = terra::project(hb_80_raster_cchb_5S, CRS("EPSG:4326"))
writeRaster(hb_80_raster_cchb_5S_wgs84, "results/hb_80_raster_cchb_5S_wgs84.tiff")

hb_80_raster_ccpl_5S <- terra::setValues(seltemp_mean_EEA, NA)
hb_80_raster_ccpl_5S[as.numeric(hb_80_allvars_5S$Assemblage)] = hb_80_allvars_5S$cluster.coefficient.LL
writeRaster(hb_80_raster_ccpl_5S, "results/hb_80_raster_ccpl_5S.tiff")

hb_80_raster_ccpl_5S_wgs84 = terra::project(hb_80_raster_ccpl_5S, CRS("EPSG:4326"))
writeRaster(hb_80_raster_ccpl_5S_wgs84, "results/hb_80_raster_ccpl_5S_wgs84.tiff")

## Together
hb_80_raster_5S <- rast(list(hb_80_raster_qTD_5S, hb_80_raster_qPD_5S,
                             hb_80_raster_qFD_5S, hb_80_raster_connectance_5S,
                             hb_80_raster_compdiv_5S, hb_80_raster_cc_5S,
                             hb_80_raster_NODF_5S, hb_80_raster_H2_5S,
                             hb_80_raster_cchb_5S, hb_80_raster_ccpl_5S))
names(hb_80_raster_5S) <- c("TD","PD","FD", "Connectance", "Compartment diversity",
                            "Cluster coefficient", "NODF", "H2", "Cluster coefficient - HB",
                            "Cluster coefficient - Plant")


ggplot() + 
  geom_spatvector(data = whem_sf_EEA) +
  geom_raster(data = hb_80_allvars_5S[hb_80_allvars_5S$variable == "qTD",],
              aes(x = cell_long, y = cell_lat, fill = value))

# Environmental variables ----
download.file(
  "https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_2.5m_bio.zip",
  destfile = "results/wc2.1_2.5m_bio.zip"
)
unzip("results/wc2.1_2.5m_bio.zip", exdir = "results/worldclim_bio")

temp_raw   <- rast("results/worldclim_bio/wc2.1_2.5m_bio_1.tif")   # Annual Mean Temp
precip_raw <- rast("results/worldclim_bio/wc2.1_2.5m_bio_12.tif")  # Annual Precip

temp_aligned   <- resample(crop(temp_raw, ext(selelev_mean)), selelev_mean, method = "bilinear")
precip_aligned <- resample(crop(precip_raw, ext(selelev_mean)), selelev_mean, method = "bilinear")

writeRaster(temp_aligned,   "results/env_temp_mean_wgs84.tiff",   overwrite = TRUE)
writeRaster(precip_aligned, "results/env_precip_mean_wgs84.tiff", overwrite = TRUE)

# SPEIbase v2.9/2.10, 12-month SPEI, monthly global 0.5° grid, 1901-2024
download.file(
  "https://spei.csic.es/spei_database/nc/spei12.nc",
  destfile = "results/spei12.nc"
)

spei_stack <- rast("results/spei12.nc")

# Layer names are dates -- subset to Jan 2015 - Dec 2024, then average
dates <- time(spei_stack)
idx <- which(dates >= as.Date("2015-01-01") & dates <= as.Date("2024-12-31"))
spei_mean_raw <- mean(spei_stack[[idx]], na.rm = TRUE)

spei_aligned <- resample(crop(spei_mean_raw, ext(selelev_mean)), selelev_mean, method = "bilinear")
writeRaster(spei_aligned, "results/env_spei_mean_wgs84.tiff", overwrite = TRUE)

# Population 
download.file(
  "https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_POP_GLOBE_R2023A/GHS_POP_E2020_GLOBE_R2023A_4326_30ss/V1-0/GHS_POP_E2020_GLOBE_R2023A_4326_30ss_V1_0.zip",
  destfile = "results/GHS_pop_2020.zip"
)
unzip("results/GHS_pop_2020.zip", exdir = "results/ghs_population2020")

pop_raw   <- rast("results/ghs_population2020/GHS_POP_E2020_GLOBE_R2023A_4326_30ss_V1_0.tif")  


pop_aligned   <- resample(crop(pop_raw, ext(selelev_mean)), selelev_mean, method = "bilinear")
plot(pop_aligned)
writeRaster(pop_aligned,   "results/pop_2020_wgs84.tiff",   overwrite = TRUE)

# NPP
# List all four tiles and mosaic them into one raster
npp_tiles <- list.files(path = "results/",
  pattern = "^env_npp_mean-.*\\.tif$", 
  full.names = TRUE)
npp_tiles 

npp_list <- lapply(npp_tiles, rast)
npp_raw <- do.call(mosaic, npp_list)  # or use merge() if mosaic() errors -- see note below

plot(npp_raw)  # quick visual check that it looks like one continuous Americas raster, no gaps/seams

npp_aligned <- resample(crop(npp_raw, ext(selelev_mean)), selelev_mean, method = "bilinear")
writeRaster(npp_aligned, "results/env_npp_mean_wgs84.tiff", overwrite = TRUE)

# Fire
fire_raw = rast("results/env_fire_occurrence.tif")
plot(fire_raw) 

fire_aligned <- resample(crop(fire_raw, ext(selelev_mean)), selelev_mean, method = "bilinear")
writeRaster(fire_aligned, "results/env_fire_mean_wgs84.tiff", overwrite = TRUE)
