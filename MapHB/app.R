#
# Shiny web application — Hummingbird biodiversity maps (Americas)
#
# Find out more about building applications with Shiny here:
#    https://shiny.posit.co/
#

library(shiny)
library(shinydashboard)
library(leaflet)
library(raster)
library(terra)
library(ggplot2)
library(ape)             # phylogenetic tree handling / plotting
library(bipartite)       # hummingbird-plant interaction network plot
library(aws.s3)          # S3-compatible client, used here against Cloudflare R2
library(arrow)           # reads the Parquet version of the photo/interaction data

# ---- R2 (Cloudflare) connection --------------------------------------------
# Credentials are read from environment variables:
#   - Locally: set in ~/.Renviron (never committed to git)
#   - On Posit Connect: set as environment variables on the deployed app
r2_key    <- Sys.getenv("AWS_ACCESS_KEY_ID")
r2_secret <- Sys.getenv("AWS_SECRET_ACCESS_KEY")
r2_base   <- "91a6e033360b81fa6d1f2bd942507685.r2.cloudflarestorage.com"
r2_bucket <- "pri-2026-hb-data"

if (r2_key == "" || r2_secret == "") {
  stop("AWS_ACCESS_KEY_ID / AWS_SECRET_ACCESS_KEY are not set. ",
       "Set them in ~/.Renviron locally, or as environment variables on Posit Connect.")
}

# ---- Paths ------------------------------------------------------------
# app.R lives inside the MapHB subfolder, one level below the project
# root (which contains `results/`). Rather than assuming the working
# directory is always MapHB (true for "Run App" and Posit Connect, but
# NOT true if someone runs source("MapHB/app.R") from one level up), we
# search upward from the current working directory for a folder that
# actually contains "results/hb_ptree.txt" -- a small file we know must
# exist at the project root. This makes path resolution independent of
# how/from-where the app happens to be launched.
find_project_root <- function(marker = file.path("results", "hb_ptree.txt"),
                              max_levels = 6) {
  dir <- normalizePath(getwd(), mustWork = TRUE)
  for (i in seq_len(max_levels)) {
    if (file.exists(file.path(dir, marker))) return(dir)
    parent <- dirname(dir)
    if (parent == dir) break  # reached filesystem root, stop
    dir <- parent
  }
  stop("Could not locate project root (folder containing '", marker,
       "') by searching upward from ", getwd(), ". ",
       "Check that the app is being run from within the project folder tree.")
}

project_root <- find_project_root()
message("app.R working dir: ", getwd())
message("Resolved project root: ", project_root)

# Local cache directory for files pulled from R2. Downloads only happen once
# per app instance -- subsequent app restarts on the same machine/container
# reuse the cached copy instead of re-downloading from R2 every time.
cache_dir <- file.path(project_root, "results_cache")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
message("Cache dir resolves to: ", cache_dir)

# Downloads `object_key` from the R2 bucket into cache_dir, unless it's
# already been downloaded, and returns the local path to use with rast()/
# read_parquet()/etc.
#
# NOTE: uses get_object() + writeBin() rather than aws.s3::save_object().
# save_object() routes its response through internal parsing helpers that
# can misfire on larger binary payloads (observed failing specifically on
# a ~7.9 MB file while smaller ~450 KB files succeeded, with an unrelated
# -looking error thrown from path_to_connection()). get_object() returns
# the raw bytes directly with no such parsing step, so we write them to
# disk ourselves -- simpler and size-independent.
fetch_from_r2 <- function(object_key) {
  local_path <- file.path(cache_dir, object_key)
  if (!file.exists(local_path)) {
    message("Downloading from R2: ", object_key)
    raw_bytes <- get_object(
      object   = object_key,
      bucket   = r2_bucket,
      key      = r2_key,
      secret   = r2_secret,
      region   = "",
      base_url = r2_base,
      use_https = TRUE,
      check_region = FALSE
    )
    writeBin(raw_bytes, local_path)
  }
  local_path
}

# ---- Extent (Americas) -------------------------------------------------
americas <- ext(-170, -30, -60, 85)

# ---- Load & crop rasters (now pulled from R2, cached locally) -------------
# Loads one raster from R2 with clear diagnostics if anything goes wrong,
# instead of a bare rast()/crop() call whose error message doesn't say
# which file was responsible.
load_raster <- function(object_key) {
  local_path <- fetch_from_r2(object_key)
  if (!file.exists(local_path)) {
    stop("Downloaded file is missing on disk: ", local_path)
  }
  sz <- file.info(local_path)$size
  message("  -> ", object_key, " (", sz, " bytes) at ", local_path)
  r <- tryCatch(
    rast(local_path),
    error = function(e) {
      stop("Failed to open '", object_key, "' at '", local_path,
           "' with terra::rast(): ", conditionMessage(e))
    }
  )
  crop(r, americas)
}

TD           <- load_raster("hb_raster_qTD_est95_wgs84.tiff")
PD           <- load_raster("hb_raster_qPD_3S_est95_wgs84.tiff")
FD           <- load_raster("hb_raster_qFD_3S_est95_wgs84.tiff")
connectance  <- load_raster("hb_80_raster_connectance_5S_wgs84.tiff")
compdiv      <- load_raster("hb_80_raster_compdiv_5S_wgs84.tiff")
cc           <- load_raster("hb_80_raster_cc_5S_wgs84.tiff")
NODF         <- load_raster("hb_80_raster_NODF_5S_wgs84.tiff")
H2           <- load_raster("hb_80_raster_H2_5S_wgs84.tiff")
cchb         <- load_raster("hb_80_raster_cchb_5S_wgs84.tiff")
ccpl         <- load_raster("hb_80_raster_ccpl_5S_wgs84.tiff")
selelev_mean <- load_raster("selelev_mean.tiff")

message("All 11 rasters loaded and cropped successfully.")

# ---- Environmental layers (2015-2024 means/totals, replacing the old live
# FIRMS/NDVI/SPEI WMS tiles with static pre-computed rasters) --------------
env_fire   <- load_raster("env_fire_mean_wgs84.tiff")
env_npp    <- load_raster("env_npp_mean_wgs84.tiff")
env_spei   <- load_raster("env_spei_mean_wgs84.tiff")
env_temp   <- load_raster("env_temp_mean_wgs84.tiff")
env_precip <- load_raster("env_precip_mean_wgs84.tiff")
env_pop    <- load_raster("pop_2020_wgs84.tiff")

message("All 6 environmental layers loaded and cropped successfully.")

env_rasters <- list(
  fire = env_fire, npp = env_npp, spei = env_spei,
  temp = env_temp, precip = env_precip, pop = env_pop
)

env_labels <- c(
  fire   = "Fire occurrence (2015-2024)",
  npp    = "NPP mean (2015-2024)",
  spei   = "Drought index (SPEI, 2015-2024 mean)",
  temp   = "Mean Annual Temperature",
  precip = "Mean Annual Precipitation",
  pop    = "Population density (2020)"
)

env_palette_colors <- list(
  fire   = c("#FFFFFF", "#B22222"),               # white -> firebrick red
  npp    = c("#FFFFFF", "#1B7837"),               # white -> dark green
  spei   = c("#A0522D", "#FFFFFF", "#1E90FF"),    # brown (dry) -> white -> blue (wet)
  temp   = c("#2166AC", "#FFFFFF", "#B2182B"),    # blue (cold) -> white -> red (hot)
  precip = c("#FFFFFF", "#08519C"),               # white -> dark blue
  pop    = c("#FFFFFF", "#4B0082")                # white -> indigo
)

# Precompute each layer's non-NA cell values once; reused both as the
# palette domain below and as this layer's legend values later, so we
# don't re-extract potentially large rasters' full value vectors every
# time a layer is toggled on/off in a session.
env_values <- lapply(env_rasters, function(r) values(r, na.rm = TRUE))

env_palettes <- lapply(names(env_rasters), function(v) {
  colorNumeric(palette = env_palette_colors[[v]], domain = env_values[[v]], na.color = "transparent")
})
names(env_palettes) <- names(env_rasters)

# ---- Photo / interaction data ---------------------------------------------
# Now stored as Parquet on R2 (smaller + faster than the raw CSV) instead of
# results/hb_df.csv. Still cached to .rds after the first parse per app
# instance, same as before -- delete hb_photos_cache.rds if the source
# Parquet file on R2 is ever updated.
photos_cache_path <- file.path(project_root, "results", "hb_photos_cache.rds")
dir.create(dirname(photos_cache_path), showWarnings = FALSE, recursive = TRUE)

if (file.exists(photos_cache_path)) {
  cached <- readRDS(photos_cache_path)
  photos        <- cached$photos
  photos_unique <- cached$photos_unique
} else {
  photos <- as.data.frame(arrow::read_parquet(fetch_from_r2("hb_df.parquet")))
  
  # Assign each photo row to a cell using the SAME grid as the biodiversity
  # rasters (rather than the dataframe's own `nth_cell`, which was built in a
  # different, unknown projected CRS). This guarantees the map click and the
  # photo data always refer to the same cell.
  photos$cell_id <- terra::cellFromXY(TD, cbind(photos$decimalLongitude, photos$decimalLatitude))
  
  # Pull the correct photo URL out of the (possibly multi-URL) `identifier`
  # field, using the trailing "_N" in imageFileName as the index into it.
  get_photo_url <- function(identifier, imageFileName) {
    urls <- trimws(strsplit(identifier, ";")[[1]])
    idx  <- suppressWarnings(as.integer(sub(".*_(\\d+)\\.[a-zA-Z]+$", "\\1", imageFileName)))
    if (is.na(idx) || idx < 1 || idx > length(urls)) urls[1] else urls[idx]
  }
  photos$photo_url <- mapply(get_photo_url, photos$identifier, photos$imageFileName)
  
  # One marker per unique photo (a photo can appear multiple times if several
  # candidate plant IDs were scored against it).
  photos_unique <- photos[!duplicated(photos$imageFileName), ]
  
  saveRDS(list(photos = photos, photos_unique = photos_unique), photos_cache_path)
}

# ---- Phylogenetic tree (hummingbirds) -------------------------------------
# Small text file -- left in git as-is, not moved to R2.
hb_tree <- ape::read.tree(file.path(project_root, "results", "hb_ptree.txt"))

# ---- Raster registry (name -> raster object) -----------------------------
raster_list <- list(
  TD = TD, PD = PD, FD = FD,
  connectance = connectance, NODF = NODF, H2 = H2
)

# Variables shown in the click-popup summary table, in this order.
# Filtered against names(raster_list) so that trimming/renaming the raster
# registry above can never again point the click handler at a raster that
# doesn't exist (this was the cause of the "x = NULL" extract() error).
popup_vars <- intersect(c("TD", "PD", "FD", "connectance", "NODF", "H2", "cc"),
                        names(raster_list))

popup_labels <- c(
  TD = "Taxonomic Diversity", PD = "Phylogenetic Diversity", FD = "Functional Diversity",
  connectance = "Connectance", NODF = "NODF", H2 = "H2", cc = "Cluster Coefficient"
)

# ---- Color palettes (edit hex codes / add entries as needed) -------------
palette_colors <- list(
  TD          = c("#FFFFFF", "#0073D1"),  # white -> blue
  PD          = c("#FFFFFF", "#CC002B"),  # white -> coral red
  FD          = c("#FFFFFF", "#E89E00"),  # white -> yellow
  connectance = c("#FFFFFF", "#009180"),  # white -> turquoise green
  NODF        = c("#FFFFFF", "#B51AFF"),  # white -> purple
  H2          = c("#FFFFFF", "#FF571F")   # white -> orange
)

# Precompute each raster's non-NA cell values once; reused both as the
# palette domain below and as the biological legend's values later, so we
# don't re-extract the full value vector every time the variable changes.
raster_values <- lapply(raster_list, function(r) values(r, na.rm = TRUE))

# Precompute one colorNumeric palette function per variable (built once,
# from the actual cell values, not from the raster object itself).
palettes <- lapply(names(raster_list), function(v) {
  colorNumeric(palette = palette_colors[[v]], domain = raster_values[[v]], na.color = "transparent")
})
names(palettes) <- names(raster_list)

variable_choices <- list(
  "Taxonomic Diversity (SC = 95 %)"                    = "TD",
  "Phylogenetic Diversity (SC = 95 %)"                 = "PD",
  "Functional Diversity (SC = 95 %)"                   = "FD",
  "Network: connectance"                               = "connectance",
  "Network: NODF"                                      = "NODF",
  "Network: H2"                                        = "H2"
)

# ---- UI -------------------------------------------------------------------
ui <- dashboardPage(
  dashboardHeader(disable = TRUE),
  dashboardSidebar(
    radioButtons("variable", "Select Variable", choices = variable_choices)
  ),
  dashboardBody(
    leafletOutput(outputId = "map", height = 600),
    br(),
    fluidRow(
      column(6,
             h4("Values at last clicked pixel"),
             tableOutput("clickTable")
      ),
      column(6,
             h4("Hummingbird subtree (this cell)"),
             plotOutput("phyloPlot", height = 400)
      )
    ),
    br(),
    # Full-width row: the network plot has many rotated species-name labels
    # that need more horizontal room than a 1/3-width column can offer.
    # Height bumped well up from the other panels -- text size in R plots
    # is roughly fixed in absolute (point/pixel) terms, so a taller canvas
    # means the same label occupies a smaller fraction of the [0,1] plot
    # space, giving edge labels more room without shrinking the font.
    fluidRow(
      column(12,
             h4("Hummingbird - plant network (this cell)"),
             plotOutput("networkPlot", height = 650)
      )
    )
  ),
  title = "Interactive Maps"
)

# ---- Server -----------------------------------------------------------
server <- function(input, output, session) {
  
  selectedRaster <- reactive({
    req(input$variable)
    raster_list[[input$variable]]
  })
  
  selectedPal <- reactive({
    req(input$variable)
    palettes[[input$variable]]
  })
  
  # Stores the values of ALL variables at the last clicked point
  clickValues <- reactiveVal(NULL)
  
  # Stores the photo rows belonging to the last clicked cell
  clickCellPhotos <- reactiveVal(NULL)
  
  # Tracks whether the (expensive) photo layer has been added yet, so we
  # only build it once per session -- the first time the user turns it on.
  photosLoaded <- reactiveVal(FALSE)
  
  # Tracks, per environmental layer, whether it has been added to the map
  # yet. Like the photo layer above, each env raster is only rendered
  # (including its costly project = TRUE reprojection) the first time its
  # checkbox is actually turned on, instead of all 6 being rendered
  # unconditionally at app startup. This keeps the first renderLeaflet()
  # call cheap regardless of how many environmental layers exist.
  envLoaded <- reactiveValues(fire = FALSE, npp = FALSE, spei = FALSE,
                              temp = FALSE, precip = FALSE, pop = FALSE)
  
  # When a photo marker is clicked, the underlying click event also reaches
  # the generic map click handler below. This flag lets the marker-click
  # handler tell the map-click handler "this click was already handled by
  # a marker -- skip the raster popup," so clicking a photo shows ONLY the
  # photo thumbnail/info, never the raster (TD/PD/FD/network) popup.
  suppressMapClick <- reactiveVal(FALSE)
  
  observeEvent(input$map_groups, {
    # ignoreInit = TRUE: only respond to the user actually checking a box.
    # Without this, this observer would also fire once during the initial
    # reactive flush using whatever value input$map_groups happens to hold
    # at that moment -- which is exactly the eager-loading behavior we're
    # trying to avoid.
    if ("Photos" %in% input$map_groups && !photosLoaded()) {
      leafletProxy("map") |>
        addCircleMarkers(
          data = photos_unique,
          lng = ~decimalLongitude, lat = ~decimalLatitude,
          # Bumped from radius = 4 / stroke = FALSE: a 4px fill-only dot was
          # a very small, low-contrast click target, making it easy to miss
          # the marker and hit the underlying map instead (which shows the
          # raster-value popup rather than the photo popup). A larger radius
          # plus a thin stroke gives a bigger, more visible, easier-to-hit
          # target without changing the map's overall visual weight.
          radius = 7, stroke = TRUE, weight = 1, color = "#7A0019",
          fillOpacity = 0.75, fillColor = "#CC002B",
          clusterOptions = markerClusterOptions(),
          group = "Photos",
          layerId = ~imageFileName  # used to look up the row on click, below
        )
      photosLoaded(TRUE)
    }
    
    # Lazily add each environmental raster the first time its group name
    # appears in the checked overlay list, and keep that layer's legend
    # (bottom-right) in sync with its checkbox state on every toggle. The
    # raster itself only needs to be added once -- Leaflet's own layer
    # control automatically shows/hides it after that -- but the legend
    # isn't tied to group visibility automatically, so it has to be
    # explicitly added/removed here each time the box is checked/unchecked.
    for (v in names(env_rasters)) {
      lbl <- env_labels[[v]]
      legend_id <- paste0("envLegend_", v)
      checked <- lbl %in% input$map_groups
      
      if (checked) {
        if (!envLoaded[[v]]) {
          leafletProxy("map") |>
            addRasterImage(env_rasters[[v]], colors = env_palettes[[v]], opacity = 0.8,
                           group = lbl, project = TRUE)
          envLoaded[[v]] <- TRUE
        }
        # addLegend() with an existing layerId replaces that control in
        # place, so calling this every time is safe and won't stack
        # duplicate legends.
        leafletProxy("map") |>
          addLegend(position = "bottomright", pal = env_palettes[[v]],
                    values = env_values[[v]], title = lbl, layerId = legend_id)
      } else if (envLoaded[[v]]) {
        leafletProxy("map") |> removeControl(layerId = legend_id)
      }
    }
  }, ignoreInit = TRUE)
  
  # Popup content is only ever built for the ONE marker that was clicked,
  # not pre-built for all ~150K+ markers up front.
  observeEvent(input$map_marker_click, {
    suppressMapClick(TRUE)  # tell the map-click handler below to skip
    
    click <- input$map_marker_click
    req(click$id)
    row <- photos_unique[photos_unique$imageFileName == click$id, ][1, ]
    req(nrow(row) == 1)
    
    popup_html <- paste0(
      "<img src='", row$photo_url, "' width='150'><br/>",
      "<b>", gsub("_", " ", row$hbSpecies), "</b><br/>",
      "on <i>", row$speciesPlant, "</i><br/>",
      "<a href='", row$photo_url, "' target='_blank'>Open full image</a>"
    )
    
    leafletProxy("map") |>
      addPopups(lng = click$lng, lat = click$lat, popup = popup_html)
  })
  
  # ---- Static map: built ONCE. Do not put input$variable-dependent
  # content here, or the whole map (and its zoom/pan state) gets rebuilt
  # every time the user changes layers.
  output$map <- renderLeaflet({
    leaflet(options = leafletOptions(
      worldCopyJump = FALSE,
      # Static minZoom, set at map initialization rather than computed via
      # JS after render -- the previous onRender()/setTimeout() approach
      # depended on Shiny's dashboard layout having finished sizing the map
      # container at the moment the JS ran, which proved unreliable (map
      # was still showing the whole world in testing). A fixed value here
      # has no such timing dependency and is guaranteed to apply. 3 keeps
      # the Americas roughly filling the view on typical screen widths --
      # raise to 4 if a sliver of Europe/Africa is still visible on very
      # wide screens, or lower to 2 if it feels too zoomed-in on narrow ones.
      minZoom = 3,
      # Declared here too (redundant with setMaxBounds() below, which also
      # still restricts panning) so the bound is baked into the map's
      # initial options rather than depending on a later chained call.
      maxBounds = list(list(-60, -170), list(85, -30)),
      maxBoundsViscosity = 1
    )) |>
      addProviderTiles(providers$CartoDB.Positron,
                       options = providerTileOptions(noWrap = TRUE)) |>
      
      # ---- Environmental layers (static, pre-computed 2015-2024 rasters --
      # replaces the old live FIRMS/NDVI/SPEI WMS tiles). Each is NOT added
      # here anymore -- they start unchecked, and are lazily added (with
      # their project = TRUE reprojection) only when the user actually
      # checks the corresponding box, via the observeEvent(input$map_groups)
      # handler above. This keeps first-render cost independent of how many
      # environmental layers exist, without changing anything the user sees:
      # the checkboxes below still list all 6, still start unchecked, and
      # still show identical data once toggled on.
      addLayersControl(
        overlayGroups = c("Selected raster", "Photos", unname(env_labels)),
        options = layersControlOptions(collapsed = FALSE)
      ) |>
      # Without an actual layer added yet, Leaflet's control still defaults
      # these checkboxes to "checked" unless explicitly told otherwise --
      # this hideGroup() call is what makes them render unchecked at
      # startup (matching reality: nothing has been added for them yet).
      # This also matters functionally, not just visually: an unchecked
      # box means input$map_groups won't list these groups on the initial
      # reactive flush, so the lazy-load observer below won't mistakenly
      # fire for all of them at once during startup.
      hideGroup(c("Photos", unname(env_labels))) |>
      
      # Restrict panning/zooming to the Americas
      setMaxBounds(lng1 = -170, lat1 = -60, lng2 = -30, lat2 = 85) |>
      fitBounds(lng1 = -170, lat1 = -60, lng2 = -30, lat2 = 85)
  })
  
  # ---- Dynamic layer: runs on every variable change, but only touches
  # the raster + legend, leaving the user's current pan/zoom untouched.
  observeEvent(input$variable, {
    pal <- selectedPal()
    r   <- selectedRaster()
    
    leafletProxy("map") |>
      clearGroup("Selected raster") |>
      addRasterImage(r, colors = pal, opacity = 0.8, group = "Selected raster",
                     project = TRUE) |>
      # Bottom-left, kept separate from the environmental legends (bottom-
      # right, added/removed in the map_groups observer above) so both can
      # be visible at once. Uses a fixed layerId so each call replaces the
      # previous bio legend in place instead of stacking duplicates -- and
      # no longer calls clearControls(), which would also have wiped out
      # any currently-active environmental legends.
      addLegend(position = "bottomleft", pal = pal, values = raster_values[[input$variable]],
                title = names(which(variable_choices == input$variable)), layerId = "bioLegend")
  }, ignoreNULL = TRUE)
  
  # ---- Click handler: extract ALL popup_vars at the clicked point ---------
  observeEvent(input$map_click, {
    # A photo marker click already showed its own popup -- don't also show
    # the raster-value popup for this same click.
    if (isTRUE(suppressMapClick())) {
      suppressMapClick(FALSE)
      return(invisible(NULL))
    }
    
    click <- input$map_click
    pt <- terra::vect(matrix(c(click$lng, click$lat), ncol = 2),
                      type = "points", crs = "EPSG:4326")
    
    vals <- sapply(popup_vars, function(v) {
      out <- terra::extract(raster_list[[v]], pt)[1, 2]
      if (is.null(out) || length(out) == 0) NA else round(out, 3)
    })
    
    # Same extraction, for the environmental layers, appended below the
    # biological variables in the same clickTable box. Extracted for all
    # 6 env layers regardless of which are currently toggled visible on
    # the map -- the click table is a data summary independent of what's
    # drawn, same as it already was for the biological variables.
    env_vals <- sapply(names(env_rasters), function(v) {
      out <- terra::extract(env_rasters[[v]], pt)[1, 2]
      if (is.null(out) || length(out) == 0) NA else round(out, 3)
    })
    
    clickValues(rbind(
      data.frame(Variable = popup_labels[popup_vars], Value = vals, row.names = NULL),
      data.frame(Variable = "\u2014 Environmental variables \u2014", Value = NA, row.names = NULL),
      data.frame(Variable = env_labels[names(env_rasters)], Value = env_vals, row.names = NULL)
    ))
    
    # The on-map popup bubble still shows only the biological variables,
    # to keep it from becoming a 13-line wall of text right on the map --
    # the full set (biological + environmental) is available just below
    # in the clickTable box. Easy to extend to the popup too if wanted.
    popup_html <- paste0(
      "<b>", popup_labels[popup_vars], ":</b> ", vals, collapse = "<br/>"
    )
    
    leafletProxy("map") |>
      clearPopups() |>
      addPopups(lng = click$lng, lat = click$lat, popup = popup_html)
    
    # Same grid used to build the photo layer's cell_id, so this always
    # matches whatever cell the user just clicked.
    cell_id <- terra::cellFromXY(TD, cbind(click$lng, click$lat))
    clickCellPhotos(photos[!is.na(photos$cell_id) & photos$cell_id == cell_id, ])
  })
  
  output$clickTable <- renderTable({
    req(clickValues())
    clickValues()
  })
  
  # ---- Phylogenetic subtree for species observed in the clicked cell ------
  output$phyloPlot <- renderPlot({
    cellData <- clickCellPhotos()
    req(cellData)
    
    species_here <- intersect(unique(cellData$hbSpecies), hb_tree$tip.label)
    
    if (length(species_here) == 0) {
      plot.new()
      text(0.5, 0.5, "No hummingbird species recorded in this cell")
    } else if (length(species_here) == 1) {
      plot.new()
      text(0.5, 0.5, gsub("_", " ", species_here))
    } else {
      subtree <- ape::keep.tip(hb_tree, species_here)
      plot(subtree, main = NULL, cex = 0.9)
    }
  })
  
  # ---- Bipartite hummingbird-plant network for the clicked cell -----------
  output$networkPlot <- renderPlot({
    cellData <- clickCellPhotos()
    req(cellData)
    
    edges <- cellData[cellData$scorePlant > 0.8, c("hbSpecies", "speciesPlant")]
    
    if (nrow(edges) == 0) {
      plot.new()
      text(0.5, 0.5, "No high-confidence interactions (score > 0.8) in this cell")
    } else {
      # NOTE: the earlier mar/oma attempts didn't fix the edge-label
      # cropping because that isn't actually a spacing problem -- R clips
      # any drawn element (including rotated text) to the plot's data
      # region by default. A label anchored at the very first/last species'
      # x-position genuinely extends past that boundary once rotated, and
      # gets hard-clipped there regardless of how much blank margin exists
      # around it. xpd=NA disables that clipping so labels can draw
      # anywhere on the device, which is the actual fix for this.
      op <- par(oma = c(0, 5, 0, 5), xpd = NA)
      on.exit(par(op), add = TRUE)
      
      # NOTE: this bipartite version (2.20) uses the older plotweb() API --
      # confirmed by the method="normal" argument, which only exists on the
      # old API. On this version, "labsize" (not "text.size") is the
      # correct label-size parameter.
      web <- table(gsub("_", " ", edges$hbSpecies), edges$speciesPlant)
      bipartite::plotweb(web, method = "normal",
                         col.high = "#009180", col.low = "#FFC20F",
                         text.rot = 90, labsize = 1.6,
                         high.lablength = 20, low.lablength = 20)
    }
  })
}

# ---- Run ------------------------------------------------------------------
shinyApp(ui = ui, server = server)