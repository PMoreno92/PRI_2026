#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(shinydashboard)
library(leaflet)
library(raster)
library(terra)
library(ggplot2)

americas <- ext(-170, -30, -60, 85)

TD = crop(rast("../results/hb_raster_qTD_est95_wgs84.tiff"), americas)
PD = crop(rast("../results/hb_raster_qPD_3S_est95_wgs84.tiff"), americas)
FD = crop(rast("../results/hb_raster_qFD_3S_est95_wgs84.tiff"), americas)
connectance = crop(rast("../results/hb_80_raster_connectance_5S_wgs84.tiff"), americas)
compdiv = crop(rast("../results/hb_80_raster_compdiv_5S_wgs84.tiff"), americas)
cc = crop(rast("../results/hb_80_raster_cc_5S_wgs84.tiff"), americas)
NODF = crop(rast("../results/hb_80_raster_NODF_5S_wgs84.tiff"), americas)
H2 = crop(rast("../results/hb_80_raster_H2_5S_wgs84.tiff"), americas)
cchb = crop(rast("../results/hb_80_raster_cchb_5S_wgs84.tiff"), americas)
ccpl = crop(rast("../results/hb_80_raster_ccpl_5S_wgs84.tiff"), americas)
selelev_mean = crop(rast("../results/selelev_mean.tiff"), americas)

# Define UI for application that draws a histogram
ui <- dashboardPage(
  dashboardHeader(disable = TRUE),
  dashboardSidebar(radioButtons("variable", "Select Variable", 
                                choices = list("Taxonomic Diversity (SC = 95 %)" = "TD",
                                               "Phylogenetic Diversity (SC = 95 %)" = "PD",
                                               "Functional Diversity (SC = 95 %)" = "FD",
                                               "Network: connectance" = "connectance",
                                               "Network: compartment diversity" = "compdiv",
                                               "Network: cluster coefficient" = "cc",
                                               "Network: NODF" = "NODF",
                                               "Network: H2" = "H2",
                                               "Network: cluster coefficient (hummingbirds)" = "cchb",
                                               "Network: cluster coefficient (plants)" = "ccpl"))),
  
  dashboardBody(leafletOutput(outputId = "map")),
  title = "Interactive Maps"
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  # Map selection (for biological variables)
  ## Add the options
  selectedRaster <- reactive({
    switch(input$variable,
           "TD" = TD,
           "PD" = PD,
           "FD" = FD,
           "connectance" = connectance,
           "compdiv" = compdiv,
           "cc" = cc,
           "NODF" = NODF,
           "H2" = H2,
           "cchb" = cchb,
           "ccpl" = ccpl)})
  
  # Add the palettes
  pal_TD <- colorNumeric(palette = c("#FFFFFF", "#0058A1"), domain = TD)
  pal_PD <- colorNumeric(palette = c("#FFFFFF", "#B80007"), domain = selectedRaster$PD)
  pal_FD <- colorNumeric(palette = c("#FFFFFF", "#591f63"), domain = selectedRaster$FD)
  pal_connectance <- colorNumeric(palette = c("#FFFFFF", "#009590"), domain = selectedRaster$connectance)
  pal_compdiv <- colorNumeric(palette = c("#FFFFFF", "#009590"), domain = selectedRaster$compdiv)
  pal_cc <- colorNumeric(palette = c("#FFFFFF", "#009590"), domain = selectedRaster$cc)
  pal_NODF <- colorNumeric(palette = c("#FFFFFF", "#009590"), domain = selectedRaster$NODF)
  pal_H2 <- colorNumeric(palette = c("#FFFFFF", "#009590"), domain = selectedRaster$H2)
  pal_cchb <- colorNumeric(palette = c("#FFFFFF", "#009590"), domain = selectedRaster$cchb)
  pal_ccpl <- colorNumeric(palette = c("#FFFFFF", "#009590"), domain = selectedRaster$ccpl)
  
  
  observeEvent(input$map_click, {
    click <- input$map_click
    pt <- terra::vect(matrix(c(click$lng, click$lat), ncol = 2),
                      type = "points", crs = "EPSG:4326")
    value <- terra::extract(selectedRaster(), pt)[1,2]
    
    leafletProxy("map") |>
      clearPopups() |>
      addPopups(lng = click$lng, lat = click$lat,
                popup = paste0("<b>Value:</b> ", round(value, 3)))})
  
  output$map <- renderLeaflet({
    
    leaflet(options = leafletOptions(worldCopyJump = FALSE)) |>
      addProviderTiles(providers$CartoDB.Positron,
                       options = providerTileOptions(noWrap = TRUE)) |>
      
    # NASA FIRMS (active fires)
    addWMSTiles(baseUrl = "https://firms.modaps.eosdis.nasa.gov/mapserver/wms",
      layers = "fires_viirs_snpp_24",
      options = WMSTileOptions(format = "image/png",
                               transparent = TRUE),
      attribution = "NASA FIRMS",
      group = "NASA FIRMS (24 h)") |>
      
    # NDVI
    addWMSTiles(baseUrl = "https://gibs.earthdata.nasa.gov/wms/epsg4326/best/wms.cgi",
                layers = "MODIS_Terra_NDVI_8Day",
                options = WMSTileOptions(format = "image/png",
                                         transparent = TRUE),
                attribution = "NASA GIBS", group = "NDVI") |>
    
    # CSIC SPEI
    addWMSTiles(baseUrl = "https://spei.csic.es/geoserver/SPEI/wms",
                layers = "spei_latest",
                options = WMSTileOptions(format = "image/png",
                                         transparent = TRUE),
                attribution = "CSIC SPEI", group = "SPEI") |>
      
      # Biodiversity raster
      addRasterImage(selectedRaster(),
                     group = "Selected raster") |>
      
      addLayersControl(overlayGroups = c("Selected raster", "NASA FIRMS (24 h)",
                                         "NDVI","SPEI"),
                       options = layersControlOptions(collapsed = FALSE)) |>
      
      addLegend(position = "bottomright", opacity = 1)
      
      fitBounds(-170, -60, -30, 85)
    
    
  })
  
}


# Run the application 
shinyApp(ui = ui, server = server)
