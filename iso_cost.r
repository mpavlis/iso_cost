library(dplyr)
library(raster)
library(rgdal)

source("~/Dropbox/git/Route/route.r")
source("~/Dropbox/git/locoh/locoh.r")
source("~/Dropbox/git/Route/sf_functions.r")
source("~/Dropbox/git/postgis/postgis_functions.r")

##### 1. points_to_graph_nodes ######################################################################

points_to_graph_nodes <- function(nodes_mat, coords_mat){
  graph_node_id <- knnx.index(nodes_mat, coords_mat, 1)
  as.integer(graph_node_id)
}

##### 2. graph_nodes ################################################################################

graph_nodes <- function(G){
  v_nodes <- do.call(cbind, tstrsplit(V(G)$name, " "))
  storage.mode(v_nodes) <- "numeric"
  v_nodes
}

##### 3. isocost_nodes ##############################################################################

isocost_nodes <- function(from_node, G, cost, to_nodes = V(G), ...){
  
  ################################# Function ###########################################
  
  get_cost <- function(G, from_node, to_nodes, cost, ...){
    # all_dist <- distances(graph = G, v = from_node, to = to_nodes, ...)
    # id_sele <- all_dist <= cost
    # cbind(as.integer(from_node),
    #       as.integer(to_nodes[id_sele]),
    #       all_dist[id_sele])
    as.numeric(distances(graph = G, v = from_node, to = to_nodes, ...))
  }
  
  ############################### Main Function ########################################
  
  # weighted <- attributes(G)$weighted
  # if (weighted == "time"){
  #   within_dist <- ceiling(attributes(G)$max_speed * 1000 * cost / 60 )
  # } else if (weighted == "distance"){
  #   within_dist <- cost * 1000
  # } else {
  #   stop("the graph has to be weighted")
  # }
  # 
  # if (clip_points){
  #   area <- st_buffer(all_nodes_as_points[from_node], within_dist)
  #   nodes_as_points <- nodes_as_points[lengths(st_within(nodes_as_points, area)) > 0, ]
  # }
  
  # get_cost(G, from_node = from_node, to_nodes = to_nodes[nodes_within], cost = cost, ...)
  
  get_cost(G = G, from_node = from_node, to_nodes = to_nodes, cost = cost, ...)
  
}

extract_road <- function(point, con, road_table, srid = 27700, geom_column = "geom", speed_column = "maxspeed", radius = 112000){
  
  st_read_db(conn = con,
             query = paste0("SELECT t1.", speed_column, ", (ST_Dump(t1.", geom_column, ")).", geom_column, " geometry FROM ", road_table, " t1, ", 
                            "ST_SetSRID(ST_Buffer(ST_Point(", point[1], ", ", point[2], "), ", radius, "), ", srid,") t2",
                            " WHERE ST_Intersects(t1.", geom_column, ", t2)"))
  
}

##### 4. combine_polys_to_raster ####################################################################

sp_to_raster <- function(sp_obj, raster_res, fun = "min", field = NULL){
  
  raster_obj <- raster()
  extent(raster_obj) <- extent(sp_obj)
  res(raster_obj) <- raster_res
  raster_obj <- rasterize(sp_obj, raster_obj, field = field, fun = fun)
  raster_obj
  
  # focal(x = raster_obj, w = weights_mat, fun=modal, NAonly = T, pad = T, na.rm = T)
  
}

##### 5. isocost_polys ##############################################################################

isocost_polys <- function(from_point, road, con, speed_limit_column = NULL, srid = 27700,
                          cost = 60, breaks = c(0,15,30,45,60), k = 20, k_max = 20,
                          raster_fun = "min", raster_res = 1000, weights_mat = matrix(1,3,3), ...){
  if (is(from_point, "vector")){
    if (length(from_point) != 2){
      stop("the from_point vector should have exactly two elements, easting and northing")
    }
    point_coords <- rbind(from_point)
  } else if (is(from_point, "matrix")){
    if (ncol(from_point) != 2 | nrow(from_point) != 1){
      stop("the from_point matrix should have a single row with columns, easting and northing")
    }
    point_coords <- from_point
  }
  
  g_road <- lines_to_graph(lines_sf = road, weighted_graph = T, speed_limit_column = speed_limit_column)
  
  nodes_as_points <- SpatialPoints(graph_nodes(g_road), proj4string = CRS(paste0("+init=epsg:",srid)))
  nodes_as_points$id <- 1:length(nodes_as_points)
  node <- points_to_graph_nodes(coordinates(nodes_as_points), point_coords)
  nodes_as_points$cost <- isocost_nodes(from_node = node, G = g_road, cost = cost, to_nodes = nodes_as_points$id)
  nodes_as_points <- nodes_as_points[nodes_as_points$cost <= cost, ]
  # nodes_lookup <- as.data.frame(isocost_nodes(from_node = node, G = G, all_nodes_as_points = nodes_as_points, cost = cost, ...))
  # nodes_lookup <- as.data.frame(isocost_nodes(from_node = node, G = g_road, cost = cost, ...))
  # names(nodes_lookup) <- c("from_node", "to_node", "cost")
  break_ids <- 1:length(breaks[breaks!=0])
  breaks <- sort(breaks)
  # nodes_lookup$breaks <- as.integer(cut(nodes_lookup$cost, breaks, include.lowest=T))
  nodes_as_points$breaks <- as.integer(cut(nodes_as_points$cost, breaks, include.lowest=T))
  
  # nodes_as_points_sele <- st_as_sf(remove.duplicates(nodes_as_points[nodes_as_points$id %in% nodes_lookup$to_node, ], rm_dupl_dist))
  # nodes_as_points_sele <- inner_join(x = nodes_as_points_sele, 
  #                                    y = nodes_lookup[,c("to_node", "cost", "breaks")],
  #                                    by = c("id" = "to_node"))
  # st_crs(nodes_as_points_sele) <- srid
  # nodes_as_points_sele$id <- as.integer(nodes_as_points_sele$id)
  
  raster_obj <- sp_to_raster(sp_obj = nodes_as_points, raster_res = raster_res, fun = modal, field = "breaks")
  raster_pnts <- as.data.frame(rasterToPoints(raster_obj))
  coordinates(raster_pnts) <- ~x+y
  proj4string(raster_pnts) <- proj4string(nodes_as_points)
  raster_pnts$id <- 1:length(raster_pnts)
  names(raster_pnts)[1] <- "breaks"
  raster_pnts$breaks <- as.integer(raster_pnts$breaks)
  
  road_nodes_table <- paste0("road_nodes_", node)
  road_nodes_shp <- paste0(road_nodes_table,".shp")
  # st_write_db(conn = con, obj = st_as_sf(raster_pnts), table = road_nodes_table, geom_name = "geom", drop = T)
  st_write(st_as_sf(raster_pnts), road_nodes_shp)
  import_or_append(con, working_dirs = getwd(), shp_names = road_nodes_shp, srid = srid, table_name = road_nodes_table)
  dbExecute(con, paste("CREATE INDEX ON", road_nodes_table, "(id)"))
  # dbExecute(con, paste("CREATE INDEX ON", road_nodes_table, "USING gist(geom)"))
  
  # create_postgis_locoh_func(con = con,
  #                           clusters_table = road_nodes_table,
  #                           geom_field = "geom",
  #                           clusters_field = "breaks",
  #                           points_id_field = "id")
  
  polys <- do.call(rbind, lapply(break_ids, get_locoh, con = con, clusters_table = road_nodes_table, clusters_field = "breaks", geom_field = "geom", points_id_field = "id", rm_holes = F, k = k, k_max = k_max))
  names(polys)[1] <- "breaks"
  
  breaks <- breaks[breaks != 0]
  check_geoms <- is.na(st_dimension(st_geometry(polys)))
  if (any(check_geoms)){
    id_rm <- which(is.na(check_geoms))
    polys <- polys[-id_rm, ]
    breaks <- breaks[-id_rm]
  }
  
  # empty_geom <- list()
  # x <- 1
  # for (i in 1:length(polys_list)){
  #   if (is.na(st_dimension(polys_list[[x]]))){
  #     empty_geom[[as.character(node)]] <- c(empty_geom[[as.character(node)]], i)
  #     polys_list[[x]] <- NULL
  #   } else {
  #     polys_list[[x]] <- st_as_text(polys_list[[x]])
  #     x <- x + 1
  #   }
  # }
  # 
  # locoh_sfc <- st_as_sfc(polys_list)
  
  # if (! is.null(empty_geom[[as.character(node)]])){
  #   breaks <- breaks[-empty_geom[[as.character(node)]]]
  # }
  # locoh_sf <- data.frame(breaks = breaks)
  # locoh_sf$geometry <- locoh_sfc
  # locoh_sf <- st_as_sf(locoh_sf)
  # st_crs(locoh_sf) <- srid
  
  locoh_raster <- sp_to_raster(sp_obj = as(polys, "Spatial"),
                               raster_res = 200,
                               fun = raster_fun,
                               field = "breaks")
  locoh_raster2 <- focal(x = locoh_raster, w = weights_mat, fun=modal, NAonly = T, pad = T, na.rm = T)
  
  # drop_postgis_locoh_func(con, clusters_table = road_nodes_table)
  
  dbExecute(con, paste("DROP TABLE IF EXISTS", road_nodes_table))
  
  unlink(paste0(road_nodes_table, "*"))
  
  st_as_sf(rasterToPolygons(x = locoh_raster2, na.rm = T, dissolve = T))
  
}
