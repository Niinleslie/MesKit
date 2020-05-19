

change_name.character <- function(list, old, new) {
    list[list %in% old] <- new
    list
}

change_name.default <- function(list, old, new) {
    nam <- names(list)
    nam[nam %in% old] <- new
    names(list) <- nam
    list
}

change_name.NULL <- function(list, old, new) {
    NULL
}


remove_new <- function(aes) {
    stringi::stri_replace_all(aes, "", regex = "(_new)*")
}

change_name <- function(list, old, new) {
    UseMethod("change_name")
}

bump_aes.Scale <- function(layer, new_aes) {
    old_aes <- layer$aesthetics[remove_new(layer$aesthetics) %in% new_aes]
    new_aes <- paste0(old_aes, "_new")
    
    layer$aesthetics[layer$aesthetics %in% old_aes] <- new_aes
    
    if (is.character(layer$guide)) {
        layer$guide <- match.fun(paste("guide_", layer$guide, sep = ""))()
    }
    layer$guide$available_aes[layer$guide$available_aes %in% old_aes] <- new_aes
    layer
}

bump_aes.Layer <- function(layer, new_aes) {
    original_aes <- new_aes
    
    old_aes <- names(layer$mapping)[remove_new(names(layer$mapping)) %in% new_aes]
    new_aes <- paste0(old_aes, "_new")
    
    old_geom <- layer$geom
    
    old_setup <- old_geom$handle_na
    new_setup <- function(self, data, params) {
        colnames(data)[colnames(data) %in% new_aes] <- original_aes
        old_setup(data, params)
    }
    
    new_geom <- ggplot2::ggproto(paste0("New", class(old_geom)[1]), old_geom,
                                 handle_na = new_setup)
    
    new_geom$default_aes <- change_name(new_geom$default_aes, old_aes, new_aes)
    new_geom$non_missing_aes <- change_name(new_geom$non_missing_aes, old_aes, new_aes)
    new_geom$required_aes <- change_name(new_geom$required_aes, old_aes, new_aes)
    new_geom$optional_aes <- change_name(new_geom$optional_aes, old_aes, new_aes)
    
    layer$geom <- new_geom
    
    old_stat <- layer$stat
    
    old_setup2 <- old_stat$handle_na
    new_setup <- function(self, data, params) {
        colnames(data)[colnames(data) %in% new_aes] <- original_aes
        old_setup2(data, params)
    }
    
    new_stat <- ggplot2::ggproto(paste0("New", class(old_stat)[1]), old_stat,
                                 handle_na = new_setup)
    
    new_stat$default_aes <- change_name(new_stat$default_aes, old_aes, new_aes)
    new_stat$non_missing_aes <- change_name(new_stat$non_missing_aes, old_aes, new_aes)
    new_stat$required_aes <- change_name(new_stat$required_aes, old_aes, new_aes)
    new_stat$optional_aes <- change_name(new_stat$optional_aes, old_aes, new_aes)
    
    layer$stat <- new_stat
    
    layer$mapping <- change_name(layer$mapping, old_aes, new_aes)
    layer
}

bump_aes.list <- function(layer, new_aes) {
    old_aes <-  names(layer)[remove_new(names(layer)) %in% new_aes]
    new_aes <- paste0(old_aes, "_new")
    
    names(layer)[names(layer) %in% old_aes] <- new_aes
    layer
}


bump_aes <- function(layer, new_aes) {
    UseMethod("bump_aes")
}

new_scale <- function(new_aes) { 
 structure(ggplot2::standardise_aes_names(new_aes), class = "new_aes")
}

new_scale_fill <- function() {
 new_scale("fill")
}

new_scale_color <- function() {
 new_scale("colour")
}

new_scale_colour <- function() {
 new_scale("colour")
}



ggplot_add.new_aes <- function(object, plot, object_name) {
 plot$layers <- lapply(plot$layers, bump_aes, new_aes = object)
 plot$scales$scales <- lapply(plot$scales$scales, bump_aes, new_aes = object)
 plot$labels <- bump_aes(plot$labels, new_aes = object)
 plot
}




