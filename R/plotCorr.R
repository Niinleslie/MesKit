plotCorr <- function(corrMat, use.circle = TRUE){
    corrMat[lower.tri(corrMat)] <- 0
    maxCorr <- max(corrMat[corrMat!=1])
    TSBs <- colnames(corrMat)

    col_fun <- circlize::colorRamp2(pretty(c(0,maxCorr),4),
                                    c("#f1eef6", "#d0d1e6", "#a6bddb", "#74a9cf", "#2b8cbe", "#023858"))
    
    p <- ComplexHeatmap::Heatmap(
        corrMat,
        col = col_fun,
        #rect_gp = grid::gpar(type = "none"),
        rect_gp = grid::gpar(type = "none", col = "white", lwd = 2),
        cell_fun = function(j, i, x, y, width, height, fill) {
            
            if (j > i) {
                if (use.circle) {
                    if(i == 1 & j != i+1) {
                        grid::grid.segments(x, y-height*0.5, x, y)
                    } else if(i != 1 & j == i+1) {
                        grid::grid.segments(x, y, x, y+height*0.5)
                    } else if(i != 1){
                        grid::grid.segments(x, y-height*0.5, x, y+height*0.5)
                    }
                    
                    if(j == ncol(corrMat) & j != i+1){
                        grid::grid.segments(x-width*0.5, y, x, y)        
                    } else if(j < ncol(corrMat) & j == i+1) {
                        grid::grid.segments(x, y, x+width*0.5, y)
                    } else if (j != i+1){
                        grid::grid.segments(x-width*0.5, y, x+width*0.5, y)
                    }
                   
                    grid::grid.circle(
                        x = x,
                        y = y,
                        r = corrMat[i, j] / 2 * min(grid::unit.c(width, height)),
                        gp = grid::gpar(
                            fill = col_fun(corrMat[i, j]),
                            col = NA
                        )
                    )
                } else {
                    grid::grid.rect(
                        x = x, y = y,
                        width = 0.99*width,
                        height = 0.99*height,
                        gp = grid::gpar(
                            col = "white",
                            fill = col_fun(corrMat[i, j])
                        )
                    )
                }
                
                grid::grid.text(
                    sprintf("%.2f", corrMat[i, j]),
                    x, y,
                    gp = grid::gpar(fontsize = 10, col = "white")
                )
            }
            
            if (j == i) {
                grid::grid.text(
                    TSBs[i],
                    x = x,
                    y = y,
                    gp = grid::gpar(fontsize = 10)
                )
            }
        },
        show_heatmap_legend = FALSE,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        width = grid::unit(12, "cm"),
        height = grid::unit(12, "cm")
    )
    
    ComplexHeatmap::draw(p,
                         heatmap_legend_list = list(
                             ComplexHeatmap::Legend(
                                 col_fun = col_fun,
                                 title = "Similarity",
                                 legend_height = grid::unit(6, "cm")
                             )
                         )
    )
    
}