plotCorr <- function(corrMat, use.circle = TRUE, title = NULL, number.cex = 8, number.col = "#c05836"){
    corrMat[lower.tri(corrMat)] <- 0
    maxCorr <- max(corrMat[corrMat!=1])
    minCorr <- min(corrMat[corrMat!=1 & corrMat!=0])
    #maxCorr <- max(corrMat)
    TSBs <- colnames(corrMat)
    
    
    #col_fun <- circlize::colorRamp2(pretty(c(0, maxCorr),6),
    #                                c("#f1eef6", "#d0d1e6", "#a6bddb", "#74a9cf", "#2b8cbe", "#023858"))
    significant_digit <- gsub(pattern =  "0\\.0*","",as.character(minCorr))
    digits <- nchar(as.character(minCorr)) - nchar(significant_digit) 
    col_fun <- circlize::colorRamp2(breaks = round( seq(0, maxCorr, length = 6),digits = digits),
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
                        r = corrMat[i, j] / maxCorr/ 2 * min(grid::unit.c(width, height))*0.7,
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
                            col = "#fdc086",
                            fill = col_fun(corrMat[i, j])
                        )
                    )
                }
                
                grid::grid.text(
                    round(corrMat[i, j],digits = 3) ,
                    x, y,
                    gp = grid::gpar(fontsize = number.cex, fontface = "bold", col = number.col)
                )
            }
            
            if (j == i) {
                grid::grid.text(
                    TSBs[i],
                    x = x,
                    y = y,
                    gp = grid::gpar(fontsize = 11, fontface = "bold")
                )
            }
        },
        column_title = title,
        column_title_gp = grid::gpar(fontsize = 13.5, fontface = "bold", col = "black"),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        width = grid::unit(12, "cm"),
        height = grid::unit(12, "cm"),
        heatmap_legend_param = list(
                title = NULL,
                col_fun = col_fun,                
                legend_height = grid::unit(4, "cm")
            )

    )
    return(p)
}
    
