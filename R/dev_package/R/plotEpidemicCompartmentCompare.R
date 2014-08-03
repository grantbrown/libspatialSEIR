plotTwoCompartments = function(compartment1, compartment2, main1 = "Epidemic Compartment",
                                main2 = "EpdiemicCompartment",zlim = NA,xlab = "Time Point",
                                ylab = "Spatial Location", col=heat.colors(20), transpose = FALSE)
{
    plot.new()
    layout.matrix = matrix(c(1,3,2,4), nrow = 2)
    layout(layout.matrix, widths = c(8,2), heights= c(2,2))
    if (transpose)
    {
        compartment1 = t(compartment1)
        compartment2 = t(compartment2)
    }
    if (length(zlim) == 1 && is.na(zlim))
    {
        zlim = c(min(min(compartment1), compartment2), max(max(compartment1), max(compartment2)))
    }


    compList = list(compartment1, compartment2)
    titleList = list(main1, main2)
    i = 1
    for (compartment in compList)
    {
        par(bty="o")
        par(xaxt="n")
        par(yaxt="n")
        par(mar=c(5, 4, 4, 1) + 0.1)
        image(compartment, main=titleList[[i]], xlab = xlab, ylab = ylab,col=col,zlim = zlim)
        axis(side = 1, at = floor(seq(1, nrow(compartment), length = 10)))

        # Make Legend
        par(bty="n")
        par(xaxt="n")
        par(yaxt="n")
        par(mar=c(5, 1, 4, 4) + 0.1)
        image(matrix(seq(min(zlim, max(zlim)), length = 20), nrow = 1))
        par(yaxt="s")
        axis(side =4, at = seq(0, 1, length = 10), labels = round(seq(min(zlim), max(zlim), length = 10), 2))
        i = i + 1
    }

}
