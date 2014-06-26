plotCompartment = function(compartment, main = "Epidemic Compartment",
                           xlab = "Time Point", ylab = "Spatial Location", 
                           col=heat.colors(20), transpose = FALSE)
{
    plot.new()
    layout.matrix = matrix(c(1,2), nrow = 1)
    layout(layout.matrix, c(8,2), c(4,4))
    par(xaxt="n")
    par(yaxt="n")
    if (transpose)
    {
        compartment = t(compartment)
    }
    image(compartment, main=main, xlab = xlab, ylab = ylab,col=col,)
    axis(side = 1, at = floor(seq(1, nrow(compartment), length = 10)))

    # Make Legend
    par(bty="n")
    par(xaxt="n")
    par(yaxt="n")
    par(mar=c(5, 0, 4, 4) + 0.1)
    image(matrix(seq(min(compartment, max(compartment)), length = 20), nrow = 1))
    par(yaxt="s")
    axis(side =4, at = seq(0, 1, length = 10), labels = round(seq(min(compartment, max(compartment)), length = 10)))
}
