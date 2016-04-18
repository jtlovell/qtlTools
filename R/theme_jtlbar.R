#' @import ggplot2
#' @export
theme_jtlbar <- function(base_size = 10, base_family = "")
{
  txt.8 <- element_text(size = 8, colour = "black", face = "plain")
  txt.10 <- element_text(size = 10, colour = "black", face = "plain")
  axrange <- list(y=c(50,90),x=c(2,5))
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      legend.key = element_blank(),
      strip.text = txt.10,

      text = txt.8,
      plot.title = txt.10,

      axis.title = txt.10,
      axis.text = txt.8,

      legend.title = txt.10,
      legend.text = txt.8 ,

      strip.background = element_blank(),

      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_line(color="grey", size=.5),

      panel.background = element_rect(color="black", fill=NA),
      panel.border = element_blank(),

      axis.ticks.length = unit(.15, "cm"),
      axis.ticks = element_line(size=.5, color="black"))

}
