#' @aliases theme_jtl theme_jtlbar
#' @title simple custom themes
#'
#' @description
#' \code{theme_jtl}  two themes for ggplot similar to theme bw, but better (i think). One is
#' appropriate for most types of plots (theme_jtl). The other is better for barplots (theme_jtlbar)
#' @export
theme_jtl <- function()
{
  base_size = 10
  base_family = ""
  txt.8 <- element_text(size = 8, colour = "black", face = "plain")
  txt.10 <- element_text(size = 10, colour = "black", face = "plain")

  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      legend.key = element_blank(),
      strip.background = element_blank(),

      strip.text = txt.10,
      text = txt.8,
      plot.title = txt.10,

      axis.title = txt.10,
      axis.text = txt.8,

      legend.title = txt.10,
      legend.text = txt.8 ,

      strip.background = element_blank(),

      panel.grid.major= element_blank(),
      panel.grid.minor = element_blank(),

      panel.background = element_blank(),
      panel.border = element_rect(color="black", size=.5, fill=NA),

      axis.ticks.length = unit(.15, "cm"),
      axis.ticks = element_line(size=.5, color="black"))
}

theme_jtlbar <- function()
{
  base_size = 10
  base_family = ""
  txt.8 <- element_text(size = 8,
                        colour = "black",
                        face = "plain")
  txt.10 <- element_text(size = 10,
                         colour = "black",
                         face = "plain")
  txt.12 <- element_text(size = 12,
                         colour = "black",
                         face = "plain")
  axrange <- list(y = c(50, 90), x = c(2, 5))
  theme_classic(base_size = base_size,
                base_family = base_family) +
    theme(legend.key = element_blank(),
          strip.text = txt.10,

          text = txt.8,
          plot.title = txt.12,
          axis.title = txt.10,
          axis.text = txt.8,
          legend.title = txt.10,
          legend.text = txt.8,

          strip.background = element_blank(),

          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_line(color = "grey",   size = 0.5),

          panel.background = element_blank(),
          panel.border = element_blank(),

          axis.ticks.length = unit(0.15, "cm"),
          axis.ticks = element_line(size = 0.5,  color = "grey"))
}
