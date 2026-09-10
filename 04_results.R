library(dplyr)
library(tidyr)
library(ggplot2)



# perf.df %>% 
#   filter(condition == 19 & lvl.redun == 0.9) %>% 
#   head(30)


# Average across the 30 iterations -> one point per condition × sigma × solution
# plot.df = data.df %>%
#   group_by(condition, sigma.type, solution) %>%
#   summarize(across(all_of(metric.cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
#   pivot_longer(all_of(metric.cols), names_to = "measure", values_to = "value")


metric.cols = c("correlation", "wd.bias", "apl.bias", "strgth.bias",
                "expctinflu.bias", "closeness.bias", "betweenness.bias", "RMSE")

truth.labels = c(
  redundant = "Common Factor Generating Model",
  composite = "Composite Generating Model",
  pseudo    = "Unique Generating Model"
)


perf.df %>%
  filter(abs(strgth.bias) > 2) %>%
  count(solution, sigma.type, lvl.redun, lat.cor, sparsity, sort = TRUE)



# Heywood cases

perf.df = perf.df %>% 
  mutate(mag = abs(strgth.bias) + 1e-6) # add a new column called mag 
  # Take absolute of strgth.bias and add 1e-6 because we are going to put values on log scale
  

perf.df %>%
  mutate(mag = abs(strgth.bias) + 1e-6) %>%
  ggplot(aes(mag)) +
  geom_histogram(bins = 100) +
  scale_x_log10() + # each step is multiplication by 10. Log axis stretches the small values so they are visible
  labs(x = "|strength bias| (log scale)") 

quantile(abs(perf.df$strgth.bias), c(.90, .95, .99, .995, .999), na.rm = TRUE)
# Reports the value below which the designated percentages of the magnitudes fall.


v <- sort(abs(perf.df$strgth.bias[is.finite(perf.df$strgth.bias)]))
tail_v <- v[v > quantile(v, .98)]
tail_v[which.max(tail_v[-1] / head(tail_v, -1))] 









# dat = perf.df %>%
#   pivot_longer(all_of(metric.cols), names_to = "measure", values_to = "value") %>%
#   filter(measure == "strgth.rbias", sigma.type == "redundant") %>%
#   group_by(across(all_of(c("solution", "sparsity", "p", "lvl.redun")))) %>%
#   summarize(
#     mean = mean(value, na.rm = TRUE),
#      se   = sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value))),
#     .groups = "drop"
#   )



# Emorie Beck, PhD, data viz theme 
my_theme = function(){
  theme_light() + 
    theme(
      legend.position = "right"
      , legend.title = element_text(face = "bold", size = rel(1.6))
      , legend.text = element_text(face = "italic", size = rel(1.6))
      , axis.text = element_text(face = "bold", size = rel(1.6), color = "black")
      , axis.title = element_text(face = "bold", size = rel(1.6))
      , plot.title = element_text(face = "bold", size = rel(1.6), hjust = .5)
      , plot.subtitle = element_text(face = "italic", size = rel(1.6), hjust = .5)
      , strip.text = element_text(face = "bold", size = rel(1.6), color = "white")
      , strip.background = element_rect(fill = "black")
    )
}





# 
# plot.data = function(m, sig, row.var = "p", col.var = "lvl.redun") {
# 
#   d = data.df %>%
#     pivot_longer(all_of(metric.cols), names_to = "measure", values_to = "value") %>%
#     filter(measure == m, sigma.type == sig) %>%
#     group_by(across(all_of(c("lat.cor", "solution", "sparsity", row.var, col.var)))) %>%
#     summarize(
#       mean = mean(value, na.rm = TRUE),
#       se   = sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value))),
#       .groups = "drop"
#     )
#   d
# }
# 
# for (sig in c("redundant")) {
#   print(plot.data(m = "apl.rbias", sig = sig), n = 36+72)
# }


cutoff = 0.2

# One figure = one measure × one sigma.type
plot.measure = function(m, sig, row.var = "p", col.var = "lvl.redun") {
  
  d = perf.df %>%
    pivot_longer(all_of(metric.cols), names_to = "measure", values_to = "value") %>%
    filter(measure == m, sigma.type == sig) %>%
    # filter(abs(value) < cutoff) %>%
    mutate(value = ifelse(is.finite(value), value, NA_real_)) %>% 
    group_by(across(all_of(c("solution", "sparsity", row.var, col.var)))) %>%
    summarize(
      mean = mean(value, na.rm = TRUE),
      se   = sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value))),
      .groups = "drop"
    ) %>% 
    mutate(lvl.redun = factor(lvl.redun, 
                              levels = c(0.7, 0.8, 0.9),
                              labels = c("Low Redundancy",
                                         "Medium Redundancy",
                                         "High Redundancy")))
  
  p.out = ggplot(d, aes(x = sparsity, y = mean, color = solution, group = solution,
                        shape = solution)) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 2.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                  linewidth = .5,
                  width = .06) +
    scale_color_manual(name = "Solution",
      values = c(
      "composite" = "#E69F00",  # orange
      "lnm"       = "#009E73",  # bluish green
      "removal"   = "#0072B2"   # blue
    )) +
    scale_shape_manual(name = "Solution", 
                       values = c(
      "composite" = 15,
      "lnm" = 16, 
      "removal" = 17)) +
    guides(color = guide_legend(title = "Solution")) +
    facet_grid(rows = vars(.data[[row.var]]),
               cols = vars(.data[[col.var]])) +
    labs(title = truth.labels[[sig]],
         x = "Sparsity",
         y = "Bias",
          color = "Solution") +
    my_theme() + 
    theme(strip.background = element_rect(fill = "black"),
          strip.text = element_text(color = "white", face = "bold"))
  p.out
}



# One measure, one sigma:
# plot.measure("betweenness.rbias", "pseudo")
# 
# Three figures for one measure (one per sigma.type):

for (sig in c("redundant", "composite", "pseudo")) {
  png(filename = sprintf("strength_%s_sparsityV8.png", sig),
      width = 12, height = 8, res = 300, units = "in")
  print(plot.measure("strgth.bias", sig))
  dev.off()
}



# metric.cols = c("correlation", "wd.bias", "apl.bias", "strgth.bias",
#                 "expctinflu.bias", "closeness.bias", "betweenness.bias", "RMSE")

# # Everything: all 8 measures × 3 sigma types = 24 figures
# for (m in metric.cols) {
#   for (sig in c("redundant", "composite", "pseudo")) {
#     print(plot.measure(m, sig))
#   }
# }
# unique(data.df$solution)


# perf.df %>%
#   filter(sigma.type == "redundant", solution %in% c("composite","removal")) %>%
#   group_by(p, lvl.redun, lat.cor, solution) %>%
#   summarise(
#     n      = sum(!is.na(wd.rbias)),
#     mean   = mean(wd.rbias, na.rm = TRUE),
#     sd     = sd(wd.rbias,   na.rm = TRUE),
#     cv     = sd / abs(mean),                       # spread relative to the mean
#     .groups = "drop"
#   ) %>%
#   arrange(desc(sd)) %>%
#   head(12)

