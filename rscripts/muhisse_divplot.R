library(gghisse)
library(ggtree)

tree_flip <- function (ggtree_object, show_tip_labels, tree_layout, tree_direction, 
          time_axis_ticks, agemax, tip_size, offset) 
{
  if (tree_layout %in% c("rectangular", "slanted")) {
    if (tree_direction == "up") {
      ggtree_object <- ggtree_object + coord_flip() + 
        scale_x_continuous(expand = c(add = c(0, 0.01)), 
                           breaks = pretty(c(0, agemax), n = time_axis_ticks), 
                           labels = rev(pretty(c(0, agemax), n = time_axis_ticks))) + 
        theme(axis.line.y = element_line(), axis.ticks.y = element_line(), 
              axis.text.y = element_text(size = 12)) + labs(x = "Time (Ma)")
      if (show_tip_labels) {
        plot_data <- ggtree_object$data %>% filter(.data$isTip == 
                                                     TRUE)
        ggtree_object <- ggtree_object + scale_y_continuous(expand = c(add = c(0.02, 
                                                                               0.02)), position = "right", breaks = plot_data$y, 
                                                            labels = plot_data$label) + theme(legend.position = "bottom", 
                                                                                              axis.ticks.x = element_blank(), axis.line.x = element_blank(), 
                                                                                              axis.text.x = element_text(angle = 90, hjust = 0, 
                                                                                                                         size = 3))
      }
    }
    if (tree_direction == "down") {
      ggtree_object <- ggtree_object + coord_flip() + 
        scale_x_continuous(trans = "reverse", expand = c(add = c(0, 
                                                                 0.01)), breaks = pretty(c(0, agemax), n = time_axis_ticks), 
                           labels = rev(pretty(c(0, agemax), n = time_axis_ticks))) + 
        theme(axis.line.y = element_line(), axis.ticks.y = element_line(), 
              axis.text.y = element_text(size = 12)) + labs(x = "Time (Ma)")
      if (show_tip_labels) {
        plot_data <- ggtree_object$data %>% filter(.data$isTip == 
                                                     TRUE)
        ggtree_object <- ggtree_object + scale_y_continuous(expand = c(add = c(0.02, 
                                                                               0.02)), position = "left", breaks = plot_data$y, 
                                                            labels = plot_data$label) + theme(legend.position = "top", 
                                                                                              axis.ticks.x = element_blank(), axis.line.x = element_blank(), 
                                                                                              axis.text.x = element_text(angle = 270, hjust = 0, 
                                                                                                                         size = 3))
      }
    }
    if (tree_direction == "left") {
      ggtree_object <- ggtree_object + scale_x_continuous(trans = "reverse", 
                                                          expand = c(add = c(0, 0.01)), breaks = pretty(c(0, 
                                                                                                          agemax), n = time_axis_ticks), labels = rev(pretty(c(0, 
                                                                                                                                                               agemax), n = time_axis_ticks))) + theme(axis.line.x = element_line(), 
                                                                                                                                                                                                       axis.ticks.x = element_line(), axis.text.x = element_text(size = 12)) + 
        labs(x = "Time (Ma)")
      if (show_tip_labels) {
        plot_data <- ggtree_object$data %>% filter(.data$isTip == 
                                                     TRUE)
        ggtree_object <- ggtree_object + scale_y_continuous(expand = c(add = c(0.02, 
                                                                               0.02)), position = "left", breaks = plot_data$y, 
                                                            labels = plot_data$label) + theme(legend.position = "top", 
                                                                                              axis.ticks.y = element_blank(), axis.line.y = element_blank(), 
                                                                                              axis.text.y = element_text(hjust = 1, size = 3))
      }
    }
    if (tree_direction == "right") {
      ggtree_object <- ggtree_object + scale_x_continuous(expand = c(add = c(0, 
                                                                             0.01)), breaks = pretty(c(0, agemax), n = time_axis_ticks), 
                                                          labels = rev(pretty(c(0, agemax), n = time_axis_ticks))) + 
        theme(axis.line.x = element_line(), axis.ticks.x = element_line(), 
              axis.text.x = element_text(size = 12)) + labs(x = "Time (Ma)")
      if (show_tip_labels) {
        plot_data <- ggtree_object$data %>% filter(.data$isTip == 
                                                     TRUE)
        ggtree_object <- ggtree_object + scale_y_continuous(expand = c(add = c(0.02, 
                                                                               0.02)), position = "right", breaks = plot_data$y, 
                                                            labels = plot_data$label) + theme(legend.position = "top", 
                                                                                              axis.ticks.y = element_blank(), axis.line.y = element_blank(), 
                                                                                              axis.text.y = element_text(hjust = 0, size = 3))
      }
    }
  }
  if (tree_layout %in% c("circular", "fan", "radial")) {
    maxx <- ggtree_object$data %>% top_n(n = 1, wt = .data$x) %>% 
      dplyr::select(.data$x) %>% unlist %>% unname %>% unique
    maxx <- round(maxx, 1)
    ntip <- ggtree_object$data %>% filter(.data$isTip == 
                                            TRUE) %>% nrow()
    ntip <- ntip + 10
    pretty_points <- maxx - c(maxx, pretty(c(maxx:0), n = time_axis_ticks))
    pp <- tibble(x = rev(pretty_points), y = 0) %>% filter(.data$x <= 
                                                             maxx, .data$x > 0) %>% dplyr::mutate(label = rev(.data$x) - 
                                                                                             min(.data$x))
    ggtree_object <- ggtree_object + geom_vline(data = pp, 
                                                aes(xintercept = .data$x), size = 0.2, color = "darkgrey") + 
      geom_text(data = pp, aes(x = .data$x + 0.1, y = ntip + 
                                 2, label = .data$label), size = 4, inherit.aes = FALSE)
    if (show_tip_labels) {
      ggtree_object <- ggtree_object + geom_tiplab2(inherit.aes = FALSE, 
                                                    size = tip_size, offset = offset)
    }
  }
  return(ggtree_object)
}



m_rate_recon_direction <- function (processed_recon, show_tip_labels = FALSE, parameter = "turnover", 
          discrete = FALSE, breaks = seq(0.3, 0.6, 1), colors = c("red", 
                                                                  "blue", "orange", "green"), plot_as_waiting_time = FALSE, 
          tree_layout = "rectangular", tree_direction = "right", time_axis_ticks = 10, 
          open_angle = 10, right = T, size = 1, tip_size = 1, offset = 0) 
{
  if (!tree_layout %in% c("rectangular", "circular", "slanted", 
                          "fan", "radial")) {
    stop("The dplyr::selected tree layout is not supported.")
  }
  tree <- processed_recon$tree_data@phylo
  datas <- processed_recon$tree_data@data
  agemax <- tree %>% branching.times() %>% max()
  if (plot_as_waiting_time) {
    datas <- mutate(datas, wanted = 1/!!as.name(parameter))
  }
  else {
    datas <- dplyr::mutate(datas, wanted = !!as.name(parameter))
  }
  ggg <- ggtree(
    tr = tree,
    layout = tree_layout,
    size = size,
    open.angle = open_angle,
    ladderize = TRUE,
    right = right
  ) + theme(
    legend.position = "right",
    legend.key.size = unit(x = 0.5, units = "cm"),
    legend.margin = margin(0,
                           0, 0, 0), legend.background = element_blank())
  if (discrete) {
    max_rate <- datas %>% dplyr::select(.data$wanted) %>% unlist %>% 
      unname %>% max
    if (!0 %in% breaks) {
      breaks <- c(0, breaks)
    }
    if (all(max_rate > breaks)) {
      breaks <- c(breaks, max_rate)
    }
    message("Cutting distribution of rate with these breaks:\n")
    print(breaks)
    param <- datas %>% dplyr::select(.data$wanted) %>% unlist %>% 
      unname %>% cut(breaks = breaks)
    ggg <- ggg + aes(color = param) + scale_color_manual(values = colors, name = parameter) + guides(color = guide_legend(override.aes = list(size = 4)))
  }
  else {
    param <- datas %>% dplyr::select(.data$wanted) %>% unlist %>% 
      unname
    ggg <- ggg + aes(color = param) + scale_color_gradient(low = colors[1], 
                                                           high = colors[2], name = parameter)
  }
  ggg <- tree_flip(ggtree_object = ggg, show_tip_labels = show_tip_labels, 
                   tree_layout = tree_layout, tree_direction = tree_direction, 
                   time_axis_ticks = time_axis_ticks, agemax = agemax, tip_size = tip_size, offset = offset)
  return(ggg + theme(plot.margin = unit(rep(0.1, 4), "in")))
}








m_trait_recon_direction <-
  function (processed_recon,
            show_tip_labels = FALSE,
            cutoff = c(0.2, 0.2),
            states_of_first_character,
            states_of_second_character,
            tree_layout = "rectangular",
            tree_direction = "right",
            time_axis_ticks = 10,
            open_angle = 10,
            colors = viridis(n = 9),
            right = T, 
            size = 1, 
            tip_size = 1, 
            offset = 0) {
    if (!tree_layout %in% c("rectangular", "circular", "slanted", 
                          "fan", "radial")) {
    stop("The dplyr::selected tree layout is not supported.")
  }
  tree <- processed_recon$tree_data@phylo
  agemax <- tree %>% branching.times() %>% max()
  ss <- processed_recon$tree_data@data %>% dplyr::mutate(prob_0x_named = case_when(prob_0x >= 
                                                                              1 - cutoff[1] ~ states_of_first_character[1], prob_0x <= 
                                                                              cutoff[1] ~ states_of_first_character[2], TRUE ~ paste(states_of_first_character[1], 
                                                                                                                                     "/", states_of_first_character[2], " uncertain", sep = ""))) %>% 
    dplyr::mutate(prob_x0_named = case_when(prob_x0 >= 1 - cutoff[2] ~ states_of_second_character[1], prob_x0 <= cutoff[2] ~ states_of_second_character[2], TRUE ~ paste(states_of_second_character[1], "/", states_of_second_character[2], " uncertain", sep = ""))) 
  
 ss[,15] <- paste(ss$prob_0x_named, ss$prob_x0_named, sep = "-")
 colnames(ss)[15] <- "four_state"
  
  message("Categories after discretizing with the provided cutoff:\n")
  ss.cnt <- ss %>% ungroup %>% group_by(.data$four_state) %>% 
    add_tally() %>% dplyr::select(.data$four_state, n) %>% distinct
  print(ss.cnt)
  nstat <- ss %>% dplyr::select(.data$four_state) %>% distinct %>% 
    nrow
  ggg <- ggtree(tr = tree, layout = tree_layout, size = size, ladderize = TRUE,
                right = right, open.angle = open_angle, aes(color = ss$four_state)) + 
    scale_color_manual(name = "", values = colors[1:nstat]) + 
    theme(legend.position = "right", legend.key.size = unit(x = 0.5, 
                                                            units = "cm"), legend.margin = margin(0, 0, 0, 0), 
          legend.background = element_blank(), plot.background = element_blank()) + 
    guides(color = guide_legend(override.aes = list(size = 4)))
  ggg <- tree_flip(ggtree_object = ggg, show_tip_labels = show_tip_labels, 
                   tree_layout = tree_layout, tree_direction = tree_direction, 
                   time_axis_ticks = time_axis_ticks, agemax = agemax, tip_size = tip_size, offset = offset)
  return(ggg + theme(plot.margin = unit(rep(0.1, 4), "in")))
}







m_ridgelines_adjust <- function(processed_recon, states_names = c("00", "01", "10", "11"), parameter = "turnover", plot_as_waiting_time = FALSE, fill_colors = rep(NA, 4), line_colors = viridis(n = 4), state_order = c("00", "01", "10", "11"), plot_min = F) {
  message("Recoding and renaming character states. The elements 1:4 of the vector `character_states_names` are assumed to match the states 00, 01, 10, 11.\n")
  ss <- utilhisse:::m_prep_df(processed_recon = processed_recon, states_names = states_names, parameter = parameter)
  message("Summarising grouped by character state\n")
  wanted <- as.name(parameter)
  ss_var <- ss %>% group_by(.data$four_state) %>% dplyr::select(.data$four_state, 
                                                         .data$wanted) %>% distinct()
  if (nrow(ss_var) == length(states_names)) {
    print(ss_var)
    stop("Looks like there is no variation within the observed states. Probably because the model has no hidden states (e.g. MuSSE?). No point in calculating densities for point estimates.")
  }
  if (plot_as_waiting_time) {
    ss <- dplyr::mutate(ss, wanted = 1/!!wanted)
    ss.sum <- ss %>% group_by(.data$four_state) %>% utilhisse:::summ()
  }
  else {
    ss <- dplyr::mutate(ss, wanted = !!wanted)
    ss.sum <- ss %>% group_by(.data$four_state) %>% utilhisse:::summ()
  }
  max_rate <- ss %>% dplyr::select(.data$wanted) %>% top_n(1, wt = .data$wanted) %>% unlist %>% unname %>% unique
  
  min_rate <- ss %>% dplyr::select(.data$wanted) %>% top_n(-1, wt = .data$wanted) %>% unlist %>% unname %>% unique
  print(ss.sum)
  
  message("\nPlotting\n\n")
  
  if (plot_min == F){
    ss$four_state <- factor(ss$four_state, levels = state_order)
    ggplot(data = ss, aes(x = .data$wanted, y = .data$four_state, 
                          fill = .data$four_state, colour = .data$four_state)) + 
      ggridges::geom_density_ridges(alpha = 0.75, size = 0.75, jittered_points = TRUE, 
                                    scale = 0.83, rel_min_height = 0.01, point_shape = "|", 
                                    point_size = 1, position = position_nudge(y = rep(-0.2, 4))) + geom_errorbarh(data = ss.sum, position = position_nudge(y = rep(-0.3, 4)), aes(xmin = .data$Mean - .data$SD, xmax = .data$Mean + .data$SD, y = .data$four_state, colour = .data$four_state), height = 0.05, inherit.aes = FALSE) + geom_point(data = ss.sum, pch = 21, position = position_nudge(y = rep(-0.3, 4)), aes(y = .data$four_state, x = .data$Mean, colour = .data$four_state), size = 3, inherit.aes = FALSE) + scale_x_continuous(breaks = pretty(x = c(min_rate, max_rate), n = 10)) + scale_fill_manual(name = "", values = fill_colors) + scale_colour_manual(name = "", values = line_colors) + labs(y = "", x = parameter) + theme_classic() + theme(legend.position = "none") + utilhisse:::theme_xy
  
  } else {
    ss$four_state <- factor(ss$four_state, levels = state_order)
    ggplot(data = ss, aes(x = .data$wanted, y = .data$four_state, 
                          fill = .data$four_state, colour = .data$four_state)) + 
      ggridges::geom_density_ridges(alpha = 0.75, size = 0.75, jittered_points = TRUE, 
                                    scale = 0.83, rel_min_height = 0.01, point_shape = "|", 
                                    point_size = 1, position = position_nudge(y = rep(-0.2, 
                                                                                      4))) + geom_errorbarh(data = ss.sum, position = position_nudge(y = rep(-0.3, 
                                                                                                                                                             4)), aes(xmin = .data$Mean - .data$SD, xmax = .data$Mean + .data$SD, y = .data$four_state, colour = .data$four_state), 
                                                                                                height = 0.05, inherit.aes = FALSE) + geom_point(data = ss.sum, 
                                                                                                                                                 pch = 21, position = position_nudge(y = rep(-0.3, 4)), 
                                                                                                                                                 aes(y = .data$four_state, x = .data$Mean, colour = .data$four_state), 
                                                                                                                                                 size = 3, inherit.aes = FALSE) + scale_x_continuous(breaks = pretty(x = c(min_rate, 
                                                                                                                                                                                                                           max_rate), n = 10)) + scale_fill_manual(name = "", values = fill_colors) + 
    scale_colour_manual(name = "", values = line_colors) + 
    labs(y = "", x = parameter) + theme_classic() + theme(legend.position = "none") + 
    utilhisse:::theme_xy
  }
}







