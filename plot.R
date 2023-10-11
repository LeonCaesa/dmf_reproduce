# [plotting utils]
t <- list(
  family = "Arial",
  size = 12)

addSmallLegend <- function(myPlot, pointSize = 0.5, textSize = 3, spaceLegend = 0.1) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize),
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}
remove_geom <- function(ggplot2_object, geom_type) {
  # Delete layers that match the requested type.
  layers <- lapply(ggplot2_object$layers, function(x) {
    if (class(x$geom)[1] == geom_type) {
      NULL
    } else {
      x
    }})
  # Delete the unwanted layers.
  layers <- layers[!sapply(layers, is.null)]
  ggplot2_object$layers <- layers
  ggplot2_object}

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
scaleFUN <- function(x) sprintf("%.2f", x)

pc_plotUpdate<-function(label, fit_result, plot_title){
  label_string = label
  label = as.numeric(factor(label))
  Rank_temp = dim(fit_result$V)[2]
  col_names = c('PC1', 'PC2', 'PC3', 'PC4')[1:Rank_temp]

  # [identify]
  dmf_identified = dmf_identify(tcrossprod(fit_result$L, fit_result$V), dim(fit_result$L)[2])
  fit_result$L = dmf_identified$L
  fit_result$V = dmf_identified$V

  principle_result = data.frame(fit_result$V[,1:Rank_temp])
  colnames(principle_result) = col_names

  principle_result = principle_result[order(label),]
  sep_list = rep(0, length(unique(label))-1)
  sep_list[1] = sum(label==1)
  for (j in 2:length(unique(label))){
    sep_list[j] = sep_list[j-1] + sum(label==j)}

  if(dim(fit_result$V)[2] ==1){
    principle_result = data.frame(cbind( index =1:length(principle_result),  value = principle_result, variable = label))
    return(ggplot(principle_result) +
             ggtitle(plot_title) +
             theme_bw(base_size = 5)+
             geom_point( aes(x= index, y = value, colour = variable)) +
             geom_vline(xintercept= sep_list , linetype = 'dashed' ) +
             annotate("text", x = sep_list, y =  mean(principle_result, na.rm = TRUE),
                      angle = -90, label = unique(label_string),
                      vjust = 1.2) + scale_y_continuous(labels=scaleFUN) +
             xlab('Patient Index') + ylab('PC Score')+ theme_light() +
             theme( plot.title = element_text(hjust = 0.5), legend.position="bottom")+
             guides(colour=guide_legend(title="PC Index"))
    )}

  long_pc = melt(principle_result)
  long_pc$index = rep(1:dim(principle_result)[1], dim(principle_result)[2])
  return(ggplot(long_pc) +
           ggtitle(plot_title) +
           theme_bw(base_size = 8)+
           geom_point( aes(x= index, y = value, colour = variable)) +
           geom_vline(xintercept= sep_list, linetype = 'dashed' ) +
           annotate("text", x = sep_list-5, y = 0.01, angle = 0, label = unique(label_string),
                    vjust = 1.2, size = 3) + scale_y_continuous(labels=scaleFUN) +
           xlab('Patient Index') + ylab('Principle Score')+ theme_light()+
           theme( plot.title = element_text(hjust = 0.5),legend.position="bottom") +
           guides(colour=guide_legend(title="PC Index"))
  )
}


plot_3dUpdate = function(pc_result, xrange, yrange, zrange, opacity = 1, font_size =12){
  return(plot_ly(pc_result, x=~X1, y=~X2, z=~X3, opacity =opacity,
                 type="scatter3d", mode="markers", color=~label,
                 marker = list(size = 4, alpha = 0.5)) %>%
           layout(font=t,
          scene = list(xaxis = list(title = 'PC1',
                                  range = xrange,
                                  titlefont = list(size = font_size-8)
                                  ),
                     yaxis = list(title = 'PC2',
                                  range = yrange,
                                  titlefont = list(size = font_size-8)),
                     zaxis = list(title = 'PC3',
                                  range = zrange,
                                  titlefont = list(size = font_size-8))),
                  margin = list(
                    l = 0,r = 0,b = 0,t = 0,pad = 1),
                legend = list(orientation = "h",   # show entries horizontally
                           xanchor = "center",  # use center of legend as anchor
                           x = 0.5,
                           font = list(size = font_size)
                           )))}





GOF_test <- function(G, glm_family, fit_result, Y, glm_weights, phi_hat = 1, phi_estimate = FALSE, return_fit =FALSE){
  if(glm_family$family=='binomial'){
    Y = Y/glm_weights}

  if((glm_family$family== "Gamma")&(phi_estimate ==TRUE)){
    warning("An MLE estimate of gamma 1/v is used")
    phi_hat = 1/mean((Y-fit_result$mu_hat)^2)}
  if((glm_family$family== "gaussian")&(phi_estimate ==TRUE)){
    warning("An MLE estimate of gaussian sigma^2 is used")
    phi_hat = mean((Y-fit_result$mu_hat)^2)}

  n = dim(Y)[1]
  d = dim(Y)[2]

  eta_hat = data.frame(eta = c(glm_family$linkfun(fit_result$mu_hat)))
  cut_quantle = seq(0, 1, 1/G)
  eta_g = eta_hat %>% group_by() %>% mutate(G = cut(eta, quantile(eta, cut_quantle),
                                                    labels = 1:G,
                                                    include.lowest = TRUE))


  resids = matrix(Y-fit_result$mu_hat)/sqrt(n*d)
  variances = matrix(glm_family$variance(fit_result$mu_hat)/ glm_weights)/(n*d)
  s_vec = matrix(0, nrow = G)
  d_diag = rep(0, G)
  for (g in 1:G){
    s_vec[g] = sum(resids[eta_g$G==g])
    d_diag[g] = phi_hat * sum(variances[eta_g$G==g])
  }
  qqplot(s_vec/sqrt(d_diag), rnorm(length(s_vec)))

  if(return_fit ==TRUE){
    return(s_vec/sqrt(d_diag))

  }

  GOF = crossprod(t(crossprod(s_vec, solve(diag(d_diag)))), s_vec)[1]
  return(GOF)
}


# Remarks:
#
# The export is done using the automated testing framework [Selenium](https://
# de.wikipedia.org/wiki/Selenium) which results in opening a browser window
# (Google Chrome) that might has to be closed by hand. Other than Plotly's
# own `export()` function this one also allows to set the `width` and `height`
# of the exported plot (in the former it's hardcoded to 800x600 pixels). If
# `incl_PDF_copy`/`incl_PNG_copy` is set to `TRUE`, the exported SVG additionally
# gets converted to a PDF/PNG using the R package [`rsvg`](https://github.com/
# jeroen/rsvg/tree/40576ac326621b40224db344b09158f4ff717433) which relies on
# [`librsvg`](https://de.wikipedia.org/wiki/Librsvg). On Linux distributions
# the development package of `librsvg` must be installed. On macOS the required
# dependency (`librsvg`) can be installed using [Homebrew](https://brew.sh/).
# Optional PNG auto-cropping is done using the `imager` R package.

ensure_package <- Vectorize(
  FUN =
    function(package,
             load = TRUE)
    {
      installed_packages <- rownames(installed.packages())

      if ( !(package %in% installed_packages) )
      {
        install.packages(package,
                         repos = "https://cloud.r-project.org/")
      }
      if ( load ) library(package, character.only = TRUE)
    }
)

  remove_geom <- function(ggplot2_object, geom_type) {
    # Delete layers that match the requested type.
    layers <- lapply(ggplot2_object$layers, function(x) {
      if (class(x$geom)[1] == geom_type) {
        NULL
      } else {
        x
      }})
    # Delete the unwanted layers.
    layers <- layers[!sapply(layers, is.null)]
    ggplot2_object$layers <- layers
    ggplot2_object}
export_plotly2SVG <- function(plotly_graph,
                              filename = NULL,
                              parent_path = paste0(getwd(), "/output"),
                              width = 800,
                              height = 600,
                              remove_title = FALSE,
                              font_family = "Arial",
                              incl_PDF_copy = FALSE,
                              incl_PNG_copy = FALSE,
                              png_scaling_factor = 1.8,
                              autocrop_png = TRUE)
{
  ensure_package(package = c("dplyr",
                             "plotly",
                             "readr",
                             "RSelenium",
                             "rsvg",
                             "stringr"),
                 load = FALSE)

  ensure_package("magrittr")

  # remove trailing slash in `parent_path`
  parent_path %<>% normalizePath()

  # ensure `parent_path` exists
  if ( !dir.exists(parent_path) ) dir.create(path = parent_path,
                                             recursive = TRUE)

  # generate sensible filename
  if ( is.null(filename) )
  {
    auto_name <- deparse(substitute(plotly_graph))

    filename <- dplyr::if_else(
      condition = auto_name == ".",
      true = "plotly_graph.svg",
      false = paste0(deparse(substitute(plotly_graph)), ".svg")
    )
  } else
  {
    filename %<>%
      stringr::str_replace(pattern = "([^\\.svg])$",
                           replacement = "\\1.svg")
  }

  filepath <- paste0(parent_path, "/", filename)

  # delete old SVG file
  if ( file.exists(filepath) )
  {
    unlink(x = filepath)
  }

  if ( remove_title )
  {
    plotly_graph %<>%
      plotly::layout(title = "",
                     margin = list(t = 0))
  }

  if ( !is.null(font_family) )
  {
    plotly_graph %<>%
      plotly::layout(font = list(family = font_family))
  }

  # temporarily export plot to a HTML file
  tempfile <- tempfile(pattern = "plotly_temp_",
                       tmpdir = parent_path,
                       fileext = ".html")

  export_plotly2HTML(plotly_graph = plotly_graph,
                     filename = basename(tempfile),
                     parent_path = parent_path)

  on.exit(unlink(tempfile),
          add = TRUE)

  # get <div> ID of exported htmlwidget
  htmlwidget_id <-
    stringr::str_extract(string = readr::read_file(file = tempfile),
                         pattern = "(?<=<div id=\")htmlwidget-[^\"]+")

  # initialize Chrome as RSelenium driver
  selenium_driver <-
    RSelenium::rsDriver(browser = "chrome",
                        extraCapabilities = list(
                          chromeOptions = list(
                            prefs = list(
                              "profile.default_content_settings.popups" = 0L,
                              "download.prompt_for_download" = FALSE,
                              "download.default_directory" = parent_path
                            )
                          )
                        ),
                        verbose = FALSE)

  # navigate to temporary HTML file
  selenium_driver$client$navigate(url = paste0("file://", normalizePath(tempfile)))

  # download plot as SVG using the native
  # [`Plotly.downloadImage`](https://plot.ly/javascript/plotlyjs-function-reference/#plotlydownloadimage) function
  selenium_driver$client$executeScript(
    script = paste0("Plotly.downloadImage(document.getElementById('", htmlwidget_id, "'), ",
                    "{format: 'svg', width: ", width, ", height: ", height, ", filename: '",
                    tools::file_path_sans_ext(x = filename), "'});"),
    args = list(NULL)
  )

  # wait for SVG to be saved to disk
  Sys.sleep(time = 1)

  # convert to PNG
  if ( incl_PNG_copy )
  {
    filepath_png <- paste0(tools::file_path_sans_ext(parent_path), ".png")

    rsvg::rsvg_png(svg = filepath,
                   file = filepath_png,
                   width = png_scaling_factor * width,
                   height = png_scaling_factor * height)

    if ( autocrop_png ) autocrop_png(path_to_png = filepath_png)
  }

  # convert to PDF
  if ( incl_PDF_copy )
  {
    rsvg::rsvg_pdf(svg = filepath,
                   file = paste0(tools::file_path_sans_ext(parent_path), ".pdf"))
  }


}

export_plotly2HTML <- function(plotly_graph,
                               filename = NULL,
                               parent_path = paste0(getwd(), "/output"),
                               selfcontained = FALSE,
                               libdir = "plotly_files",
                               disable_legend_toggling = NULL,
                               # you can provide a link to webfonts to be used like this:
                               # add_web_font = "https://fonts.googleapis.com/css?family=Work+Sans:200,300,400,600,700"
                               add_web_font = NULL)
{
  ensure_package(package = c("checkmate",
                             "dplyr",
                             "readr",
                             "stringr"),
                 load = FALSE)

  ensure_package("magrittr")

  # remove trailing slash in `parent_path`
  parent_path %<>% normalizePath()

  # ensure `parent_path` exists
  if ( !dir.exists(parent_path) ) dir.create(path = parent_path,
                                             recursive = TRUE)

  # generate sensible filename
  if ( is.null(filename) )
  {
    auto_name <- deparse(substitute(plotly_graph))

    filename <- dplyr::if_else(
      condition = auto_name == ".",
      true = "plotly_graph.html",
      false = paste0(deparse(substitute(plotly_graph)), ".html")
    )
  }

  filepath <- paste0(parent_path, "/", filename)

  htmlwidgets::saveWidget(
    widget = plotly_graph,
    file = filepath,
    selfcontained = selfcontained,
    libdir = libdir
  )

  if ( !is.null(disable_legend_toggling) )
  {
    test_char <- checkmate::check_choice(x = disable_legend_toggling,
                                         choices = "all")
    test_num <- checkmate::check_numeric(x = disable_legend_toggling,
                                         lower = 1,
                                         upper = length(plotly_graph$x$attrs),
                                         min.len = 1,
                                         max.len = length(plotly_graph$x$attrs),
                                         unique = TRUE,
                                         any.missing = FALSE,
                                         all.missing = FALSE)

    if ( !isTRUE(test_char) & !isTRUE(test_num) )
    {
      stop("Invalid argument provided: disable_legend_toggling\n",
           ifelse(!isTRUE(test_char) & is.character(disable_legend_toggling),
                  paste0(test_char, ". Or alternatively can also be a vector of integers >= 1 and <= number of traces."),
                  ""),
           ifelse(!isTRUE(test_num) & is.numeric(disable_legend_toggling),
                  paste0(test_num, ". Or alternatively can also be \"all\"."),
                  ""))

    } else if ( isTRUE(test_char) )
    {
      css_rules <-
        c("",
          "/* hides the svg dom element that has the click handler responsible for toggling */",
          ".legend .traces .legendtoggle {",
          "  display: none;",
          "}",
          "/* just for presentation: shows the default cursor instead of the text cursor */",
          ".legend .traces .legendtext {",
          "  cursor: default;",
          "}",
          "")
    } else
    {
      disable_legend_toggling %<>% as.integer()

      css_rules <-
        c("",
          "/* hides the svg dom element that has the click handler responsible for toggling */")

      for ( i in disable_legend_toggling )
      {
        css_rules %<>%
          c(paste0(".legend .groups:nth-of-type(", i, ") .legendtoggle",
                   dplyr::if_else(i == last(disable_legend_toggling),
                                  " {",
                                  ","), " "))
      }

      css_rules %<>%
        c("  display: none;",
          "}",
          "/* just for presentation: shows the default cursor instead of the text cursor */")

      for ( i in disable_legend_toggling )
      {
        css_rules %<>%
          c(paste0(".legend .groups:nth-of-type(", i, ") .legendtext",
                   dplyr::if_else(i == last(disable_legend_toggling),
                                  " {",
                                  ","), " "))
      }

      css_rules %<>%
        c("  cursor: default;",
          "}",
          "")
    }

    # write modified .css file
    plotly_dir <-
      list.dirs(path = paste0(parent_path, "/", libdir),
                full.names = TRUE,
                recursive = FALSE) %>%
      stringr::str_subset(pattern = "plotlyjs")

    readr::read_lines(file = paste0(plotly_dir, "/plotly-htmlwidgets.css")) %>%
      c(css_rules) %>%
      readr::write_lines(path = paste0(plotly_dir, "/plotly_htmlwidgets_custom.css"),
                         append = FALSE)

    # modify dependency path in HTML file
    readr::read_file(file = filepath) %>%
      stringr::str_replace(pattern = "plotly-htmlwidgets\\.css",
                           replacement = "plotly_htmlwidgets_custom.css") %>%
      readr::write_file(path = filepath,
                        append = FALSE)
  }

  if ( !is.null(add_web_font) )
  {
    webfont_tag <-
      "<link href=\"" %>%
      paste0(checkmate::assert_character(x = add_web_font,
                                         pattern = "^https?://\\w.*",
                                         ignore.case = TRUE,
                                         any.missing = FALSE,
                                         all.missing = FALSE,
                                         unique = TRUE)) %>%
      paste0("\" rel=\"stylesheet\" />")

    readr::read_file(file = filepath) %>%
      stringr::str_replace(pattern = "<link href=",
                           replacement = paste0(webfont_tag, "\n<link href=")) %>%
      readr::write_file(path = filepath,
                        append = FALSE)
  }
}

autocrop_png <- function(path_to_png)
{
  ensure_package("magrittr")
  ensure_package(package = "imager",
                 load = FALSE)

  imager::load.image(file = path_to_png) %>%
    imager::autocrop() %>%
    imager::pad(nPix = 4,
                axes = "xy",
                pos = 0) %>%
    imager::flatten.alpha() %>%
    imager::save.image(file = path_to_png)
}
