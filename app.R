# STA2ASM Week 4 — Continuous distributions and the Central Limit Theorem
# Standalone Shiny app; requires only the shiny package.
library(shiny)

PAL <- list(teal = "#236A6E", blue = "#174F8A", coral = "#C0504D",
            green = "#548235", purple = "#7030A0", lightblue = "#9DC3E6",
            lightcoral = "#F4A6A6", lightgreen = "#A9D8B2",
            lightpurple = "#C5B5E8", binomial = "#FFFF00",
            poisson = "#C65D7B", negbin = "#8E70B5")
fmt <- function(x, digits = 4) formatC(x, format = "f", digits = digits)
math <- function(s) withMathJax(HTML(paste0("\\(", s, "\\)")))

# One specification for all distribution-driven panels. The traditional NB
# variable counts total trials up to the rth success (not R's failure count).
make_spec <- function(dist, mu = 0, sd = 1, df = 5, trials = 20,
                      prob = .3, rate = 2.6, size = 3, nbmu = 7, theta = 3) {
  switch(dist,
    norm = list(label = "Normal", discrete = FALSE, mean = mu, var = sd^2,
      dens = function(x) dnorm(x, mu, sd), cdf = function(x) pnorm(x, mu, sd),
      quant = function(p) qnorm(p, mu, sd), gen = function(n) rnorm(n, mu, sd),
      latex = sprintf("X \\sim N(%s,%s^2)", fmt(mu, 2), fmt(sd, 2)),
      code = sprintf("mean = %s, sd = %s", fmt(mu, 2), fmt(sd, 2)),
      col = PAL$blue, fill = PAL$lightblue),
    t = list(label = "t", discrete = FALSE, mean = if (df > 1) 0 else NA_real_,
      var = if (df > 2) df/(df-2) else Inf,
      dens = function(x) dt(x, df), cdf = function(x) pt(x, df),
      quant = function(p) qt(p, df), gen = function(n) rt(n, df),
      latex = sprintf("X \\sim t_{%d}", df), code = sprintf("df = %d", df),
      col = PAL$purple, fill = PAL$lightpurple),
    chisq = list(label = "Chi-squared", discrete = FALSE, mean = df, var = 2*df,
      dens = function(x) dchisq(x, df), cdf = function(x) pchisq(x, df),
      quant = function(p) qchisq(p, df), gen = function(n) rchisq(n, df),
      latex = sprintf("X \\sim \\chi^2_{%d}", df), code = sprintf("df = %d", df),
      col = PAL$coral, fill = PAL$lightcoral),
    lnorm = list(label = "Log-normal", discrete = FALSE,
      mean = exp(mu + sd^2/2), var = (exp(sd^2)-1)*exp(2*mu+sd^2),
      dens = function(x) dlnorm(x, mu, sd), cdf = function(x) plnorm(x, mu, sd),
      quant = function(p) qlnorm(p, mu, sd), gen = function(n) rlnorm(n, mu, sd),
      latex = sprintf("\\log(X) \\sim N(%s,%s^2)", fmt(mu, 2), fmt(sd, 2)),
      code = sprintf("meanlog = %s, sdlog = %s", fmt(mu, 2), fmt(sd, 2)),
      col = PAL$green, fill = PAL$lightgreen),
    bern = list(label = "Bernoulli", discrete = TRUE, mean = prob, var = prob*(1-prob),
      dens = function(x) dbinom(x, 1, prob), cdf = function(x) pbinom(x, 1, prob),
      quant = function(p) qbinom(p, 1, prob), gen = function(n) rbinom(n, 1, prob),
      latex = sprintf("X \\sim \\operatorname{Bernoulli}(%.2f)", prob),
      code = sprintf("size = 1, prob = %s", fmt(prob, 2)),
      col = "#D6A300", fill = PAL$binomial),
    bin = list(label = "Binomial", discrete = TRUE, mean = trials*prob,
      var = trials*prob*(1-prob),
      dens = function(x) dbinom(x, trials, prob), cdf = function(x) pbinom(x, trials, prob),
      quant = function(p) qbinom(p, trials, prob), gen = function(n) rbinom(n, trials, prob),
      latex = sprintf("X \\sim \\operatorname{Binomial}(%d,%.2f)", trials, prob),
      code = sprintf("size = %d, prob = %s", trials, fmt(prob, 2)),
      col = "#D6A300", fill = PAL$binomial),
    pois = list(label = "Poisson", discrete = TRUE, mean = rate, var = rate,
      dens = function(x) dpois(x, rate), cdf = function(x) ppois(x, rate),
      quant = function(p) qpois(p, rate), gen = function(n) rpois(n, rate),
      latex = sprintf("X \\sim \\operatorname{Poisson}(%.2f)", rate),
      code = sprintf("lambda = %s", fmt(rate, 2)),
      col = PAL$poisson, fill = PAL$poisson),
    nbt = list(label = "Negative binomial (total trials)", discrete = TRUE,
      mean = size/prob, var = size*(1-prob)/prob^2,
      dens = function(x) dnbinom(x-size, size, prob),
      cdf = function(x) pnbinom(x-size, size, prob),
      quant = function(p) qnbinom(p, size, prob)+size,
      gen = function(n) rnbinom(n, size, prob)+size,
      latex = sprintf("X=F+%d,\\quad F\\sim\\operatorname{NB}(%d,%.2f)",size,size,prob),
      code = sprintf("size = %d, prob = %s", size, fmt(prob, 2)),
      col = PAL$negbin, fill = PAL$negbin),
    nba = list(label = "Negative binomial (mean–dispersion)", discrete = TRUE,
      mean = nbmu, var = nbmu+nbmu^2/theta,
      dens = function(x) dnbinom(x, size = theta, mu = nbmu),
      cdf = function(x) pnbinom(x, size = theta, mu = nbmu),
      quant = function(p) qnbinom(p, size = theta, mu = nbmu),
      gen = function(n) rnbinom(n, size = theta, mu = nbmu),
      latex = sprintf("F\\sim\\operatorname{NB}(\\mu=%.2f,\\theta=%.2f)",nbmu,theta),
      code = sprintf("size = %s, mu = %s",fmt(theta,2),fmt(nbmu,2)),
      col = PAL$negbin, fill = PAL$negbin)
  )
}

# Use a stable display window for t curves so changes in degrees of freedom
# show differences in shape rather than repeatedly rescaling the horizontal axis.
# Theoretical probabilities and quantiles are always calculated on full support.
view_range <- function(s, lower = .005, upper = .995) {
  if (identical(s$label, "t")) return(c(-5, 5))
  # A very large log-scale SD produces a long tail that can conceal the
  # visible body of the density. Focus the display on the main density.
  if (identical(s$label, "Log-normal")) {
    median_x <- s$quant(.5)
    sdlog <- log(s$quant(.841344746) / median_x)
    mode_x <- median_x * exp(-sdlog^2)
    return(c(0, min(s$quant(.995), 100 * mode_x)))
  }
  q <- suppressWarnings(s$quant(c(lower, upper)))
  if (any(!is.finite(q)) || diff(q) <= 0) q <- c(-4, 4)
  if (s$discrete) c(floor(q[1]), ceiling(q[2])) else q
}

# Density of chi-squared(1) is unbounded at zero; cap the display height
# without altering the distribution or its probability calculations.
density_plot_limit <- function(s, x, y) {
  if (identical(s$label, "Chi-squared") && is.finite(s$var) && s$mean < 2)
    return(2)
  valid <- y[is.finite(y)]
  if (!length(valid)) return(1)
  max(valid) * 1.14
}

# Preserve a useful reference scale while parameters change. Expand the window
# only when the current distribution would otherwise be substantially clipped.
# Reset axes explicitly to fit the current distribution again.
axis_limits <- function(s, baseline = NULL, marks = numeric(0), cdf = FALSE) {
  current <- view_range(s)
  if (is.null(baseline)) baseline <- current
  if (identical(s$label, "t")) {
    lim <- c(-5, 5)
  } else {
    lim <- baseline
    span <- diff(lim)
    if (!is.finite(span) || span <= 0) span <- 1
    # Only enlarge when the central 99% is not already visible.
    if (current[1] < lim[1]) lim[1] <- current[1] - .06 * span
    if (current[2] > lim[2]) lim[2] <- current[2] + .06 * span
    visible_marks <- marks[is.finite(marks)]
    if (length(visible_marks)) {
      if (min(visible_marks) < lim[1]) lim[1] <- min(visible_marks) - .06 * span
      if (max(visible_marks) > lim[2]) lim[2] <- max(visible_marks) + .06 * span
    }
  }
  if (diff(lim) <= 0) lim <- lim + c(-1, 1)
  xx <- seq(lim[1], lim[2], length.out = 1200)
  if (identical(s$label, "Chi-squared")) xx <- pmax(xx, 1e-3)
  yy <- s$dens(xx)
  peak <- density_plot_limit(s, xx, yy)
  base_peak <- if (is.null(attr(baseline, "peak"))) peak else attr(baseline, "peak")
  # Avoid y-axis shrinkage; still allow expansion when a curve becomes taller.
  ymax <- if (cdf) 1.02 else max(peak, base_peak)
  list(x = lim, y = c(0, ymax))
}

plot_density <- function(s, mode = "none", a = 0, b = 1, p = .975,
                         main = "Probability density", cdf = FALSE,
                         limits = NULL) {
  if (mode == "quantile") a <- s$quant(p)
  if (mode == "between" && b < a) { tmp <- a; a <- b; b <- tmp }
  marks <- if (mode == "between") c(a,b) else if (mode != "none") a else numeric(0)
  if (is.null(limits)) limits <- axis_limits(s, marks=marks, cdf=cdf)
  lim <- limits$x
  x <- seq(lim[1], lim[2], length.out = 850)
  if (identical(s$label, "Chi-squared")) x <- pmax(x, 1e-3)
  y <- if (cdf) s$cdf(x) else s$dens(x)
  y[!is.finite(y)] <- NA_real_
  plot(x, y, type = "n", xlab = "Value", ylab = if (cdf) "Cumulative probability" else "Density",
       main = main, las = 1, xlim = lim, ylim = limits$y)
  if (!cdf && mode != "none") {
    sel <- switch(mode, left = x <= a, right = x >= a,
                  between = x >= a & x <= b, quantile = x <= a, rep(FALSE,length(x)))
    if (any(sel)) polygon(c(x[sel][1], x[sel], tail(x[sel],1)),
                          c(0, y[sel], 0), col = adjustcolor(s$fill, alpha.f=.63), border=NA)
  }
  lines(x, y, col = s$col, lwd = 3)
  for (v in marks) if (is.finite(v) && v >= lim[1] && v <= lim[2]) {
    abline(v = v, col = PAL$coral, lty = 2, lwd = 1.7)
    if (cdf) points(v, s$cdf(v), pch=21, bg=PAL$lightcoral, col=PAL$coral, cex=1.5)
  }
  if (identical(s$label, "Chi-squared") && !cdf && s$mean == 1)
    mtext("Density near zero is clipped for readability", side=3, line=0.15, cex=.7, col=PAL$teal)
  if (identical(s$label, "Log-normal") && !cdf &&
      lim[2] < s$quant(.995) * .99)
    mtext("Long right tail extends beyond this view", side=3,
          line=0.15, cex=.7, col=PAL$teal)
  if (identical(s$label, "t") && any(marks < lim[1] | marks > lim[2]))
    mtext("Selected cutoff is outside the displayed range", side=3, line=0.15, cex=.7, col=PAL$teal)
}

# Discrete exact probability and continuity-corrected (or uncorrected) normal approximation.
approx_probability <- function(dist, trials, prob, lambda, mode, a, b, cc) {
  mean <- if (dist == "bin") trials*prob else lambda
  sd <- sqrt(if (dist == "bin") trials*prob*(1-prob) else lambda)
  cdf <- if (dist == "bin") function(k) pbinom(k, trials, prob) else function(k) ppois(k, lambda)
  offset <- if (cc) .5 else 0
  k <- floor(a); h <- floor(b)
  if (mode == "left") {
    exact <- cdf(k); approx <- pnorm(k+offset, mean, sd)
    bounds <- c(-Inf,k+offset)
  } else if (mode == "right") {
    exact <- 1-cdf(k); approx <- 1-pnorm(k+offset, mean, sd)
    bounds <- c(k+offset,Inf)
  } else {
    lo <- min(k,h); hi <- max(k,h)
    exact <- cdf(hi)-cdf(lo-1)
    bounds <- c(lo-offset,hi+offset)
    approx <- pnorm(bounds[2],mean,sd)-pnorm(bounds[1],mean,sd)
  }
  list(exact=exact, approx=approx, mean=mean, sd=sd, bounds=bounds)
}

# Render MathJax in control labels; native dropdown options use Unicode.
math_label <- function(text, symbol) {
  withMathJax(HTML(paste0(text, " \\(", symbol, "\\)")))
}

# Distribution-specific parameter controls, reused by Explore, Simulation and CLT.
parameter_ui <- function(prefix, dist) {
  id <- function(s) paste0(prefix,"_",s)
  switch(dist,
    norm = tagList(numericInput(id("mu"), math_label("Mean", "\\mu"), 0, step=.5),
                   numericInput(id("sd"), math_label("Standard deviation", "\\sigma"), 1, min=.1, max=8, step=.1)),
    t = numericInput(id("df"), math_label("Degrees of freedom", "\\nu"), 5, min=1, max=500, step=1),
    chisq = numericInput(id("df"), math_label("Degrees of freedom", "\\nu"), 5, min=1, max=500, step=1),
    lnorm = tagList(numericInput(id("mu"), math_label("Mean on log scale", "\\mu"), 0, step=.25),
                    numericInput(id("sd"), math_label("Standard deviation on log scale", "\\sigma"), 1, min=.1, max=3, step=.1)),
    bern = sliderInput(id("prob"), math_label("Probability of success", "p"), min=.01,max=.99,value=.3,step=.01),
    bin = tagList(numericInput(id("trials"), math_label("Number of trials", "n"), 20,min=1,max=300,step=1),
                  sliderInput(id("prob"), math_label("Probability of success", "p"),min=.01,max=.99,value=.3,step=.01)),
    pois = numericInput(id("rate"), math_label("Expected count", "\\lambda"), 2.6,min=.05,max=100,step=.5),
    nbt = tagList(numericInput(id("size"), math_label("Required successes", "r"), 3,min=1,max=100,step=1),
                  sliderInput(id("prob"), math_label("Probability of success", "p"),min=.01,max=.99,value=.3,step=.01)),
    nba = tagList(numericInput(id("nbmu"), math_label("Mean", "\\mu"), 7,min=.05,max=100,step=.5),
                  numericInput(id("theta"), math_label("Dispersion", "\\theta"), 3,min=.1,max=100,step=.5))
  )
}

cont_choices <- c("Normal"="norm", "t"="t", "Chi-squared"="chisq", "Log-normal"="lnorm")
all_choices <- c("Bernoulli"="bern", "Binomial"="bin", "Poisson"="pois",
                 "Negative binomial (total trials)"="nbt",
                 "Negative binomial (mean–dispersion)"="nba", cont_choices)

ui <- fluidPage(
  tags$head(tags$style(HTML('
    body {font-family:"Helvetica Neue",Arial,sans-serif; color:#333; background:white;}
    h1,h2,h3,h4 {color:#236A6E; font-weight:500;}
    .well {background:#eef5f5; border:1px solid #cfdfe0; border-left:4px solid #4d8386;}
    .navbar-default {background:white;border:0;border-bottom:1px solid #dce6e6;border-radius:0;box-shadow:none;margin-bottom:16px;}
    .navbar-default .navbar-nav>li>a {color:#4a5b5c; padding:12px 14px 9px;border-bottom:3px solid transparent;}
    .navbar-default .navbar-nav>.active>a,.navbar-default .navbar-nav>.active>a:hover {color:#236A6E;font-weight:600;background:transparent;border-bottom:3px solid #236A6E;}
    .btn-primary {background:#236A6E;border-color:#236A6E;}
    .btn-primary:hover {background:#174F54;border-color:#174F54;}
    .result {background:#eef5f5;border-left:4px solid #236A6E;padding:12px 15px;margin:12px 0;}
    .note {background:#EEF8EF;border-left:4px solid #548235;padding:10px 14px;margin-top:12px;}
    .subtle {color:#58696a;font-size:13px;}
    pre {white-space:pre-wrap;}
    .shiny-input-container {max-width:100%;}
  '))),
  titlePanel("STA2ASM — Week 4: Continuous distributions"),
  withMathJax(),
  p(class = "subtle", "Explore continuous distributions, normal approximations, simulation and the Central Limit Theorem using the Week 4 readings."),
  navbarPage(title = NULL, id = "main_tab", collapsible = TRUE,
  tabPanel("Explore distributions", value="explore",
    fluidRow(column(12,h2("Explore continuous distributions"),p("Change parameters, calculate probabilities and compare density and cumulative distribution functions."))),
    fluidRow(column(3, wellPanel(
      selectInput("ex_dist","Distribution",cont_choices),
      uiOutput("ex_params"),
      selectInput("ex_mode","Explore",c("Left tail: P(X ≤ a)"="left", "Right tail: P(X > a)"="right",
                                        "Between: P(a ≤ X ≤ b)"="between", "Quantile: find x from p"="quantile")),
      conditionalPanel("input.ex_mode != 'quantile'",numericInput("ex_a",math_label("Value", "a"),1.96,step=.5)),
      conditionalPanel("input.ex_mode == 'between'",numericInput("ex_b",math_label("Value", "b"),2,step=.5)),
      conditionalPanel("input.ex_mode == 'quantile'",numericInput("ex_p",math_label("Cumulative probability", "p"),.975,min=.001,max=.999,step=.025)),
      actionButton("ex_reset","Restore reading defaults"),
      actionButton("ex_axes","Reset axes"))),
      column(9, uiOutput("ex_formula"),
        conditionalPanel("input.ex_dist == 'norm'",
          fluidRow(column(6,plotOutput("ex_pdf",height="290px")),column(6,plotOutput("ex_standard",height="290px"))),
          div(class="subtle","For a normal distribution, the shaded probability is identical after standardisation."),
          fluidRow(column(6, plotOutput("ex_cdf", height="245px")))),
        conditionalPanel("input.ex_dist != 'norm'",
          fluidRow(column(6,plotOutput("ex_pdf2",height="310px")),column(6,plotOutput("ex_cdf2",height="310px")))),
        uiOutput("ex_result")))),
  tabPanel("Normal approximations",value="approx",
    fluidRow(column(12,h2("Normal approximations to discrete distributions"))),
    fluidRow(column(3,wellPanel(
      selectInput("ap_dist","Exact distribution",c("Binomial"="bin","Poisson"="pois")),
      conditionalPanel("input.ap_dist == 'bin'",numericInput("ap_n",math_label("Number of trials", "n"),20,min=1,max=300),
                       sliderInput("ap_p",math_label("Probability of success", "p"),.01,.99,.3,step=.01)),
      conditionalPanel("input.ap_dist == 'pois'",numericInput("ap_lambda",math_label("Expected count", "\\lambda"),5,min=.1,max=100,step=.5)),
      selectInput("ap_mode","Probability",c("P(X ≤ k)"="left","P(X > k)"="right","P(a ≤ X ≤ b)"="between")),
      numericInput("ap_k",withMathJax(HTML("Value \\(k\\) (or \\(a\\))")),5,step=1),
      conditionalPanel("input.ap_mode == 'between'",numericInput("ap_b",math_label("Value", "b"),10,step=1)),
      checkboxInput("ap_cc","Use continuity correction",TRUE),
      actionButton("ap_reset","Restore reading defaults"))),
      column(9,uiOutput("ap_formula"),plotOutput("ap_plot",height="400px"),uiOutput("ap_result")))),
  tabPanel(title = HTML("\\(t\\) versus Normal"),value="tcompare",
    fluidRow(column(12,withMathJax(HTML("<h2>How the \\(t\\)-distribution approaches the standard normal</h2>")))),
    fluidRow(column(3,wellPanel(
      sliderInput("tc_df_slider",math_label("Degrees of freedom", "\\nu"),min=1,max=200,value=5,step=1),
      numericInput("tc_df",math_label("Or enter degrees of freedom", "\\nu"),5,min=1,max=10000,step=1),
      selectInput("tc_preset","Cumulative probability",c("0.025"=".025","0.975"=".975","Other"="other"),selected=".975"),
      conditionalPanel("input.tc_preset == 'other'",numericInput("tc_p",math_label("Probability", "p"),.95,min=.001,max=.999,step=.025)),
      actionButton("tc_reset","Restore reading defaults"))),
      column(9,plotOutput("tc_plot",height="420px"),uiOutput("tc_result")))),
  tabPanel("Simulation and CLT",value="clt",
    fluidRow(column(12,h2("Simulation and the Central Limit Theorem"),
      p("Explore individual observations, the sampling distribution of the mean, and confidence-interval coverage."))),
    fluidRow(column(3,wellPanel(
      selectInput("clt_dist","Underlying distribution",all_choices,selected="chisq"),
      uiOutput("clt_params"),
      numericInput("clt_n",math_label("Observations per sample", "n"),value=10,min=1,max=500,step=1),
      numericInput("clt_m",math_label("Number of simulated samples", "m"),1000,min=50,max=5000,step=50),
      actionButton("clt_rerun","Generate new simulation"),
      actionButton("clt_reset","Restore reading defaults"))),
      column(9,uiOutput("clt_formula"),
        fluidRow(column(6,plotOutput("clt_population",height="340px")),
                 column(6,plotOutput("clt_means",height="340px"))),
        div(class="subtle","Left: one random sample of individual observations and the theoretical distribution. Right: sample means from repeated samples, with an empirically fitted normal curve."),
        uiOutput("clt_note"),
        tags$hr(),h3("95% confidence-interval coverage"),
        withMathJax(HTML("<p>For each simulated sample, calculate a usual one-sample \\(t\\)-interval and check whether it contains the population mean.</p>")),
        conditionalPanel("input.clt_dist == 'lnorm'",
          selectInput("clt_ci_type","Interval and population parameter",
            c("Original-scale 𝑡-interval: arithmetic mean"="ordinary",
              "Back-transformed log-scale interval: median"="back_median",
              "Back-transformed log-scale interval: arithmetic mean"="back_mean"))),
        selectInput("clt_show","Number of intervals to display",choices=c(20,50,100),selected=50),
        uiOutput("clt_coverage"),
        plotOutput("clt_intervals",height="400px"),
        div(class="subtle","Green intervals contain the specified population parameter; coral intervals do not. The plot shows only a subset, but coverage uses all simulated samples."),
        uiOutput("clt_ci_note"))))
  )
)

server <- function(input, output, session) {
  spec_from <- function(prefix, dist) {
    val <- function(n, default) { x <- input[[paste0(prefix,"_",n)]]; if (is.null(x) || !is.finite(x)) default else x }
    make_spec(dist, mu=val("mu",0),sd=max(.01,val("sd",1)),
      df=as.integer(max(1,val("df",5))),trials=as.integer(max(1,val("trials",20))),
      prob=min(.999,max(.001,val("prob",.3))),rate=max(.001,val("rate",2.6)),
      size=as.integer(max(1,val("size",3))),nbmu=max(.001,val("nbmu",7)),theta=max(.01,val("theta",3)))
  }
  output$ex_params <- renderUI(parameter_ui("ex",input$ex_dist))
  output$clt_params <- renderUI(parameter_ui("clt",input$clt_dist))
  ex <- reactive(spec_from("ex",input$ex_dist))
  clt_spec <- reactive(spec_from("clt",input$clt_dist))
  output$ex_formula <- renderUI(math(ex()$latex))
  output$clt_formula <- renderUI(math(clt_spec()$latex))

  observeEvent(input$ex_reset, {
    updateSelectInput(session,"ex_dist",selected="norm")
    updateSelectInput(session,"ex_mode",selected="left")
    updateNumericInput(session,"ex_a",value=1.96)
    updateNumericInput(session,"ex_b",value=2)
    updateNumericInput(session,"ex_p",value=.975)
  })
  # Reference axes are established on entering a distribution, and reset on
  # demand. Their baseline does not change when its parameters change.
  ex_baseline <- reactiveVal(NULL)
  ex_axis_state <- reactiveVal(NULL)
  observeEvent(input$ex_dist, {
    s <- ex()
    lim <- view_range(s)
    xx <- seq(lim[1],lim[2],length.out=1000)
    if (identical(s$label,"Chi-squared")) xx <- pmax(xx,1e-3)
    attr(lim,"peak") <- density_plot_limit(s,xx,s$dens(xx))
    ex_baseline(lim)
    ex_axis_state(NULL)
  }, ignoreInit=FALSE)
  observeEvent(input$ex_axes, {
    s <- ex(); lim <- view_range(s)
    xx <- seq(lim[1],lim[2],length.out=1000)
    if (identical(s$label,"Chi-squared")) xx <- pmax(xx,1e-3)
    attr(lim,"peak") <- density_plot_limit(s,xx,s$dens(xx))
    ex_baseline(lim)
    ex_axis_state(NULL)
  })
  observeEvent(input$ex_reset, { ex_axis_state(NULL) })
  ex_limits <- reactive({
    s <- ex(); baseline <- ex_baseline()
    if (is.null(baseline)) baseline <- view_range(s)
    marks <- if (input$ex_mode=="between") c(input$ex_a,input$ex_b) else if (input$ex_mode=="quantile") ex()$quant(input$ex_p) else input$ex_a
    needed <- axis_limits(s,baseline,marks)
    old <- isolate(ex_axis_state())
    if (!is.null(old)) {
      needed$x <- c(min(needed$x[1],old$x[1]),max(needed$x[2],old$x[2]))
      needed$y[2] <- max(needed$y[2],old$y[2])
    }
    if (is.null(old) || !identical(old,needed)) ex_axis_state(needed)
    needed
  })
  ex_mode <- reactive(input$ex_mode)
  ex_a <- reactive(if (ex_mode()=="quantile") ex()$quant(input$ex_p) else input$ex_a)
  ex_b <- reactive(input$ex_b)
  draw_ex <- function(cdf=FALSE, standard=FALSE) {
    s <- if (standard) make_spec("norm") else ex()
    a <- ex_a(); b <- ex_b()
    if (standard) {
      mu <- ex()$mean; sd <- sqrt(ex()$var)
      a <- (a-mu)/sd; b <- (b-mu)/sd
    }
    limits <- ex_limits()
    if (standard) {
      # Use exactly corresponding scales in both normal panels.
      limits$x <- (limits$x - mu) / sd
      limits$y <- limits$y * sd
    }
    if (cdf) limits$y <- c(0, 1.02)
    plot_density(s,mode=if (ex_mode()=="quantile") "left" else ex_mode(),
                 a=a,b=b,p=input$ex_p,cdf=cdf,limits=limits,
                 main=if (standard) "Standard normal (Z)" else if (cdf) "Cumulative distribution function" else "Density and shaded probability")
  }
  output$ex_pdf <- renderPlot(draw_ex())
  output$ex_standard <- renderPlot(draw_ex(standard=TRUE))
  output$ex_cdf <- renderPlot(draw_ex(cdf=TRUE))
  output$ex_pdf2 <- renderPlot(draw_ex())
  output$ex_cdf2 <- renderPlot(draw_ex(cdf=TRUE))
  output$ex_result <- renderUI({
    s <- ex(); a <- ex_a(); b <- ex_b(); mode <- ex_mode()
    if (mode == "quantile") {
      val <- a
      expr <- sprintf("q%s(%s, %s)", c(norm="norm",t="t",chisq="chisq",lnorm="lnorm")[[input$ex_dist]],fmt(input$ex_p,3),s$code)
      desc <- sprintf("%.1fth percentile: %s",100*input$ex_p,fmt(val))
    } else {
      val <- switch(mode,left=s$cdf(a),right=1-s$cdf(a),between=abs(s$cdf(b)-s$cdf(a)))
      base <- c(norm="norm",t="t",chisq="chisq",lnorm="lnorm")[[input$ex_dist]]
      expr <- switch(mode,left=sprintf("p%s(%s, %s)",base,fmt(a,3),s$code),
                     right=sprintf("1 - p%s(%s, %s)",base,fmt(a,3),s$code),
                     between=sprintf("p%s(%s, %s) - p%s(%s, %s)",base,fmt(max(a,b),3),s$code,base,fmt(min(a,b),3),s$code))
      desc <- sprintf("Probability: %s",fmt(val))
    }
    div(class="result",strong(desc),tags$br(),tags$code(expr),
        if (input$ex_dist=="norm") tagList(tags$br(),
          math(sprintf("Z=(X-%s)/%s",fmt(s$mean,2),fmt(sqrt(s$var),2)))) )
  })

  observeEvent(input$ap_reset, {
    updateSelectInput(session,"ap_dist",selected="bin")
    updateNumericInput(session,"ap_n",value=20)
    updateSliderInput(session,"ap_p",value=.3)
    updateNumericInput(session,"ap_lambda",value=5)
    updateSelectInput(session,"ap_mode",selected="left")
    updateNumericInput(session,"ap_k",value=5)
    updateNumericInput(session,"ap_b",value=10)
    updateCheckboxInput(session,"ap_cc",value=TRUE)
  })
  ap <- reactive(approx_probability(input$ap_dist,input$ap_n,input$ap_p,input$ap_lambda,
                                     input$ap_mode,input$ap_k,input$ap_b,input$ap_cc))
  output$ap_formula <- renderUI({
    a <- ap(); math(sprintf("Y\\sim N(%.3f,%.3f^2)",a$mean,a$sd))
  })
  output$ap_plot <- renderPlot({
    a <- ap(); s <- if (input$ap_dist=="bin") make_spec("bin",trials=input$ap_n,prob=input$ap_p) else make_spec("pois",rate=input$ap_lambda)
    xlim <- view_range(s,.0005,.9995); xx <- seq(max(0,floor(xlim[1])),ceiling(xlim[2]))
    if (length(xx)>300) xx <- seq(0,300)
    pp <- s$dens(xx)
    xxcurve <- seq(min(xx)-.5,max(xx)+.5,length.out=900)
    ny <- dnorm(xxcurve,a$mean,a$sd)
    ymax <- max(c(pp,ny))*1.18
    plot(range(xxcurve),c(0,ymax),type="n",xlab="Count",ylab="Probability / density",
         main="Exact probabilities and normal approximation",las=1)
    points(xx, pp, pch = 21, bg = s$fill, col = "black", cex = 2.0, lwd = 1.1)
    good <- xxcurve>=a$bounds[1] & xxcurve<=a$bounds[2]
    if (any(good)) polygon(c(xxcurve[good][1],xxcurve[good],tail(xxcurve[good],1)),
                           c(0,ny[good],0),col=adjustcolor(PAL$lightgreen,alpha.f=.6),border=NA)
    lines(xxcurve,ny,col=PAL$blue,lwd=2.8)
    legend("topright",legend=c("Exact discrete probabilities","Normal approximation"),
           col=c("black",PAL$blue),pch=c(21,NA),pt.bg=c(s$fill,NA),
           lty=c(NA,1),lwd=c(NA,3),pt.cex=c(2,NA),bty="n",cex=.9)
  })
  output$ap_result <- renderUI({
    a <- ap(); div(class="result",strong(sprintf("Exact probability: %s",fmt(a$exact))),tags$br(),
      strong(sprintf("Normal approximation: %s",fmt(a$approx))),tags$br(),
      sprintf("Absolute difference: %s",fmt(abs(a$exact-a$approx))),tags$br(),
      if (input$ap_dist=="bin") sprintf("np = %s; n(1-p) = %s",fmt(input$ap_n*input$ap_p,2),fmt(input$ap_n*(1-input$ap_p),2))
      else sprintf("Expected count = %s",fmt(input$ap_lambda,2)),tags$br(),
      if (input$ap_cc) "Normal bounds include the 0.5 continuity correction." else "Continuity correction switched off.")
  })

  # Synchronise t controls; slider gives quick access, numeric entry allows any df.
  observeEvent(input$tc_df_slider, { if (!is.null(input$tc_df) && input$tc_df!=input$tc_df_slider) updateNumericInput(session,"tc_df",value=input$tc_df_slider) })
  observeEvent(input$tc_df, { if (!is.null(input$tc_df) && input$tc_df<=200 && input$tc_df>=1 && input$tc_df!=input$tc_df_slider) updateSliderInput(session,"tc_df_slider",value=input$tc_df) })
  observeEvent(input$tc_reset, {
    updateNumericInput(session,"tc_df",value=5)
    updateSliderInput(session,"tc_df_slider",value=5)
    updateSelectInput(session,"tc_preset",selected=".975")
    updateNumericInput(session,"tc_p",value=.95)
  })
  tc_df <- reactive(as.integer(max(1,input$tc_df)))
  tc_p <- reactive(if (input$tc_preset=="other") input$tc_p else as.numeric(input$tc_preset))
  output$tc_plot <- renderPlot({
    df <- tc_df(); p <- tc_p(); tq <- qt(p,df); zq <- qnorm(p)
    lim <- 5; x <- seq(-lim,lim,length.out=1000)
    yt <- dt(x,df); yz <- dnorm(x)
    plot(x,yz,type="l",col=PAL$blue,lwd=3,xlab="Value",ylab="Density",las=1,
         main="t and standard normal density curves",ylim=c(0,max(yz,yt)*1.12))
    lines(x,yt,col=PAL$purple,lwd=2.7)
    abline(v=tq,col=PAL$purple,lty=2,lwd=2)
    abline(v=zq,col=PAL$blue,lty=3,lwd=2)
    if (abs(tq) > lim) mtext("t quantile lies outside the displayed range",side=3,line=.1,cex=.75,col=PAL$teal)
    legend("topright",legend=c("Standard normal","t distribution","t quantile","Normal quantile"),
           col=c(PAL$blue,PAL$purple,PAL$purple,PAL$blue),lty=c(1,1,2,3),lwd=2,bty="n",cex=.85)
  })
  output$tc_result <- renderUI({
    df <- tc_df(); p <- tc_p()
    div(class="result",math(sprintf("t_{%d,%.3f}=%s,\\qquad z_{%.3f}=%s",df,p,fmt(qt(p,df)),p,fmt(qnorm(p)))),
        sprintf("Absolute difference: %s",fmt(abs(qt(p,df)-qnorm(p)))),tags$br(),
        tags$code(sprintf("qt(%s, df = %d)",fmt(p,3),df)),"  |  ",tags$code(sprintf("qnorm(%s)",fmt(p,3))))
  })

  observeEvent(input$clt_reset, {
    updateSelectInput(session,"clt_dist",selected="chisq")
    updateNumericInput(session,"clt_n",value=10)
    updateNumericInput(session,"clt_m",value=1000)
  })
  # All three visualisations use the same collection of m independent samples.
  # The first column supplies the individual-observation plot.
  clt_data <- reactive({
    input$clt_rerun
    s <- clt_spec()
    n <- as.integer(input$clt_n)
    m <- as.integer(input$clt_m)
    req(is.finite(n), n >= 1, is.finite(m), m >= 1)
    observations <- matrix(s$gen(n * m), nrow = n, ncol = m)
    means <- colMeans(observations)
    sds <- if (n >= 2) apply(observations, 2, sd) else rep(NA_real_, m)
    # For n >= 2, intervals with zero sample variance are point intervals.
    margin <- if (n >= 2) qt(.975, df = n - 1) * sds / sqrt(n) else rep(NA_real_, m)
    ordinary <- cbind(lower = means - margin, upper = means + margin)
    result <- list(spec = s, n = n, m = m, first = observations[, 1],
                   means = means, ordinary = ordinary)
    if (identical(input$clt_dist, "lnorm") && n >= 2) {
      logs <- log(observations)
      log_means <- colMeans(logs)
      log_sds <- apply(logs, 2, sd)
      log_margin <- qt(.975, df = n - 1) * log_sds / sqrt(n)
      result$back <- cbind(lower = exp(log_means - log_margin),
                           upper = exp(log_means + log_margin))
    }
    result
  })

  output$clt_population <- renderPlot({
    d <- clt_data(); s <- d$spec; x <- d$first
    if (s$discrete) {
      # Relative-frequency bars and theoretical probability masses share a scale.
      observed <- table(x)
      support <- sort(unique(c(as.numeric(names(observed)),
                               seq(floor(s$quant(.005)), ceiling(s$quant(.995))))))
      support <- support[is.finite(support)]
      if (length(support) > 220) {
        support <- seq(floor(s$quant(.005)), ceiling(s$quant(.995)), length.out = 220)
        support <- unique(round(support))
      }
      actual <- as.numeric(table(factor(x, levels = support))) / length(x)
      theory <- s$dens(support)
      ylim <- c(0, max(c(actual, theory), na.rm = TRUE) * 1.18)
      plot(range(support) + c(-.7, .7), ylim, type = "n", las = 1,
           main = "One random sample", xlab = "Individual observation", ylab = "Relative frequency / probability")
      points(support, theory, pch = 21, bg = s$fill, col = "black", cex = 2.0, lwd = 1.1)
      points(support, actual, pch = 23, bg = PAL$teal, col = "black", cex = 1.35, lwd = 1.1)
      legend("topright", c("Observed proportion", "Theoretical probability"),
             pch = c(23, 21), pt.bg = c(PAL$teal, s$fill),
             pt.cex = c(1.35, 2), bty = "n", cex = .78)
    } else {
      h <- hist(x, plot = FALSE)
      lim <- range(c(h$breaks, view_range(s)))
      # Avoid extreme right tails overwhelming the displayed density.
      if (identical(s$label, "Log-normal")) lim <- c(0, max(h$breaks, view_range(s)[2]))
      xx <- seq(max(if (identical(s$label, "t") || identical(s$label, "Normal")) lim[1] else 1e-3,
                        lim[1]), lim[2], length.out = 850)
      yy <- s$dens(xx); yy[!is.finite(yy)] <- NA_real_
      ylim <- c(0, max(max(h$density, na.rm = TRUE),
                        density_plot_limit(s, xx, yy)) * 1.08)
      hist(x, breaks = h$breaks, freq = FALSE,
           col = adjustcolor(s$fill, alpha.f = .75), border = "black",
           main = "One random sample", xlab = "Individual observation",
           ylab = "Density", las = 1, xlim = lim, ylim = ylim)
      lines(xx, yy, col = s$col, lwd = 3)
      legend("topright", c("Sample", "Theoretical PDF"),
             fill = c(s$fill, NA), border = c("black", NA),
             lty = c(NA, 1), col = c(NA, s$col), lwd = c(NA, 3),
             bty = "n", cex = .78)
    }
  })

  output$clt_means <- renderPlot({
    x <- clt_data()$means
    h <- hist(x, plot = FALSE)
    fitted_sd <- sd(x)
    xx <- seq(min(h$breaks), max(h$breaks), length.out = 600)
    normal_y <- if (is.finite(fitted_sd) && fitted_sd > 0)
      dnorm(xx, mean(x), fitted_sd) else rep(0, length(xx))
    hist(x, breaks = h$breaks, freq = FALSE,
         col = adjustcolor(PAL$lightgreen, alpha.f = .75), border = "black",
         main = "Sampling distribution of the mean", xlab = "Sample mean",
         ylab = "Density", las = 1,
         ylim = c(0, max(c(h$density, normal_y), na.rm = TRUE) * 1.12))
    if (is.finite(fitted_sd) && fitted_sd > 0) {
      lines(xx, normal_y, col = PAL$blue, lwd = 3)
      legend("topright", "Fitted normal curve", col = PAL$blue,
             lty = 1, lwd = 3, bty = "n", cex = .8)
    }
  })

  output$clt_note <- renderUI({
    s <- clt_data()$spec
    if (identical(s$label, "t") && !is.finite(s$mean))
      div(class = "note", strong("t with 1 degree of freedom: "),
          "The population mean and variance do not exist. The sample mean remains Cauchy-distributed regardless of n; the usual CLT does not apply.")
    else if (identical(s$label, "t") && !is.finite(s$var))
      div(class = "note", strong("t with 2 degrees of freedom: "),
          "The population mean exists but the variance is infinite. The usual finite-variance CLT and standard-error formula do not apply.")
    else NULL
  })

  ci_selected <- reactive({
    d <- clt_data()
    kind <- if (identical(input$clt_dist, "lnorm")) input$clt_ci_type else "ordinary"
    if (is.null(kind)) kind <- "ordinary"
    target <- if (identical(kind, "back_median")) d$spec$quant(.5) else d$spec$mean
    intervals <- if (identical(kind, "ordinary")) d$ordinary else d$back
    list(kind = kind, target = target, intervals = intervals, d = d)
  })

  output$clt_coverage <- renderUI({
    c <- ci_selected()
    if (c$d$n < 2)
      return(div(class = "note", "At least two observations per sample are needed to calculate the usual t-interval."))
    if (!is.finite(c$target))
      return(div(class = "note", "The population mean is undefined for t with 1 degree of freedom; mean coverage is not defined."))
    covered <- c$intervals[, 1] <= c$target & c$target <= c$intervals[, 2]
    coverage <- mean(covered)
    label <- switch(c$kind, ordinary = "Original-scale \\(t\\)-interval for the mean",
                    back_median = "Back-transformed interval for the median",
                    back_mean = "Back-transformed interval checked against the mean")
    div(class = "result", withMathJax(HTML(paste0("<strong>", label, "</strong>"))), tags$br(),
        sprintf("Estimated coverage: %.1f%% (%d of %d intervals); nominal coverage: 95%%.",
                100 * coverage, sum(covered), length(covered)),
        if (identical(input$clt_dist, "lnorm")) tagList(tags$br(),
          withMathJax(HTML(sprintf("Original-scale \\(t\\)-interval for mean: %.1f%%; back-transformed interval for median: %.1f%%; back-transformed interval checked against mean: %.1f%%.",
            100 * mean(c$d$ordinary[, 1] <= c$d$spec$mean & c$d$spec$mean <= c$d$ordinary[, 2]),
            100 * mean(c$d$back[, 1] <= c$d$spec$quant(.5) & c$d$spec$quant(.5) <= c$d$back[, 2]),
            100 * mean(c$d$back[, 1] <= c$d$spec$mean & c$d$spec$mean <= c$d$back[, 2]))))))
  })

  output$clt_intervals <- renderPlot({
    c <- ci_selected()
    if (c$d$n < 2 || !is.finite(c$target)) {
      plot.new(); text(.5, .5, "Confidence intervals are not available for this setting.")
      return(invisible(NULL))
    }
    k <- min(as.integer(input$clt_show), nrow(c$intervals))
    bounds <- c$intervals[seq_len(k), , drop = FALSE]
    covered <- bounds[, 1] <= c$target & c$target <= bounds[, 2]
    finite_bounds <- as.vector(bounds[is.finite(bounds)])
    lim <- range(c(finite_bounds, c$target))
    if (diff(lim) == 0) lim <- lim + c(-1, 1)
    pad <- .06 * diff(lim)
    plot(c(lim[1] - pad, lim[2] + pad), c(.5, k + .5), type = "n",
         xlab = if (identical(c$kind, "back_median")) "Population median / confidence interval" else "Population mean / confidence interval", ylab = "Simulated sample",
         main = "95% confidence intervals from repeated samples", las = 1)
    abline(v = c$target, col = PAL$blue, lty = 2, lwd = 2)
    parameter_name <- if (identical(c$kind, "back_median")) "median" else "mean"
    # A zero-variance sample produces a zero-width CI; make the point visible.
    cols <- ifelse(covered, PAL$green, PAL$coral)
    segments(bounds[, 1], seq_len(k), bounds[, 2], seq_len(k), col = cols, lwd = 2)
    points(rowMeans(bounds), seq_len(k), col = cols, pch = 16,
           cex = ifelse(bounds[, 1] == bounds[, 2], 1.25, .5))
    legend("topright", c(paste("Contains", parameter_name),
                         paste("Misses", parameter_name),
                         paste("Population", parameter_name)),
           col = c(PAL$green, PAL$coral, PAL$blue), lty = c(1, 1, 2),
           lwd = 2, bty = "n", cex = .8)
  })

  output$clt_ci_note <- renderUI({
    d <- clt_data()
    if (d$n == 1)
      div(class = "note", "With n = 1, a usual one-sample t confidence interval cannot be calculated.")
    else if (d$spec$discrete && any(d$ordinary[, 1] == d$ordinary[, 2], na.rm = TRUE))
      div(class = "note", withMathJax(HTML("Some simulated samples have zero sample variance, so their usual \\(t\\)-interval collapses to a single point. This is common for small samples from discrete distributions, particularly with extreme probabilities. These point intervals are shown as larger dots and included in the coverage calculation.")))
    else if (identical(input$clt_dist, "lnorm"))
      div(class = "note", withMathJax(HTML("Exponentiating a \\(t\\)-interval for the log-scale mean targets the original-scale median, not the arithmetic mean. The same back-transformed interval may sometimes contain the arithmetic mean, but its coverage for that parameter is not generally 95%.")))
    else if (identical(input$clt_dist, "t") && !is.finite(d$spec$var))
      div(class = "note", "The usual finite-variance CLT does not apply here. Interpret the simulated coverage cautiously.")
    else NULL
  })
}

shinyApp(ui,server)
