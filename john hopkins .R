# ==============================================================================
# 0. DEPENDENCY MANAGEMENT
# ==============================================================================
required_pkgs <- c("tidyverse", "lubridate", "zoo", "tseries", "forecast", 
                   "ggplot2", "plotly", "shiny", "bslib", "patchwork")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs) > 0) install.packages(new_pkgs, quiet = TRUE)

library(tidyverse)
library(lubridate)
library(zoo)
library(tseries)
library(forecast)
library(ggplot2)
library(plotly)
library(shiny)
library(bslib)
library(patchwork)

if (!is.null(dev.list())) graphics.off()

# ==============================================================================
# 1. SAFE DATA INGESTION & DATA CLEANING FUNCTION
# ==============================================================================
load_jhu_data <- function() {
  cat("\n--- Loading Johns Hopkins Datasets ---\n")
  
  jhu_base_url <- "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/"
  url_cases     <- paste0(jhu_base_url, "time_series_covid19_confirmed_global.csv")
  url_deaths    <- paste0(jhu_base_url, "time_series_covid19_deaths_global.csv")
  url_recovered <- paste0(jhu_base_url, "time_series_covid19_recovered_global.csv")
  
  tidy_jhu <- function(url, metric_name) {
    df <- tryCatch({
      read_csv(url, show_col_types = FALSE)
    }, error = function(e) {
      stop(paste("Failed to download data from:", url))
    })
    
    # Standardize column headers
    colnames(df)[1:4] <- c("province", "country", "lat", "long")
    
    df %>%
      pivot_longer(
        cols = -(1:4),
        names_to = "date",
        values_to = metric_name
      ) %>%
      mutate(
        date = mdy(date),
        !!sym(metric_name) := as.numeric(!!sym(metric_name))
      ) %>%
      filter(!is.na(date)) %>%
      group_by(country, date) %>%
      summarize(!!sym(metric_name) := sum(!!sym(metric_name), na.rm = TRUE), .groups = "drop")
  }
  
  df_cases     <- tidy_jhu(url_cases, "total_cases")
  df_deaths    <- tidy_jhu(url_deaths, "total_deaths")
  df_recovered <- tidy_jhu(url_recovered, "total_recovered")
  
  df_clean <- df_cases %>%
    full_join(df_deaths, by = c("country", "date")) %>%
    full_join(df_recovered, by = c("country", "date")) %>%
    arrange(country, date) %>%
    group_by(country) %>%
    mutate(
      total_cases     = na.locf(total_cases, na.rm = FALSE),
      total_deaths    = na.locf(total_deaths, na.rm = FALSE),
      total_cases     = ifelse(is.na(total_cases), 0, total_cases),
      total_deaths    = ifelse(is.na(total_deaths), 0, total_deaths),
      
      total_recovered = ifelse(date <= "2021-08-05", na.locf(total_recovered, na.rm = FALSE), NA_real_),
      total_recovered = ifelse(date <= "2021-08-05" & is.na(total_recovered), 0, total_recovered),
      
      total_active    = ifelse(!is.na(total_recovered), total_cases - total_deaths - total_recovered, NA_real_),
      total_active    = ifelse(total_active < 0, 0, total_active),
      
      new_cases     = total_cases - lag(total_cases, default = 0),
      new_deaths    = total_deaths - lag(total_deaths, default = 0),
      new_recovered = ifelse(!is.na(total_recovered), total_recovered - lag(total_recovered, default = 0), NA_real_),
      
      new_cases     = ifelse(new_cases < 0, 0, new_cases),
      new_deaths    = ifelse(new_deaths < 0, 0, new_deaths),
      new_recovered = ifelse(!is.na(new_recovered) & new_recovered < 0, 0, new_recovered)
    ) %>%
    ungroup() %>%
    filter(date >= "2020-03-01" & date <= "2023-03-09")
  
  # Outlier cleaning
  clean_outliers <- function(series) {
    if (all(is.na(series)) || length(series) < 14) return(series)
    ts_obj <- ts(ifelse(is.na(series), 0, series), frequency = 7)
    cleaned <- suppressWarnings(as.numeric(tsclean(ts_obj)))
    if (any(is.na(series))) cleaned[is.na(series)] <- NA_real_
    return(cleaned)
  }
  
  df_treated <- df_clean %>%
    group_by(country) %>%
    mutate(
      new_cases_cleaned     = clean_outliers(new_cases),
      new_deaths_cleaned    = clean_outliers(new_deaths),
      new_recovered_cleaned = clean_outliers(new_recovered)
    ) %>%
    ungroup()
  
  # Feature Engineering
  df_features <- df_treated %>%
    group_by(country) %>%
    mutate(
      cases_ma7     = rollmean(new_cases_cleaned, k = 7, fill = NA, align = "right"),
      deaths_ma7    = rollmean(new_deaths_cleaned, k = 7, fill = NA, align = "right"),
      recovered_ma7 = rollmean(new_recovered_cleaned, k = 7, fill = NA, align = "right"),
      
      cases_lag7  = lag(new_cases_cleaned, 7),
      growth_rate = ifelse(!is.na(cases_lag7) & cases_lag7 > 0, ((cases_ma7 - lag(cases_ma7, 7)) / lag(cases_ma7, 7)) * 100, 0),
      
      case_fatality_rate = ifelse(total_cases > 0, (total_deaths / total_cases) * 100, 0),
      recovery_rate      = ifelse(!is.na(total_recovered) & total_cases > 0, (total_recovered / total_cases) * 100, NA_real_)
    ) %>%
    ungroup()
  
  return(df_features)
}

# Run data ingestion
df_features <- load_jhu_data()

# ==============================================================================
# 2. SHINY DASHBOARD INTERFACE
# ==============================================================================
ui <- fluidPage(
  theme = bs_theme(bootswatch = "flatly"),
  titlePanel("COVID-19 Surveillance & Forecasting Platform"),
  
  fluidRow(
    column(
      width = 4,
      div(
        class = "card text-white bg-primary mb-3 text-center",
        div(class = "card-header", strong("Total Confirmed Cases")),
        div(class = "card-body", h3(textOutput("kpi_confirmed", inline = TRUE)))
      )
    ),
    column(
      width = 4,
      div(
        class = "card text-white bg-danger mb-3 text-center",
        div(class = "card-header", strong("Total Deaths")),
        div(class = "card-body", h3(textOutput("kpi_deaths", inline = TRUE)))
      )
    ),
    column(
      width = 4,
      div(
        class = "card text-white bg-success mb-3 text-center",
        div(class = "card-header", strong("Total Recovered")),
        div(class = "card-body", h3(textOutput("kpi_recovered", inline = TRUE)))
      )
    )
  ),
  hr(),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput("country_select", "Select Country for Time-Series:", 
                  choices = sort(unique(df_features$country)), 
                  selected = if("India" %in% df_features$country) "India" else unique(df_features$country)[1]),
      dateRangeInput("date_range", "Select Date Range:",
                     start = "2020-03-01", end = "2023-03-09",
                     min = "2020-01-22", max = "2023-03-09"),
      sliderInput("fc_days", "Forecast Horizon (Days):", min = 7, max = 60, value = 30, step = 7),
      hr(),
      helpText("Data source: JHU CSSE repository.")
    ),
    
    mainPanel(
      width = 9,
      tabsetPanel(
        tabPanel("Epidemic Dynamics",
                 fluidRow(
                   column(12, plotlyOutput("threeStreamsPlot", height = "350px")),
                   column(12, plotlyOutput("activeVsTotalPlot", height = "300px"))
                 )
        ),
        tabPanel("Fatality & Recovery Rates",
                 plotlyOutput("ratesPlot", height = "400px"),
                 plotlyOutput("growthRatePlot", height = "280px")
        ),
        tabPanel("ARIMA Forecast Model",
                 plotlyOutput("arimaForecastPlot", height = "450px"),
                 h4("Model Diagnostics & Summary"),
                 verbatimTextOutput("modelMetrics")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  
  country_data <- reactive({
    req(input$country_select)
    df_features %>%
      filter(country == input$country_select) %>%
      filter(date >= input$date_range[1] & date <= input$date_range[2])
  })
  
  kpi_metrics <- reactive({
    df <- country_data()
    req(nrow(df) > 0)
    
    latest_confirmed <- max(df$total_cases, na.rm = TRUE)
    latest_deaths    <- max(df$total_deaths, na.rm = TRUE)
    
    rec_vals <- df$total_recovered[!is.na(df$total_recovered)]
    latest_recovered <- if (length(rec_vals) > 0) max(rec_vals, na.rm = TRUE) else NA_real_
    
    list(
      confirmed = latest_confirmed,
      deaths    = latest_deaths,
      recovered = latest_recovered
    )
  })
  
  output$kpi_confirmed <- renderText({
    vals <- kpi_metrics()
    format(vals$confirmed, big.mark = ",")
  })
  
  output$kpi_deaths <- renderText({
    vals <- kpi_metrics()
    format(vals$deaths, big.mark = ",")
  })
  
  output$kpi_recovered <- renderText({
    vals <- kpi_metrics()
    if (is.na(vals$recovered)) {
      "N/A (Halted Aug 2021)"
    } else {
      format(vals$recovered, big.mark = ",")
    }
  })
  
  output$threeStreamsPlot <- renderPlotly({
    df <- country_data()
    p <- ggplot(df, aes(x = date)) +
      geom_line(aes(y = cases_ma7, color = "Confirmed Cases (7d MA)"), linewidth = 1) +
      geom_line(data = df %>% filter(!is.na(recovered_ma7)), 
                aes(y = recovered_ma7, color = "Recoveries (7d MA)"), linewidth = 0.9) +
      geom_line(aes(y = deaths_ma7, color = "Deaths (7d MA)"), linewidth = 0.9) +
      scale_color_manual(values = c(
        "Confirmed Cases (7d MA)" = "#1f78b4",
        "Recoveries (7d MA)"      = "#33a02c",
        "Deaths (7d MA)"          = "#e31a1c"
      )) +
      labs(title = paste("Daily Epidemiological Streams:", input$country_select), 
           x = "Date", y = "Daily Count", color = "Stream") +
      theme_minimal()
    ggplotly(p)
  })
  
  output$activeVsTotalPlot <- renderPlotly({
    df <- country_data()
    p <- ggplot(df, aes(x = date)) +
      geom_area(aes(y = total_cases, fill = "Total Confirmed"), alpha = 0.3) +
      geom_area(data = df %>% filter(!is.na(total_active)), 
                aes(y = total_active, fill = "Active Cases (until Aug 2021)"), alpha = 0.6) +
      scale_fill_manual(values = c("Total Confirmed" = "#a6cee3", "Active Cases (until Aug 2021)" = "#fb9a99")) +
      labs(title = "Cumulative Trajectory: Total vs. Active Cases", 
           x = "Date", y = "Total Population Impacted", fill = "Category") +
      theme_minimal()
    ggplotly(p)
  })
  
  output$ratesPlot <- renderPlotly({
    df <- country_data()
    p <- ggplot(df, aes(x = date)) +
      geom_line(aes(y = case_fatality_rate, color = "Case Fatality Rate (%)"), linewidth = 0.9) +
      geom_line(data = df %>% filter(!is.na(recovery_rate)), 
                aes(y = recovery_rate, color = "Recovery Rate (%) [Active until Aug 2021]"), linewidth = 0.9) +
      scale_color_manual(values = c(
        "Case Fatality Rate (%)" = "#e31a1c", 
        "Recovery Rate (%) [Active until Aug 2021]" = "#33a02c"
      )) +
      labs(title = "Clinical Outcomes: Fatality Rate vs. Recovery Rate (%)", 
           x = "Date", y = "Percentage (%)", color = "Metric") +
      theme_minimal()
    ggplotly(p)
  })
  
  output$growthRatePlot <- renderPlotly({
    df <- country_data()
    p <- ggplot(df, aes(x = date, y = growth_rate)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      geom_line(color = "#ff7f00", linewidth = 0.8) +
      labs(title = "Week-over-Week Growth Rate (%)", x = "Date", y = "% Change") +
      theme_minimal()
    ggplotly(p)
  })
  
  model_fit <- reactive({
    df <- country_data() %>% filter(!is.na(cases_ma7))
    req(nrow(df) >= 14)
    
    clean_series <- na.omit(df$cases_ma7)
    ts_target <- ts(clean_series, frequency = 7)
    
    fit <- auto.arima(ts_target, seasonal = TRUE, stepwise = TRUE, approximation = TRUE)
    fc <- forecast(fit, h = input$fc_days)
    
    list(model = fit, forecast = fc, df = df)
  })
  
  output$arimaForecastPlot <- renderPlotly({
    res <- model_fit()
    df  <- res$df
    fc  <- res$forecast
    
    hist_dates <- df$date[!is.na(df$cases_ma7)]
    hist_cases <- round(na.omit(df$cases_ma7))
    
    last_date <- max(hist_dates)
    fc_dates  <- seq.Date(from = last_date + 1, by = "day", length.out = input$fc_days)
    
    fc_mean   <- round(as.numeric(fc$mean))
    lower_95  <- pmax(0, round(as.numeric(fc$lower[, 2])))
    upper_95  <- round(as.numeric(fc$upper[, 2]))
    lower_80  <- pmax(0, round(as.numeric(fc$lower[, 1])))
    upper_80  <- round(as.numeric(fc$upper[, 1]))
    
    plot_ly() %>%
      add_ribbons(
        x = fc_dates, ymin = lower_95, ymax = upper_95,
        name = "95% Confidence Interval", fillcolor = "rgba(179, 205, 227, 0.4)",
        line = list(color = "transparent"), hoverinfo = "text",
        text = paste0("<b>Date:</b> ", fc_dates, "<br><b>95% CI:</b> [", format(lower_95, big.mark = ","), " - ", format(upper_95, big.mark = ","), "]")
      ) %>%
      add_ribbons(
        x = fc_dates, ymin = lower_80, ymax = upper_80,
        name = "80% Confidence Interval", fillcolor = "rgba(128, 177, 211, 0.5)",
        line = list(color = "transparent"), hoverinfo = "text",
        text = paste0("<b>Date:</b> ", fc_dates, "<br><b>80% CI:</b> [", format(lower_80, big.mark = ","), " - ", format(upper_80, big.mark = ","), "]")
      ) %>%
      add_lines(
        x = hist_dates, y = hist_cases, name = "Historical Cases (7d MA)",
        line = list(color = "#1f78b4", width = 2), hoverinfo = "text",
        text = paste0("<b>Date:</b> ", hist_dates, "<br><b>Cases (7d MA):</b> ", format(hist_cases, big.mark = ","))
      ) %>%
      add_lines(
        x = fc_dates, y = fc_mean, name = "ARIMA Point Forecast",
        line = list(color = "#e31a1c", width = 2.5, dash = "dash"), hoverinfo = "text",
        text = paste0("<b>Forecast Date:</b> ", fc_dates, "<br><b>Projected Cases:</b> ", format(fc_mean, big.mark = ","))
      ) %>%
      layout(
        title = list(text = paste(input$fc_days, "Day Forward Forecast for", input$country_select)),
        xaxis = list(title = "Date", showgrid = TRUE),
        yaxis = list(title = "Daily Cases", showgrid = TRUE),
        hovermode = "x unified",
        legend = list(orientation = "h", y = -0.2)
      )
  })
  
  output$modelMetrics <- renderPrint({
    res <- model_fit()
    summary(res$model)
  })
}

shinyApp(ui = ui, server = server)