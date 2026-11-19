# =============================================================================
# Time of concentration (original vs. adjusted Kirpich) as a function of velocity
# -----------------------------------------------------------------------------
# Reads the data from a CSV file, fits linear regressions TC ~ Velocity for
# both methods, builds a scatter plot labelled with the reach (TC) names, and
# produces residual diagnostic plots.
#
# How to use:
#   1. Open the project at the repository root folder (or setwd() to it).
#   2. Adjust the PARAMETERS below, if needed.
#   3. Run: source("R/analise_tc_kirpich.R")
#
# Outputs:
#   figuras/    -> main plot and diagnostic plots (PNG)
#   resultados/ -> table with the regression coefficients (CSV)
#
# Author: Kamila Almeida dos Santos
# =============================================================================

# ---- Packages ----------------------------------------------------------------
pacotes <- c("ggplot2", "dplyr", "tidyr", "ggrepel")
faltando <- pacotes[!vapply(pacotes, requireNamespace, logical(1), quietly = TRUE)]
if (length(faltando) > 0) install.packages(faltando)  # install missing packages

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)

# ---- Parameters --------------------------------------------------------------
arquivo_dados  <- "dados/tc_bot.csv"   # input CSV file
nome_area      <- "BOT"                # used in the plot title and file names
zeros_como_na  <- TRUE                 # treat TC = 0 as missing data?
excluir_v_zero <- FALSE                # drop reaches with velocity = 0?

pasta_figuras    <- "figuras"
pasta_resultados <- "resultados"
dir.create(pasta_figuras,    showWarnings = FALSE)
dir.create(pasta_resultados, showWarnings = FALSE)

# ---- Data loading and preparation --------------------------------------------
# Expected columns: Nome_TC, Velocidade_ms,
#                   TC_min_kirpich_original, TC_min_kirpich_ajustado
dados <- read.csv(arquivo_dados, stringsAsFactors = FALSE)

colunas_tc <- c("TC_min_kirpich_original", "TC_min_kirpich_ajustado")

# Replace zeros with NA so they are excluded from the regressions
if (zeros_como_na) {
  dados <- dados %>%
    mutate(across(all_of(colunas_tc), ~ na_if(.x, 0)))
}
# Optionally remove reaches with zero velocity
if (excluir_v_zero) {
  dados <- dados %>% filter(Velocidade_ms > 0)
}

# Display labels and colours for each method
rotulos_metodo <- c(
  TC_min_kirpich_original = "TC Kirpich Original",
  TC_min_kirpich_ajustado = "TC Kirpich Ajustado"
)
cores_metodo <- c("TC Kirpich Original" = "blue",
                  "TC Kirpich Ajustado" = "red")

# Reshape to long format (one row per reach and method)
dados_longos <- dados %>%
  pivot_longer(cols = all_of(colunas_tc),
               names_to = "Tipo_de_TC", values_to = "Valor_TC") %>%
  mutate(Metodo = factor(rotulos_metodo[Tipo_de_TC], levels = rotulos_metodo)) %>%
  filter(!is.na(Valor_TC))

# ---- Linear regressions ------------------------------------------------------
# Fit TC ~ Velocity for a given TC column
ajustar_modelo <- function(coluna) {
  lm(reformulate("Velocidade_ms", response = coluna), data = dados)
}
modelos <- lapply(setNames(colunas_tc, colunas_tc), ajustar_modelo)

# Summarise a model into a one-row data frame
resumo_modelo <- function(modelo, nome) {
  s <- summary(modelo)
  data.frame(
    Metodo     = rotulos_metodo[[nome]],
    n          = nobs(modelo),
    Intercepto = coef(modelo)[[1]],
    Inclinacao = coef(modelo)[[2]],
    R2         = s$r.squared,
    p_valor    = s$coefficients[2, 4]
  )
}
tabela_resultados <- do.call(rbind, Map(resumo_modelo, modelos, names(modelos)))
rownames(tabela_resultados) <- NULL
print(tabela_resultados)

# Save the results table
write.csv(tabela_resultados,
          file.path(pasta_resultados, paste0("regressoes_", nome_area, ".csv")),
          row.names = FALSE)

# Equation text with the correct sign for the slope
texto_equacao <- function(b0, b1) {
  sprintf("y = %.1f %s %.1fx", b0, ifelse(b1 < 0, "-", "+"), abs(b1))
}

# ---- Main plot ---------------------------------------------------------------
# Axis limits with headroom for the annotations
y_max <- ceiling(max(dados_longos$Valor_TC, na.rm = TRUE) / 10) * 10 + 20
x_max <- ceiling(max(dados$Velocidade_ms, na.rm = TRUE) * 2) / 2 + 0.5

# Annotation text (method, R² and equation) and its x position
anotacoes <- tabela_resultados %>%
  mutate(
    x     = c(0, x_max * 0.55),
    rotulo = sprintf("%s\nR² = %.3f\n%s", Metodo, R2,
                     texto_equacao(Intercepto, Inclinacao))
  )

grafico <- ggplot(dados_longos,
                  aes(x = Velocidade_ms, y = Valor_TC, color = Metodo)) +
  # Regression lines with 95% confidence bands
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 1) +
  geom_point(size = 2, alpha = 0.7) +
  # Reach names, repelled to avoid overlapping labels
  geom_text_repel(aes(label = Nome_TC), size = 2.3, show.legend = FALSE,
                  max.overlaps = Inf, min.segment.length = 0,
                  segment.alpha = 0.4, seed = 42) +
  # Regression statistics in the top of the panel
  geom_text(data = anotacoes,
            aes(x = x, y = y_max, label = rotulo, color = Metodo),
            hjust = 0, vjust = 1, size = 3.3, lineheight = 1.1,
            fontface = "bold", show.legend = FALSE, inherit.aes = FALSE) +
  scale_color_manual(values = cores_metodo) +
  scale_x_continuous(breaks = seq(0, x_max, 0.5), limits = c(-0.1, x_max)) +
  scale_y_continuous(breaks = seq(0, y_max, 10), limits = c(0, y_max)) +
  labs(
    title = paste("TC Kirpich Original vs. TC Kirpich Ajustado -", nome_area),
    x     = "Velocidade (m/s)",
    y     = "Tempo de Concentração (min)",
    color = "Método de Cálculo"
  ) +
  theme_minimal() +
  theme(
    plot.title       = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title       = element_text(size = 12),
    axis.text        = element_text(size = 10),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "gray50", fill = NA, linewidth = 0.5),
    legend.position  = "bottom"
  )

print(grafico)

# Save the main plot
ggsave(file.path(pasta_figuras, paste0("tc_vs_velocidade_", nome_area, ".png")),
       grafico, width = 10, height = 7, dpi = 300, bg = "white")

# ---- Residual diagnostics ----------------------------------------------------
# Draw the four standard lm() diagnostic plots, on screen and to a PNG file
salvar_diagnostico <- function(modelo, nome) {
  titulo  <- paste("Diagnóstico -", rotulos_metodo[[nome]], "-", nome_area)
  arquivo <- file.path(pasta_figuras,
                       paste0("diagnostico_", sub("TC_min_kirpich_", "", nome),
                              "_", nome_area, ".png"))
  desenhar <- function() {
    op <- par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))
    on.exit(par(op))  # restore the plot layout afterwards
    plot(modelo)
    mtext(titulo, outer = TRUE, font = 2)
  }
  desenhar()                                   # show on screen
  png(arquivo, width = 2400, height = 2400, res = 300)
  desenhar()                                   # write to file
  dev.off()
}
invisible(Map(salvar_diagnostico, modelos, names(modelos)))

message("Done! Plots saved in '", pasta_figuras,
        "/' and results in '", pasta_resultados, "/'.")
