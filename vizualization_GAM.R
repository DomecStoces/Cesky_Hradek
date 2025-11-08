library(gratia)
library(dplyr)
library(ggplot2)
library(tibble)

## 1) Smooth estimates for s(Altitude_scaled) ----
sm <- smooth_estimates(
  mod_gam1,
  select = "s(Altitude_scaled)",  # <— use `select`, not `smooth`
  n = 200,
  unconditional = TRUE
)

## 2) Get inverse link and intercept ----
ilink <- family(mod_gam1)$linkinv
b0    <- coef(mod_gam1)[["(Intercept)"]]

## 3) Build both LINK- and RESPONSE-scale summaries ----
sm <- sm %>%
  mutate(
    eta   = b0 + .estimate,                 
    fit   = ilink(eta),                     
    lower = ilink(b0 + .estimate - 1.96 * .se),
    upper = ilink(b0 + .estimate + 1.96 * .se)
  )

## 4) Partial residuals for this smooth ----
# raw partials (may be a tibble with various schemas)
pres_raw <- gratia::partial_residuals(mod_gam1, select = "s(Altitude_scaled)")

# Normalize to a tibble with `.partial_residual` + `Altitude_scaled`
if (".partial_residual" %in% names(pres_raw)) {
  pres <- pres_raw
} else if (".partial" %in% names(pres_raw)) {
  pres <- dplyr::rename(pres_raw, .partial_residual = .partial)
} else if ("s(Altitude_scaled)" %in% names(pres_raw)) {
  pres <- tibble(
    Altitude_scaled  = model.frame(mod_gam1)$Altitude_scaled,
    .partial_residual = pres_raw[["s(Altitude_scaled)"]]
  )
} else {
  stop("Unexpected column names from partial_residuals(): ", paste(names(pres_raw), collapse = ", "))
}

# Add intercept & inverse-link to get response-scale partials
ilink <- family(mod_gam1)$linkinv
b0    <- coef(mod_gam1)[["(Intercept)"]]

pres <- pres %>%
  dplyr::mutate(
    eta_partial = b0 + .partial_residual,
    y_partial   = ilink(eta_partial)
  )

## 5B) Plot on response scale (adds intercept, inverse link)
p_resp <- ggplot() +
  geom_ribbon(
    data = sm,
    aes(x = Altitude_scaled, ymin = lower, ymax = upper),
    fill = "grey70", alpha = 0.35
  ) +
  geom_line(
    data = sm,
    aes(x = Altitude_scaled, y = fit),
    linewidth = 1.1
  ) +
  geom_jitter(
    data = pres,
    aes(x = Altitude_scaled, y = y_partial),
    width = 0.03, height = 0, size = 1.8, alpha = 0.6
  ) +
  labs(
    x = "Altitude (scaled)",
    y = "Wings CWM"
  ) +
  scale_x_continuous(breaks = seq(-2, 2, 1), minor_breaks = NULL) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(5),
    expand = expansion(mult = c(0, 0.02))
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  theme(
    panel.background = element_blank(),
    plot.background  = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line        = element_line(colour = "black", linewidth = 0.6),
    axis.ticks       = element_line(colour = "black", linewidth = 0.5),
    axis.ticks.length= unit(4, "pt"),
    axis.title       = element_text(size = 15),
    axis.text        = element_text(colour = "black", size = 11),
    plot.margin      = margin(6, 8, 6, 6)
  )
p_resp
