library(tibble)

comparisons <- tribble(
  ~label,               ~experiment, ~control,  ~test,

  # Exp1 core comparisons
  "Guy11_vs_mst7",      "Exp1",      "Guy11",   "mst7",
  "Guy11_vs_MST7WT",    "Exp1",      "Guy11",   "MST7WT",
  "MST7WT_vs_mst7",     "Exp1",      "MST7WT",  "mst7",
  "MST7WT_vs_PM",       "Exp1",      "MST7WT",  "MST7PM",
  "MST7WT_vs_PD",       "Exp1",      "MST7WT",  "MST7PD",
  "PM_vs_PD",           "Exp1",      "MST7PM",  "MST7PD",

  # Exp2 core comparison
  "Guy11M_vs_pmk1",     "Exp2",      "Guy11M",  "pmk1",

  # Cross‑experiment comparisons
  "mst7_vs_pmk1",       "Cross",     "mst7",    "pmk1",
  "Guy11_vs_Guy11M",    "Cross",     "Guy11",   "Guy11M"
)
