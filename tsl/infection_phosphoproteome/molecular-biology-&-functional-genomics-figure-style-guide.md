# Molecular Biology &amp; Functional Genomics Figure Style Guide

*For Publication-Level Figures*

---

## **1. General Principles**

- **Purpose**: Figures must prioritize **clarity**, **reproducibility**, and **accessibility**.
- **Audience**: Researchers, reviewers, and readers (including those with color vision deficiencies).
- **Avoid AI-generated aesthetics**: No overly smooth gradients, unnatural symmetry, no subheadings, or "perfect" AI-like rendering.

---

## **2. Color Guidelines**

### **2.1 Colorblind-Friendly Palettes**

- Use **scientifically validated palettes** to ensure accessibility.
- **Recommended tools**:
  - [ColorBrewer](https://colorbrewer2.org/)
  - [Viridis](https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html)

### **2.2 Palette Types**


| Type        | Use Case                         | Recommended Palettes              |
| ----------- | -------------------------------- | --------------------------------- |
| Sequential  | Continuous data (e.g., heatmaps) | Viridis, Plasma, Cividis          |
| Qualitative | Categorical data (e.g., groups)  | Okabe-Ito, Tol’s bright palette   |
| Diverging   | Data with a critical midpoint    | Coolwarm, RdBu (use with caution) |


### **2.3 Colors to Avoid**

- ❌ Red-green combinations (most common form of color blindness).
- ❌ Low-contrast colors (e.g., light yellow on white).
- ❌ Overly bright or neon colors.

---

## **3. Fonts and Text**

### **3.1 Typeface**

- **Labels and axes**: Sans-serif (e.g., Arial, Helvetica).
- **Legends and captions**: Match journal requirements (often serif, e.g., Times New Roman).

### **3.2 Font Sizes**


| Element       | Size (pt) | Style   |
| ------------- | --------- | ------- |
| Axis labels   | 10–12     | Bold    |
| Tick labels   | 8–10      | Regular |
| Legends       | 9–11      | Regular |
| Figure titles | 12–14     | Bold    |


### **3.3 Consistency**

- Use the **same font family** throughout the figure.
- Avoid mixing serif and sans-serif in the same figure.

---

## **4. Line and Shape Styles**

### **4.1 Lines**

- **Data lines**: Solid, 1–2 pt width.
- **Reference lines**: Dashed or dotted, 0.5–1 pt width.
- **Grid lines**: Light gray, dashed, 0.5 pt width (if used).

### **4.2 Markers**

- Use **distinct shapes** for multiple datasets:
  - Circle (○)
  - Square (□)
  - Triangle (△)
  - Diamond (◇)
- **Size**: 4–8 pt (adjust for visibility).

---

## **5. Layout and Composition**

### **5.1 Margins and Spacing**

- Leave **10–15% white space** around the figure.
- **Multi-panel figures**: Maintain consistent spacing (e.g., 0.5 cm between panels).

### **5.2 Aspect Ratios**

- **Scatter plots**: Square (1:1).
- **Bar charts**: 4:3 or 16:9.
- **Heatmaps**: Adjust to data dimensions (avoid distortion).

### **5.3 Alignment**

- Align **axes, labels, and legends** with the figure edges.
- Group related panels **visually** (e.g., shared axes for multi-panel figures).

---

## **6. File Formats**


| Format   | Use Case        | Notes                               |
| -------- | --------------- | ----------------------------------- |
| `.svg`   | Vector graphics | Preferred for scalability.          |
| `.pdf`   | Vector graphics | Embed fonts for journal submission. |
| `.tiff`  | Raster graphics | 300+ DPI for print.                 |
| `.png`   | Raster graphics | 150+ DPI for web.                   |
| ❌ `.jpg` | **Avoid**       | Lossy compression degrades quality. |


---

## **7. Labels and Annotations**

### **7.1 Axis Labels**

- Include **units** (e.g., "Concentration (µM)").
- Use **clear, descriptive language** (e.g., "Gene Expression (FPKM)").

### **7.2 Legends**

- Place **inside the figure** if space allows.
- Otherwise, place **below or to the side** of the figure.
- Match **marker shapes/colors** to the data.

### **7.3 Error Bars**

- Clearly label (e.g., "± SEM", "95% CI").
- Use **whiskers** for box plots (label as "1.5× IQR").

---

## **8. Data Representation**

### **8.1 Bar Charts**

- Use for **categorical comparisons**.
- Avoid **3D effects** (they distort perception).
- **Error bars**: Include if applicable.

### **8.2 Heatmaps**

- Use **hierarchical clustering** sparingly (can be hard to interpret).
- Always include a **color scale** with a clear label.

### **8.3 Scatter Plots**

- Add **trend lines** only if statistically justified.
- Use **transparency** for overlapping points.

### **8.4 Box Plots**

- Show **outliers** as individual points.
- Label **whiskers** (e.g., "1.5× IQR").

---

## **9. Accessibility**

### **9.1 Alt Text**

- Provide **descriptive text** for each figure:
  - Example: "Bar chart showing gene expression levels in Condition A (mean = 5.2) vs. Condition B (mean = 3.1). Error bars represent ± SEM."

### **9.2 Patterns and Textures**

- Use **textures** (e.g., stripes, dots) for colorblind-safe distinctions in grayscale figures.

---

## **10. Tools and Software**

### **10.1 Recommended Tools**


| Tool                | Use Case                          |
| ------------------- | --------------------------------- |
| Python (Matplotlib) | Customizable, scriptable figures. |
| Python (Seaborn)    | Statistical data visualization.   |
| Python (Plotly)     | Interactive figures (for web).    |
| R (ggplot2)         | Publication-quality plots.        |
| Inkscape            | Vector graphics editing.          |
| Adobe Illustrator   | Advanced vector editing.          |


### **10.2 Tools to Avoid**

- ❌ Default styles from Excel or PowerPoint (often non-compliant with journal standards).

---

## **11. Journal-Specific Requirements**

- **Always check** the journal’s *Instructions for Authors* for:
  - Maximum figure dimensions.
  - File format restrictions (e.g., some journals require `.tiff` for print).
  - Color vs. grayscale policies (some journals charge for color figures).

---

## **12. Pre-Submission Checklist**

- [ ] All fonts are **embedded** (for `.pdf`/`.svg`).
- [ ] Figure resolution meets journal requirements (e.g., 300 DPI for print).
- [ ] Colors are **colorblind-friendly** (test with [Color Oracle](https://colororacle.org/)).
- [ ] All axes, labels, and legends are **clearly readable** at submission size.
- [ ] File format is **journal-compliant** (e.g., `.tiff`, `.pdf`, `.svg`).
- [ ] **Alt text** is provided for each figure.
- [ ] Figure legends are **included** (either in the figure or manuscript).

---

## **13. Code Snippets for Compliant Figures**

### **13.1 Matplotlib (Python)**

```python
import matplotlib.pyplot as plt

# Set style for publication
plt.style.use('seaborn-v0_8-poster')  # Clean, scalable style
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']

# Example plot
fig, ax = plt.subplots(figsize=(6, 4), dpi=300)
ax.plot([1, 2, 3], [4, 5, 6], color='#1f77b4', linewidth=1.5)
ax.set_xlabel('X Axis (Units)', fontsize=12, fontweight='bold')
ax.set_ylabel('Y Axis (Units)', fontsize=12, fontweight='bold')
ax.tick_params(labelsize=10)

# Save as vector graphic
fig.savefig('figure.svg', format='svg', bbox_inches='tight')
```

### **13.2 ggplot2 (R)**

```r
library(ggplot2)

# Set theme for publication
theme_pub <- theme(
  text = element_text(family = "Arial", size = 12),
  axis.text = element_text(size = 10),
  axis.title = element_text(size = 12, face = "bold")
)

# Example plot
p <- ggplot(mtcars, aes(x = wt, y = mpg)) +
  geom_point(color = "#1f77b4", size = 3) +
  labs(x = "Weight (1000 lbs)", y = "Miles per Gallon") +
  theme_pub

# Save as PDF (vector)
ggsave("figure.pdf", p, width = 6, height = 4, dpi = 300)
```

---

## **14. Examples and Templates**

### **14.1 Before and After**

- **Before**: Default Excel bar chart with red-green colors and 3D effects.
- **After**: Flat 2D bar chart with Okabe-Ito palette, embedded fonts, and clear labels.

*(Add visual examples here as needed.)*

---

## **15. References**

- [ColorBrewer](https://colorbrewer2.org/)
- [Viridis: A Colorblind-Friendly Palette](https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html)
- [Okabe-Ito Color Palette](https://jfly.uni-koeln.de/color/)
- [Color Oracle (Color Blindness Simulator)](https://colororacle.org/)

---

*Last updated: \[Insert Date\]*

*Author: \[Your Name\]*