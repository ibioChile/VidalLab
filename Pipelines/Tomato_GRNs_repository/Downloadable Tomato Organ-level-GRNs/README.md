# 🍅 Downloadable Organ-level GRNs for *Solanum lycopersicum*

Welcome to the **Tomato Organ-level Gene Regulatory Networks (GRNs)** resource!  
This folder provides ready-to-use, community-accessible GRNs for *Solanum lycopersicum* (tomato), generated from high-quality transcriptomic and epigenomic datasets.

---

## 📂 Available Files

You’ll find one GRN file per organ:

- `Seed_GRN.txt`
- `Root_GRN.txt`
- `Leaf_GRN.txt`
- `Flower_GRN.txt`
- `Fruit_GRN.txt`

Each file is a **tab-separated `.txt` file** with the following structure:

| **Column** | **Description** |
|------------|-----------------|
| `TF`       | Transcription Factor |
| `Target`   | Target Gene |
| `Evidence` | Supporting regulatory evidence: may include one or more of the following → `GENIE3`, `Coex`, `OCS`, `Promo` |

> ⚠️ All edges originate from **GENIE3** predictions and are optionally supported by:
> - **Coex**: Co-expression  
> - **OCS**: Open Chromatin Sites (e.g., ATAC-seq, DNase-seq)  
> - **Promo**: Promoter motif matches

---

## 🧪 How to Use

These networks are compatible with a range of bioinformatics and visualization tools:

- 🧬 **Cytoscape** – for interactive network exploration  
- 📊 **R / Python** – for filtering, graph modeling, and custom analyses  
- 📑 **Excel / LibreOffice** – for simple browsing or filtering  

---

## 📚 Citation

If you use these GRNs in your work, please cite our article:

**Organ-level Gene Regulatory Network models enable the identification of central transcription factors in Solanum lycopersicum**
(2025-04-01)  doi: https://doi.org/10.1101/2025.03.26.645553
**[https://github.com/ibioChile/VidalLab](https://github.com/ibioChile/VidalLab)**

---

Thanks for visiting — we hope this resource accelerates your plant regulatory genomics research! 🌱  
