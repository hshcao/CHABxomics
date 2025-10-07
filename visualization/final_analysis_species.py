# -*- coding: utf-8 -*-
"""
Author:         Ziyue Yin
Last Modified:  2025-09-23
Description:
    This script analyzes the response of individual algal species to temperature peaks in Lake Taihu,
    classifying species as either 'growing' or 'suppressed' based on their abundance changes during
    temperature peak months. The analysis generates publication-ready visualizations and statistical summaries.

Output Structure:
    PDF_species/
    ├── All/                    # Full year analysis
    │   ├── growing/           # Species showing growth at temperature peak
    │   └── suppressed/        # Species showing suppression at temperature peak
    ├── 6-10/                  # June-October subset analysis
    │   ├── growing/
    │   └── suppressed/
    ├── class_counts_all.csv   # Summary statistics for full year
    ├── class_counts_6-10.csv  # Summary statistics for June-October
    └── SingleSpecies_*_5x3_A4_portrait.pdf  # Collage plots (5×3 grid)
"""

import os
import re
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings('ignore')

# ----------------------------- Configuration -----------------------------
# Data files and worksheet names
FILE_XLS      = '2019年藻类和水温_英文学名.xlsx'
SHEET_ALGAE   = '藻类计数'
SHEET_TEMP    = '水温'

# Classification parameters
ZERO_CUTOFF         = 1.0   # Abundance values ≤ this threshold treated as missing for classification
USE_MAJORITY_VOTE   = False # True: site majority voting; False: cross-site mean
SITE_MIN_K          = 3     # Minimum vote threshold for majority voting

# 6-10 species filtering: require at least one site with complete data (no missing/zero values)
REQUIRE_AT_LEAST_ONE_COMPLETE_SITE_6_10 = True

# Plotting appearance
plt.rcParams['font.sans-serif'] = ['Arial Unicode MS', 'SimHei']
plt.rcParams['axes.unicode_minus'] = False

# Deep color palette (excluding light blue/gray)
DEEP_COLORS = ['#1f77b4','#2ca02c','#9467bd','#8c564b','#17becf',
               '#7f7f7f','#bcbd22','#393b79','#637939','#5254a3',
               '#6b6ecf','#8ca252']
EPS = 1e-6

# ----------------------------- Helpers ----------------------------------
def dt_for_month(period_m):
    """Create datetime for month display: use 15th for months 1-11, 14th for December.
    This facilitates x-axis labels showing only MM-DD without year."""
    y, m = period_m.year, period_m.month
    d = 14 if m == 12 else 15
    return datetime(y, m, d)

def strip_2019_mmdd(dt: datetime) -> str:
    return dt.strftime('%m-%d')

def natural_site_sort_key(s: str):
    m = re.search(r'(\d+)', str(s))
    return (re.sub(r'\d+', '', str(s)), int(m.group(1)) if m else -1)

KNOWN_PHYLA = {
    'cyanophyta','cyanobacteria','chlorophyta','bacillariophyta','diatomea','diatoms',
    'euglenophyta','chrysophyta','cryptophyta','dinophyta','haptophyta',
    'rhodophyta','xanthophyta','ochrophyta','charophyta','myzozoa'
}
def looks_like_phylum(val: str) -> bool:
    if not isinstance(val, str):
        return False
    s = val.strip().lower()
    return (s in KNOWN_PHYLA) or s.endswith('phyta') or s.endswith('bacteria') or s in {'diatoms','green algae'}

def ensure_dirs():
    for span in ['All', '6-10']:
        for cls in ['growing', 'suppressed']:
            os.makedirs(f'./PDF_species/{span}/{cls}', exist_ok=True)

# ----------------------------- Load data ---------------------------------
xls = pd.ExcelFile(FILE_XLS)
df_algae = xls.parse(SHEET_ALGAE)   # Species abundance (sites as columns)
df_temp  = xls.parse(SHEET_TEMP)    # Water temperature (month × sites)

# Temperature table: construct YM (month period) and standardize site names
df_temp['YM'] = pd.to_datetime(
    df_temp[['年','月']].rename(columns={'年':'year','月':'month'}).assign(day=1)
).dt.to_period('M')

df_temp_long = (
    df_temp.assign(Site=df_temp['点位'].astype(str).str.replace('S','St', regex=False))
           .rename(columns={'水温':'Temp'})[['YM','Site','Temp']]
)

# --------------------- Parse Phylum / Genus / Species --------------------
# Candidate column names
CAND_PHYLUM   = ['门','Phylum','phylum','phyla','Phyla','phylum_en','Phylum_EN','Phylum (English)']
CAND_GENUS    = ['属','Genus','genus']
CAND_SPECIES  = ['种','Species','species','species_name','物种','species_en']

def choose_col(df, cand):
    for c in cand:
        if c in df.columns:
            return c
    return None

phylum_col  = choose_col(df_algae, CAND_PHYLUM)
genus_col   = choose_col(df_algae, CAND_GENUS)
species_col = choose_col(df_algae, CAND_SPECIES)

def row_to_taxa(row):
    phy = row[phylum_col] if phylum_col else None
    gen = row[genus_col] if genus_col else None
    spf = row[species_col] if species_col else None  # Original "Species" text

    # If no "Phylum" column but "Genus" column actually contains phylum (e.g., Chlorophyta)
    used_genus_as_phylum = False
    if phy is None and isinstance(gen, str) and looks_like_phylum(gen):
        phy = gen
        used_genus_as_phylum = True

    # Construct full species name and provide genus name when possible
    species_name = None
    if isinstance(spf, str) and len(spf.strip()) > 0:
        s = spf.strip()
        if ' ' in s:  # Already a "binomial name"
            species_name = s
            gen2 = s.split()[0]
            if not used_genus_as_phylum:
                gen = gen2
        else:
            # Single word: if gen is not a phylum, concatenate with gen
            if (not used_genus_as_phylum) and isinstance(gen, str) and len(gen.strip())>0 and not looks_like_phylum(gen):
                species_name = f"{gen.strip()} {s}"
            else:
                species_name = s
    if species_name is None and isinstance(gen, str) and len(gen.strip())>0 and not looks_like_phylum(gen):
        species_name = gen.strip()
    if not isinstance(species_name, str) or len(species_name.strip()) == 0:
        species_name = 'Unknown species'

    phy = phy if isinstance(phy, str) and len(phy.strip())>0 else 'Unknown'
    return pd.Series({'Phylum': phy.strip(), 'Genus': (gen if isinstance(gen, str) else ''), 'SpeciesName': species_name.strip()})

taxa = df_algae.apply(row_to_taxa, axis=1)
df_algae = pd.concat([df_algae, taxa], axis=1)

# ----------------------- Build long table & merge T ----------------------
# Time column and YM
time_col = '时间' if '时间' in df_algae.columns else None
if time_col is None:
    raise ValueError("Time column not found (e.g., '时间'). Please check source table.")
df_algae['YM'] = pd.to_datetime(df_algae[time_col]).dt.to_period('M')

# Site columns: starting with 'St' or 'S' + number
site_cols = []
for c in df_algae.columns:
    s = str(c)
    if re.match(r'^(St|S)\d+$', s, flags=re.IGNORECASE):
        if s.startswith('S') and not s.startswith('St'):
            # Standardize to St1 format
            df_algae.rename(columns={c: 'St'+s[1:]}, inplace=True)
            site_cols.append('St'+s[1:])
        else:
            site_cols.append(s)
site_cols = sorted(list(set(site_cols)), key=natural_site_sort_key)
if not site_cols:
    raise ValueError("Site columns not found (format: 'St1','St2',... or 'S1','S2',...).")

id_vars = ['Phylum','Genus','SpeciesName','YM', time_col]
df_abund = df_algae.melt(id_vars=id_vars, value_vars=site_cols,
                         var_name='Site', value_name='Abundance')

df_long = df_abund.merge(df_temp_long, on=['YM','Site'], how='left')
SITES = sorted(df_long['Site'].dropna().unique(), key=natural_site_sort_key)

# ----------------------- 6-10 subset helpers ----------------------------
months_6_10 = pd.period_range('2019-06', '2019-10', freq='M')

def site_complete_6_10(df_species, site):
    """Check if a site has complete data (no missing/zero values) for June-October period."""
    sub = df_species[df_species['Site'] == site]
    for m in months_6_10:
        row = sub[sub['YM'] == m]
        if row.empty or pd.isna(row['Abundance'].iloc[0]) or (row['Abundance'].iloc[0] <= 0):
            return False
    return True

def species_has_any_complete_site_6_10(species_name):
    df_sp = df_long[df_long['SpeciesName'] == species_name]
    for s in SITES:
        if site_complete_6_10(df_sp, s):
            return True
    return False

# -------------------- Classification on FULL timeline --------------------
def find_prev_with_data(series_by_month: pd.Series, month: pd.Period, k_back=2):
    """Search backwards up to k_back months to find the nearest month with valid data.
    Returns None if no valid data found."""
    m = month
    for _ in range(k_back):
        m = (m - 1)
        if m in series_by_month.index and pd.notna(series_by_month.loc[m]):
            return m
    return None

def classify_full_series(df_species_full: pd.DataFrame,
                         zero_cutoff: float = ZERO_CUTOFF,
                         majority_vote: bool = USE_MAJORITY_VOTE,
                         site_min_k: int = SITE_MIN_K):
    """
    Classify species as growing/suppressed based on full timeline:
      • Temperature peak month: Find month with highest average temperature across sites
      • Abundance: Values ≤ zero_cutoff treated as missing (NaN)
      • v_prev: If no data in previous month, search backwards up to 2 months for recent data
      • If majority_vote=True: Vote by site, classify as growing if ≥site_min_k sites show increase
    Returns: (classification, peak_month, info_dict)
    """
    all_months_full = sorted(df_species_full['YM'].dropna().unique())

    # Temperature mean (across sites)
    t_mean = (df_species_full.groupby('YM')['Temp']
              .mean().reindex(all_months_full))
    if t_mean.isna().all():
        return 'suppressed', None, {'reason': 'all-temp-NaN'}

    try:
        peak_idx = int(np.nanargmax(t_mean.values))
    except ValueError:
        return 'suppressed', None, {'reason': 'no-temp-peaks'}
    peak_month = all_months_full[peak_idx]

    # Abundance: threshold filtering (for classification only, does not affect plotting)
    abund = df_species_full.copy()
    abund.loc[abund['Abundance'] <= zero_cutoff, 'Abundance'] = np.nan

    if majority_vote:
        votes = 0
        total_sites = 0
        for site, g in abund.groupby('Site'):
            s = g.set_index('YM')['Abundance'].reindex(all_months_full)
            v_now = s.loc[peak_month]
            prev_m = find_prev_with_data(s, peak_month, k_back=2)
            if pd.notna(v_now) and prev_m is not None and pd.notna(s.loc[prev_m]):
                total_sites += 1
                if v_now > s.loc[prev_m]:
                    votes += 1
        if total_sites == 0:
            return 'suppressed', peak_month, {'reason': 'no-usable-sites'}
        cls = 'growing' if votes >= max(site_min_k, total_sites//2 + 1) else 'suppressed'
        return cls, peak_month, {'votes': votes, 'total_sites': total_sites}

    # Default: cross-site mean
    a_mean = (abund.groupby('YM')['Abundance']
              .mean().reindex(all_months_full))
    v_now = a_mean.loc[peak_month]
    prev_m = find_prev_with_data(a_mean, peak_month, k_back=2)
    if pd.notna(v_now) and prev_m is not None and pd.notna(a_mean.loc[prev_m]):
        cls = 'growing' if v_now > a_mean.loc[prev_m] else 'suppressed'
        return cls, peak_month, {'prev': str(prev_m)}
    else:
        # Insufficient data, conservatively classify as suppressed
        return 'suppressed', peak_month, {'reason': 'insufficient-abundance'}

def classify_species(species_name: str):
    df_sp_full = df_long[df_long['SpeciesName'] == species_name]
    return classify_full_series(df_sp_full)

# ------------------------------- Plotting --------------------------------
def plot_one_species(species_name: str, phylum: str, time_span='All', cls_hint=None,
                     save_base='./PDF_species'):
    """
    Plot individual species: temperature (red line) + log10(abundance) for each site.
    Title format: "Phylum: Species Name".
    
    Parameters:
    -----------
    time_span : str
        'All' for full year or '6-10' for June-October subset
    cls_hint : str, optional
        Pre-determined classification (growing/suppressed) to avoid recalculation
    """
    df_sp = df_long[df_long['SpeciesName'] == species_name].copy()
    if df_sp.empty:
        return None

    # Subset months
    if time_span == '6-10':
        df_sp = df_sp[df_sp['YM'].isin(months_6_10)]
    all_months = sorted(df_sp['YM'].dropna().unique())
    if len(all_months) == 0:
        return None

    # Classification label: use full timeline result
    cls = cls_hint
    if cls is None:
        cls, _pm, _info = classify_full_series(df_long[df_long['SpeciesName']==species_name])

    # x 轴
    dates = [dt_for_month(m) for m in all_months]
    x = list(range(len(all_months)))

    out_dir = os.path.join(save_base, time_span, cls)

    # Temperature reference (cross-site mean)
    t_ref = df_sp.groupby('YM')['Temp'].mean().reindex(all_months).values

    # Create plot
    figsize = (7.5, 7) if time_span == 'All' else (7.5, 7)
    fig, ax_t = plt.subplots(figsize=figsize)
    ax_a = ax_t.twinx()

    # Temperature (red line)
    ax_t.plot(x, t_ref, color='red', linewidth=2, marker='o', markersize=4, alpha=0.75)
    ax_t.set_ylabel('Temperature (°C)', color='red')
    ax_t.tick_params(axis='y', labelcolor='black')

    # Site abundance (log10)
    plotted = 0
    for site in SITES:
        sub = df_sp[df_sp['Site'] == site]
        if sub.empty:
            continue
        if time_span == '6-10' and REQUIRE_AT_LEAST_ONE_COMPLETE_SITE_6_10:
            # Only filter out incomplete site curves in 6-10 plots
            if not site_complete_6_10(df_long[df_long['SpeciesName'] == species_name], site):
                continue

        vals = []
        for m in all_months:
            r = sub[sub['YM'] == m]
            v = r['Abundance'].iloc[0] if not r.empty else 0.0
            vals.append(v)
        vals = np.array(vals, dtype=float)
        vals_log = np.log10(vals + EPS)

        ax_a.plot(x, vals_log, color=DEEP_COLORS[plotted % len(DEEP_COLORS)],
                  linewidth=2, marker='o', markersize=4, label=str(site))
        plotted += 1

    ax_a.set_ylabel(r'$\log_{10}(\mathrm{Abundance})$')
    ax_a.tick_params(axis='y', labelcolor='black')

    # x-axis labels (remove 2019)
    ax_t.set_xticks(x)
    ax_t.set_xticklabels([strip_2019_mmdd(d) for d in dates], rotation=0)

    # Title: "Phylum: Species Name"
    title = f"{phylum}: {species_name}  ({'Full Year' if time_span=='All' else 'June - October'})"
    ax_t.set_title(title, fontsize=14, fontweight='bold', pad=10)

    ax_t.grid(True, linestyle='--', alpha=0.3)
    lines2, labels2 = ax_a.get_legend_handles_labels()
    if labels2:
        # Keep legend in one row, adjust font size based on number of sites
        ncol = len(labels2)  # All legends in one row
        fontsize = max(6, 12 - len(labels2))  # Adjust font size based on number of sites
        plt.legend(lines2, labels2, loc='upper center', bbox_to_anchor=(0.5, -0.06),
                   ncol=ncol, frameon=False, fontsize=fontsize)

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.16)

    # Filename: replace colons with "_"
    safe_name = f"{phylum}_{species_name}.pdf"
    out_path = os.path.join(out_dir, safe_name)
    plt.savefig(out_path, format='pdf', dpi=300, bbox_inches='tight')
    plt.close(fig)

    print(f"[{time_span:>4}] {cls:10s} -> {out_path}")
    return cls

# --------------------------- Batch & Summaries ---------------------------
def summarize_classes(df_res: pd.DataFrame, csv_path: str):
    if df_res is None or df_res.empty:
        out = pd.DataFrame(columns=['phylum','n_growing','n_suppressed','n_total'])
        out.to_csv(csv_path, index=False)
        return out
    tmp = (df_res.groupby(['phylum','class'])
           .agg(n=('species','nunique'))
           .reset_index())
    piv = (tmp.pivot(index='phylum', columns='class', values='n')
              .fillna(0).astype(int))
    
    if 'growing' not in piv.columns:    piv['growing'] = 0
    if 'suppressed' not in piv.columns: piv['suppressed'] = 0
    piv = piv.rename(columns={'growing':'n_growing','suppressed':'n_suppressed'})
    piv['n_total'] = piv['n_growing'] + piv['n_suppressed']
    out = piv.reset_index().sort_values(['n_growing','n_total'], ascending=[False, False])
    out.to_csv(csv_path, index=False)
    return out

def batch_run(time_span: str, df_sp_list: pd.DataFrame):
    """
    Plot and classify species from given species list.
    Note: Classification based on full timeline completed once; plotting filters months by time_span.
    """
    records = []
    for _, r in df_sp_list.iterrows():
        phy, spname = r['Phylum'], r['SpeciesName']
        cls, _pm, _info = classify_species(spname)
        # Plot (using pre-determined cls to avoid recalculation)
        _ = plot_one_species(spname, phy, time_span=time_span, cls_hint=cls)
        records.append((phy, spname, cls))
    return pd.DataFrame(records, columns=['phylum','species','class'])

# =================== 5×3 Collage Export ===================
# Dependencies:
#   Vector mode: pip install pymupdf
#   Raster fallback: pip install pillow pymupdf
from pathlib import Path
from io import BytesIO

def _parse_meta_from_filename(pdf_path: Path):
    """Parse filename like 'Phylum_Species name.pdf' → (phylum, genus, species)."""
    name = pdf_path.stem
    parts = name.split('_', 1)
    if len(parts) < 2:
        return ('Unknown', 'Unknown', name)
    phylum, species = parts[0], parts[1]
    genus = species.split()[0] if species.strip() else 'Unknown'
    return (phylum, genus, species)

def _collect_items(root_dir: Path, cls: str, time_span: str = 'All'):
    """
    Collect individual PDF files from PDF_species/<time_span>/<cls>/ and mark whether
    corresponding figures exist in 6-10 directory.
    Returns: list of dict(path, phylum, genus, species, has610)
    """
    dir_target = root_dir / time_span / cls
    dir_610 = root_dir / '6-10' / cls
    items = []
    for p in sorted(dir_target.glob('*.pdf')):
        phylum, genus, species = _parse_meta_from_filename(p)
        has610 = (dir_610 / f"{phylum}_{species}.pdf").exists()
        items.append({
            'path': p, 'phylum': phylum, 'genus': genus, 'species': species,
            'has610': has610
        })
    return items

def _select_15(items, limit=15):
    """
    Select ≤15 items with rules:
      1) Prefer has610=True;
      2) Try to cover more phyla/genera;
      3) Finally sort by (phylum, genus, species) to keep same phylum/genus adjacent.
    """
    if len(items) <= limit:
        return sorted(items, key=lambda d: (d['phylum'].lower(), d['genus'].lower(), d['species'].lower()))

    preferred = [it for it in items if it['has610']]
    nonpref   = [it for it in items if not it['has610']]

    selected, chosen = [], set()

    # A. At least 1 from each phylum in preferred
    by_phy = {}
    for it in preferred:
        by_phy.setdefault(it['phylum'], []).append(it)
    for phy in sorted(by_phy.keys(), key=lambda s: s.lower()):
        for it in by_phy[phy]:
            key = it['path']
            if key not in chosen:
                selected.append(it); chosen.add(key)
                break
        if len(selected) >= limit: break

    # B. Cover more genera in preferred
    seen_pg = {(it['phylum'], it['genus']) for it in selected}
    for it in preferred:
        if len(selected) >= limit: break
        key = it['path']
        if key in chosen: continue
        if (it['phylum'], it['genus']) not in seen_pg:
            selected.append(it); chosen.add(key); seen_pg.add((it['phylum'], it['genus']))

    # C. Fill remaining from preferred
    for it in preferred:
        if len(selected) >= limit: break
        key = it['path']
        if key not in chosen:
            selected.append(it); chosen.add(key)

    # D. Cover genera from nonpref → fill remaining
    for it in nonpref:
        if len(selected) >= limit: break
        key = it['path']
        if (it['phylum'], it['genus']) not in seen_pg:
            selected.append(it); chosen.add(key); seen_pg.add((it['phylum'], it['genus']))
    for it in nonpref:
        if len(selected) >= limit: break
        key = it['path']
        if key not in chosen:
            selected.append(it); chosen.add(key)

    selected = sorted(selected, key=lambda d: (d['phylum'].lower(), d['genus'].lower(), d['species'].lower()))
    return selected[:limit]

def _draw_collage_a4_portrait_vector(pdf_paths, out_pdf_path, rows=5, cols=3,
                                     margin_pt=18, gap_pt=6):
    """
    Vector collage: Embed first page of each "small PDF" as vector in A4 portrait page, 5×3 grid.
    Best clarity; requires pymupdf (fitz).
    """
    try:
        import fitz  # PyMuPDF
    except Exception as e:
        raise RuntimeError("PyMuPDF required: pip install pymupdf") from e

    doc = fitz.open()
    page_rect = fitz.paper_rect("a4")  # portrait
    page = doc.new_page(width=page_rect.width, height=page_rect.height)
    W, H = page_rect.width, page_rect.height

    tile_w = (W - 2*margin_pt - (cols - 1)*gap_pt) / cols
    tile_h = (H - 2*margin_pt - (rows - 1)*gap_pt) / rows

    for idx in range(rows * cols):
        if idx >= len(pdf_paths):
            break
        r = idx // cols
        c = idx % cols
        x0 = margin_pt + c * (tile_w + gap_pt)
        y0 = margin_pt + r * (tile_h + gap_pt)

        with fitz.open(str(pdf_paths[idx])) as src:
            sp = src[0]
            sw, sh = sp.rect.width, sp.rect.height
            s_ar = sw / max(1e-6, sh)
            c_ar = tile_w / tile_h
            if c_ar > s_ar:
                draw_h = tile_h
                draw_w = draw_h * s_ar
            else:
                draw_w = tile_w
                draw_h = draw_w / s_ar
            dx0 = x0 + (tile_w - draw_w) / 2.0
            dy0 = y0 + (tile_h - draw_h) / 2.0
            dest = fitz.Rect(dx0, dy0, dx0 + draw_w, dy0 + draw_h)
            page.show_pdf_page(dest, src, 0)

    doc.save(str(out_pdf_path))
    doc.close()

def _draw_collage_a4_portrait_raster(pdf_paths, out_pdf_path, rows=5, cols=3,
                                     dpi=300, margin_px=60, gap_px=20):
    """
    Fallback: High DPI rasterization → composite into large image → save as PDF.
    Clear and controllable, but final result is bitmap embedding; requires Pillow + PyMuPDF.
    """
    try:
        from PIL import Image
    except Exception as e:
        raise RuntimeError("Pillow required: pip install pillow") from e
    try:
        import fitz
    except Exception as e:
        raise RuntimeError("PyMuPDF required: pip install pymupdf") from e

    # A4 portrait pixel dimensions
    A4_W_IN, A4_H_IN = 8.27, 11.69
    W = int(round(A4_W_IN * dpi))
    H = int(round(A4_H_IN * dpi))

    canvas = Image.new("RGB", (W, H), "white")

    tile_w = (W - 2*margin_px - (cols - 1)*gap_px) // cols
    tile_h = (H - 2*margin_px - (rows - 1)*gap_px) // rows

    for idx in range(rows * cols):
        if idx >= len(pdf_paths):
            break
        r = idx // cols
        c = idx % cols
        x0 = margin_px + c * (tile_w + gap_px)
        y0 = margin_px + r * (tile_h + gap_px)

        with fitz.open(str(pdf_paths[idx])) as src:
            sp = src[0]
            sw, sh = sp.rect.width, sp.rect.height
            s_ar = sw / max(1e-6, sh)
            c_ar = tile_w / tile_h
            if c_ar > s_ar:
                target_h = tile_h
                target_w = int(round(target_h * s_ar))
            else:
                target_w = tile_w
                target_h = int(round(target_w / s_ar))
            scale = min(target_w / sw, target_h / sh)
            mat = fitz.Matrix(scale, scale)
            pix = sp.get_pixmap(matrix=mat, alpha=False)
            img = Image.open(BytesIO(pix.tobytes("png")))
            dx = x0 + (tile_w - img.width) // 2
            dy = y0 + (tile_h - img.height) // 2
            canvas.paste(img, (dx, dy))

    # Save as PDF (single page)
    canvas.save(str(out_pdf_path), "PDF", resolution=dpi)

def build_collages_portrait(base_dir="./PDF_species", time_span="All", mode="vector"):
    """
    Build collage PDFs (A4 portrait, 5×3):
      For specified time_span:
      1) growing    → SingleSpecies_{time_span}_growing_5x3_A4_portrait.pdf
      2) suppressed → SingleSpecies_{time_span}_suppressed_5x3_A4_portrait.pdf
    mode:
      - "vector" (default): Vector embedding, clearest.
      - "raster": High DPI rasterization then composite.
    """
    root = Path(base_dir)
    root.mkdir(parents=True, exist_ok=True)

    for cls in ["growing", "suppressed"]:
        items = _collect_items(root, cls, time_span)
        if not items:
            print(f"[WARN] No individual PDFs found in {root/time_span/cls}, skipping collage.")
            continue
        selected = _select_15(items, limit=15)
        pdf_list = [it['path'] for it in selected]

        out_name = f"SingleSpecies_{time_span}_{cls}_5x3_A4_portrait.pdf"
        out_path = root / out_name

        if mode == "raster":
            _draw_collage_a4_portrait_raster(pdf_list, out_path, rows=5, cols=3,
                                             dpi=300, margin_px=60, gap_px=20)
        else:
            _draw_collage_a4_portrait_vector(pdf_list, out_path, rows=5, cols=3,
                                             margin_pt=18, gap_pt=6)

        print(f"[Collage Complete] {time_span:>4} {cls:10s} -> {out_path}")

def build_all_collages_portrait(base_dir="./PDF_species", mode="vector"):
    """
    Build two collage PDFs (A4 portrait, 5×3):
      1) All/growing    → SingleSpecies_All_growing_5x3_A4_portrait.pdf
      2) All/suppressed → SingleSpecies_All_suppressed_5x3_A4_portrait.pdf
    mode:
      - "vector" (default): Vector embedding, clearest.
      - "raster": High DPI rasterization then composite.
    """
    build_collages_portrait(base_dir, "All", mode)


# ------------------------------- Main ------------------------------------
def main():
    ensure_dirs()

    # All species list (by phylum-species)
    all_species = (df_long[['Phylum','SpeciesName']]
                   .drop_duplicates()
                   .reset_index(drop=True))

    # 6-10: Filter species with "at least one complete site" (can be changed to no filtering if not needed)
    if REQUIRE_AT_LEAST_ONE_COMPLETE_SITE_6_10:
        sp_6_10 = all_species[
            all_species['SpeciesName'].apply(species_has_any_complete_site_6_10)
        ].reset_index(drop=True)
    else:
        sp_6_10 = all_species.copy()

    print("\n==== Starting generation: Full Year (All) ====")
    res_all = batch_run('All', all_species)

    print("\n==== Starting generation: June-October (6-10) ====")
    res_610 = batch_run('6-10', sp_6_10)

    # Summary statistics (growing & suppressed)
    sum_all = summarize_classes(res_all, './PDF_species/class_counts_all.csv')
    sum_610 = summarize_classes(res_610, './PDF_species/class_counts_6-10.csv')

    print("\n==== Statistical Results (for manuscript) ====")
    print("【Full Year All】Number of species per phylum showing \"biomass increase vs decrease/no increase\" at temperature peak month:")
    print(sum_all.to_string(index=False))
    print("\n【June-October】Number of species per phylum showing \"biomass increase vs decrease/no increase\" at temperature peak month:")
    print(sum_610.to_string(index=False))
    print("\nAll completed.")
    
    # === Build 5×3 collage PDFs (A4 portrait) ===
    print("\n==== Building collage PDFs ====")
    # Build full year collage
    build_collages_portrait("./PDF_species", "All", mode="vector")
    # Build 6-10 month collage
    build_collages_portrait("./PDF_species", "6-10", mode="vector")

if __name__ == '__main__':
    main()