import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from PIL import Image

# ------------------------------
# LOAD DATA
# ------------------------------
@st.cache_data
def load_data():
    df = pd.read_csv(
        "genome_introgression_segments_merged_bins.tsv",
        sep="\t",
        comment="#"
    )
    return df

segments_df = load_data()

# Precompute dynamic NSS threshold base
length_threshold = segments_df['Length'] / 100000


# ------------------------------
# SIDEBAR CONTROLS
# ------------------------------
st.sidebar.title("Desert Filter Controls")

intro = st.sidebar.selectbox(
    "Introgression class",
    ["low", "medium", "high"],
    index=0
)

avail = st.sidebar.slider(
    "Minimum available_bp",
    min_value=0.0,
    max_value=1.0,
    value=0.8,
    step=0.01
)

min_len = st.sidebar.slider(
    "Minimum segment length (bp)",
    min_value=0,
    max_value=2_000_000,
    value=10_000,
    step=1000
)

nss_per_100kb = st.sidebar.slider(
    "Bins with NSS per 100kb (min = 1 per desert)",
    min_value=0.1,
    max_value=5.0,
    value=1.0,
    step=0.1
)

# ------------------------------
# ADD README EXPLANATION HERE
# ------------------------------

st.sidebar.markdown("---")   # horizontal rule for visual separation

with st.sidebar.expander("What do these filters mean?"):
    st.markdown("""
    **Introgression class**  
    Controls which type of segment is considered:  
    low, medium, or high average Neanderthal introgressed allele frequency.

    **Minimum available_bp**  
    Ensures at least this fraction of the region is *mappable*.  
    Example: `0.8` → ≥ 80% of the bases are real sequence (not masked).

    **Minimum segment length**  
    Removes fragments smaller than selected length.

    **Bin with NSS / 100kb**  
    Filters deserts by 10kb bins with Neanderthal-specific SNPs.  
    Higher threshold → stricter desert definition.
    """)

# ------------------------------
# MAIN COMPUTATION
# ------------------------------
threshold = np.maximum(length_threshold * nss_per_100kb, 1)

results = {}
all_desert_lengths = []

for chrom in range(1, 23):
    label = f"chr{chrom}"
    seg_chr = segments_df[segments_df["CHROM"] == label]

    if seg_chr.empty:
        results[label] = np.nan
        continue

    chrom_len = (seg_chr["END_POS"] - seg_chr["START_POS"] + 1).sum()

    deserts_chr = seg_chr[
        (seg_chr['introgression'] == intro) &
        (seg_chr['available_bp'] >= avail) &
        (seg_chr['Length'] >= min_len) &
        (seg_chr['bin_low_NSS'] > threshold.loc[seg_chr.index])
    ]

    all_desert_lengths.extend(deserts_chr["Length"].tolist())

    desert_len = deserts_chr["Length"].sum()
    results[label] = 100 * desert_len / chrom_len

# Make a table
out_df = pd.DataFrame({
    "CHROM": list(results.keys()),
    "Percent": list(results.values())
})

# ------------------------------
# Filtered Deserts
# ------------------------------
filtered_deserts = []   # NEW: store all deserts across chromosomes

for chrom in range(1, 23):
    label = f"chr{chrom}"
    seg_chr = segments_df[segments_df["CHROM"] == label]

    if seg_chr.empty:
        results[label] = np.nan
        continue

    chrom_len = (seg_chr["END_POS"] - seg_chr["START_POS"] + 1).sum()

    deserts_chr = seg_chr[
        (seg_chr['introgression'] == intro) &
        (seg_chr['available_bp'] >= avail) &
        (seg_chr['Length'] >= min_len) &
        (seg_chr['bin_low_NSS'] > threshold.loc[seg_chr.index])
    ]

    # add filtered segments to master list
    if not deserts_chr.empty:
        filtered_deserts.append(deserts_chr)

    all_desert_lengths.extend(deserts_chr["Length"].tolist())

    desert_len = deserts_chr["Length"].sum()
    results[label] = 100 * desert_len / chrom_len

# after loop:
if filtered_deserts:
    filtered_df = pd.concat(filtered_deserts, ignore_index=True)
else:
    filtered_df = pd.DataFrame()


# ------------------------------
# DISPLAY
# ------------------------------

img = Image.open("neanderthal.webp")
st.image(img, caption="Reconstructed Neanderthal", use_column_width=True)

st.title("Neanderthal Desert Explorer")

st.markdown(
    f"""
### Current filter:
- introgression: **{intro}**
- available_bp ≥ **{avail}**
- length ≥ **{min_len:,} bp**
- NSS threshold: **{nss_per_100kb:.2f} per 100kb (min 1)**  
    """
)

# PLOTS

## Bar plot
st.subheader("Bar plot: desert coverage per chromosome")

orange = "#FFA500"   # clean amber-orange

fig1, ax1 = plt.subplots(figsize=(10,4))
ax1.bar(out_df["CHROM"], out_df["Percent"], color=orange)

ax1.set_ylabel("% desert")
ax1.set_xlabel("Chromosome")
ax1.set_title("Desert coverage per chromosome")
plt.xticks(rotation=90)
st.pyplot(fig1)



## Histogram
st.subheader("Histogram: length distribution of filtered deserts")

fig2, ax2 = plt.subplots(figsize=(10,4))

counts, bins = np.histogram(all_desert_lengths, bins=30)
counts = counts.astype(float)

# normalize: 0 -> min count, 1 -> max count
norm = 1 - counts / counts.max()   # small -> light, big -> dark with YlOrBr
cmap = cm.get_cmap("YlGnBu")   # light to dark orange

for i in range(len(counts)):
    color = cmap(norm[i])       # bigger count -> higher norm -> darker color
    ax2.bar(
        bins[i],
        counts[i],
        width=bins[i+1] - bins[i],
        align="edge",
        color=color,
        edgecolor="black",
        linewidth=0.3,
    )

ax2.set_xlabel("Length (bp)")
ax2.set_ylabel("Count")
ax2.set_title("Distribution of desert lengths (filtered)")
plt.ticklabel_format(style="plain", axis="x")

st.pyplot(fig2)

# Table
st.subheader("Desert coverage per chromosome")
st.dataframe(out_df.style.format({"Percent": "{:.2f}"}))

st.subheader("Filtered Desert Segments")

if filtered_df.empty:
    st.info("No desert segments matched your filters.")
else:
    # scrollable table
    st.dataframe(filtered_df, height=400)

    # downloadable CSV
    csv_data = filtered_df.to_csv(index=False).encode("utf-8")
    st.download_button(
        label="⬇ Download filtered deserts as CSV",
        data=csv_data,
        file_name="filtered_deserts.csv",
        mime="text/csv"
    )


