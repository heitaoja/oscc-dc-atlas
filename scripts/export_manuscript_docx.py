# -*- coding: utf-8 -*-
"""Export publication-grade DOCX and LaTeX submission files from manuscript.md and live result tables."""

from __future__ import annotations

import csv
import re
from pathlib import Path

from docx import Document
from docx.enum.table import WD_ALIGN_VERTICAL, WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement, parse_xml
from docx.oxml.ns import nsdecls, qn
from docx.shared import Cm, Inches, Pt, RGBColor

ROOT = Path(__file__).resolve().parent.parent
MANUSCRIPT_MD = ROOT / "report" / "manuscript.md"
TAB_DIR = ROOT / "results" / "tables"
FIG_DIR = ROOT / "results" / "figures"
OUT_DOCX = ROOT / "report" / "manuscript_submission.docx"
OUT_TEX = ROOT / "report" / "manuscript_submission.tex"

REPO_URL = "https://github.com/heitaoja/oscc-dc-atlas"
RUNNING_TITLE = "Migratory DC transcriptional programs in OSCC"

GENE_SYMBOLS = sorted(
    [
        "LAMP3", "CCR7", "CCL19", "CCL22", "CCL17", "CCR4",
        "CD274", "PDCD1LG2", "CD200", "IDO1", "CD83", "CD1C",
        "CLEC10A", "FCER1A", "CLEC9A", "XCR1", "CADM1", "FSCN1",
        "CD207", "CD1A", "CD1E", "CD163", "BIRC3", "RELB",
        "NFKBIA", "TNFAIP3", "TRAF1", "TNIP1", "TREM1", "OLR1",
        "CD300E", "THBS1", "F13A1", "SOCS2", "PDCD1", "CTLA4",
        "TIGIT", "LAG3", "HAVCR2", "FOXP3",
    ],
    key=len,
    reverse=True,
)
GENE_SET = set(GENE_SYMBOLS)

REFERENCES = [
    (
        1,
        "Maier B, Leader AM, Chen ST, et al. A conserved dendritic-cell regulatory program limits antitumour immunity.",
        "Nature",
        "2020;580(7802):257–262. doi:10.1038/s41586-020-2134-y",
    ),
    (
        2,
        "Solhaug MB, Schreurs O, Schenck K, et al. Origin of langerin (CD207)-expressing antigen presenting cells in the normal oral mucosa and in oral lichen planus lesions.",
        "Eur J Oral Sci",
        "2022;130(1):e12835. doi:10.1111/eos.12835",
    ),
    (
        3,
        "Werner SL, Kearns JD, Zadorozhnaya V, et al. Encoding NF-kappaB temporal control in response to TNF: distinct roles for the negative regulators IkappaBalpha and A20.",
        "Genes Dev",
        "2008;22(15):2093–2101. doi:10.1101/gad.1680708",
    ),
    (
        4,
        "Squair JW, Gautier M, Kathe C, et al. Confronting false discoveries in single-cell differential expression.",
        "Nat Commun",
        "2021;12(1):5692. doi:10.1038/s41467-021-25960-2",
    ),
    (
        5,
        "Peng Y, Xiao L, Rong H, et al. Single-cell profiling of tumor-infiltrating TCF1/TCF7(+) T cells reveals a T lymphocyte subset associated with tertiary lymphoid structures/organs and a superior prognosis in oral cancer.",
        "Oral Oncol",
        "2021;119:105348. doi:10.1016/j.oraloncology.2021.105348",
    ),
    (
        6,
        "Chen X, Cai Y, Hu X, et al. Differential metabolic requirement governed by transcription factor c-Maf dictates innate γδT17 effector functionality in mice and humans.",
        "Sci Adv",
        "2022;8(21):eabm9120. doi:10.1126/sciadv.abm9120",
    ),
    (
        7,
        "Robinson MD, McCarthy DJ, Smyth GK. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data.",
        "Bioinformatics",
        "2010;26(1):139–140. doi:10.1093/bioinformatics/btp616",
    ),
    (
        8,
        "Ritchie ME, Phipson B, Wu D, et al. limma powers differential expression analyses for RNA-sequencing and microarray studies.",
        "Nucleic Acids Res",
        "2015;43(7):e47. doi:10.1093/nar/gkv007",
    ),
    (
        9,
        "Wu D, Smyth GK. Camera: a competitive gene set test accounting for inter-gene correlation.",
        "Nucleic Acids Res",
        "2012;40(17):e133. doi:10.1093/nar/gks461",
    ),
    (
        10,
        "Choi JH, Lee BS, Jang JY, et al. Single-cell transcriptome profiling of the stepwise progression of head and neck cancer.",
        "Nat Commun",
        "2023;14(1):1055. doi:10.1038/s41467-023-36691-x",
    ),
    (
        11,
        "Heyninck K, De Valck D, Vanden Berghe W, et al. The zinc finger protein A20 inhibits TNF-induced NF-kappaB-dependent gene expression by interfering with an RIP- or TRAF2-mediated transactivation signal and directly binds to a novel NF-kappaB-inhibiting protein ABIN.",
        "J Cell Biol",
        "1999;145(7):1471–1482. doi:10.1083/jcb.145.7.1471",
    ),
    (
        12,
        "Iellem A, Mariani M, Lang R, et al. Unique chemotactic response profile and specific expression of chemokine receptors CCR4 and CCR8 by CD4(+)CD25(+) regulatory T cells.",
        "J Exp Med",
        "2001;194(6):847–853. doi:10.1084/jem.194.6.847",
    ),
]

SUPPLEMENTARY_INDEX = [
    (
        "Supplementary Table S1",
        "Discovery atlas library composition, matched tissue fractions, and cell-type counts.",
        "results/tables/03_celltype_by_sample.csv; 09_dc_fraction_by_sample.csv; 09_dc_subtype_pooled.csv; 09_paired_tissue_descriptive.csv",
    ),
    (
        "Supplementary Table S2",
        "Discovery donor-paired differential expression across DC lineages and lineage marker panel estimates.",
        "results/tables/11_de_migDC_vs_cDC2.csv; 11_de_migDC_vs_cDC1.csv; 11_de_cDC1_vs_cDC2.csv; 11_de_migDC_vs_LClike.csv; 11_de_LClike_vs_cDC2.csv; 11_panel_model_results.csv; 11_primary_pseudobulk_counts.tsv.gz; 11_primary_pseudobulk_libraries.csv",
    ),
    (
        "Supplementary Table S3",
        "Discovery sensitivity analyses, leave-one-donor-out checks, equal-cell subsampling, and objective marker-gate reclassification.",
        "results/tables/11_sensitivity.csv; 11_leave_one_donor_out.csv; 11_balanced_cell_sensitivity.csv; 11_marker_gate_sensitivity.csv; config/external_dc_gates.csv",
    ),
    (
        "Supplementary Table S4",
        "Held-out expression program definitions, donor-level effects, limma::fry tests, and gate-based regulatory program results.",
        "config/dc_programs.csv; results/tables/11_program_donor_effects.csv; 11_program_summary.csv; 11_program_fry.csv; 11_marker_gate_programs.csv",
    ),
    (
        "Supplementary Table S5",
        "Competitive Hallmark pathway enrichment (CAMERA) on assignment-excluded discovery genes.",
        "results/tables/11_paired_hallmark_camera.csv",
    ),
    (
        "Supplementary Table S6",
        "Patient-level tumour microenvironment (TME) immune proportions, pseudobulk profiles, and Spearman covariation tests.",
        "results/tables/13_tme_cell_proportions.csv; 13_tme_donor_profiles.csv; 13_tme_ligand_receptor_correlations.csv",
    ),
    (
        "Supplementary Table S7",
        "GSE181919 external oral-cavity cohort inventory, verified GEO file ledger, and myeloid APC gate classification audit.",
        "results/tables/12_inputs.csv; 12_cohort_inventory.csv; 12_classification.csv",
    ),
    (
        "Supplementary Table S8",
        "GSE181919 external paired differential expression sensitivities and cross-cohort gene-level replication.",
        "results/tables/12_sensitivity.csv; 12_cross_cohort_summary.csv; 12_cross_cohort_genes.csv",
    ),
    (
        "Supplementary Table S9",
        "GSE181919 external held-out and expanded regulatory program summaries and limma::fry results.",
        "results/tables/12_program_summary.csv; 12_program_fry.csv",
    ),
    (
        "Supplementary Table S10",
        "Cross-cohort shared-gene donor program scores and combined Discovery + GSE181919 meta-analysis.",
        "results/tables/12_shared_program_donor_effects.csv; 12_combined_program_summary.csv",
    ),
    (
        "Supplementary Figure S1",
        "Discovery atlas UMAP across 85,164 quality-controlled cells and per-library quality-control distributions.",
        "results/figures/03_umap_celltype.png; 10_qc_cell_numbers.png; 10_qc_mito.png",
    ),
    (
        "Supplementary Figure S2",
        "Discovery DC subclustering UMAP, subtype marker dotplot, and paired-program forest plot.",
        "results/figures/05_dc_subtype_umap.png; 05_dc_subtype_dotplot.png; 11_paired_dc_programs.png",
    ),
    (
        "Supplementary Figure S3",
        "Lineage comparator gene panel across five paired contrasts and complete 46-axis TME covariation overview.",
        "results/figures/11_dc_lineage_comparators.png; 13_tme_covariation.png",
    ),
    (
        "Supplementary Figure S4",
        "Cross-cohort gene-level effect concordance in GSE181919.",
        "results/figures/12_cross_cohort_effects.png",
    ),
]


def read_csv_rows(name: str) -> list[dict[str, str]]:
    with open(TAB_DIR / name, encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def fmt_p(val: float) -> str:
    if val >= 0.001:
        return f"{val:.4g}" if val < 0.01 else f"{val:.3g}"
    s = f"{val:.2e}"
    mant, exp = s.split("e")
    exp_i = int(exp)
    sup_map = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")
    return f"{mant} × 10{str(exp_i).translate(sup_map)}"


def fmt_p_tex(val: float) -> str:
    if val >= 0.001:
        return f"{val:.4g}" if val < 0.01 else f"{val:.3g}"
    s = f"{val:.2e}"
    mant, exp = s.split("e")
    return f"${mant} \\times 10^{{{int(exp)}}}$"


def fmt_signed(val: float, digits: int = 2) -> str:
    return f"+{val:.{digits}f}" if val >= 0 else f"−{abs(val):.{digits}f}"


def strip_md_links_and_normalize(text: str) -> str:
    # Remove trailing link-only clusters that look like: " [Link](...) · [Link2](...)"
    # while replacing literature citations with numbered Vancouver citations.
    citation_replacements = [
        ("[Maier et al.](https://pubmed.ncbi.nlm.nih.gov/32269339/)", "[1]"),
        ("[Human oral mucosal langerin-expressing cells](https://pubmed.ncbi.nlm.nih.gov/34935212/)", "[2]"),
        ("[NF-κB feedback kinetics](https://pubmed.ncbi.nlm.nih.gov/18676814/)", "[3]"),
        ("[Squair et al.](https://doi.org/10.1038/s41467-021-25960-2)", "[4]"),
        ("[edgeR](https://bioconductor.org/packages/edgeR)", "[7]"),
        ("[limma](https://bioconductor.org/packages/limma)", "[8]"),
        ("[Wu and Smyth](https://pubmed.ncbi.nlm.nih.gov/22638577/)", "[9]"),
        ("[Choi et al.](https://pubmed.ncbi.nlm.nih.gov/36828832/)", "[10]"),
        ("[NF-κB temporal feedback](https://pubmed.ncbi.nlm.nih.gov/18676814/) · [A20 and ABIN](https://pubmed.ncbi.nlm.nih.gov/10385526/)", "[3,11]"),
        ("[Oral mucosal APCs](https://pubmed.ncbi.nlm.nih.gov/34935212/)", "[2]"),
        ("[Treg chemotaxis](https://pubmed.ncbi.nlm.nih.gov/11560999/)", "[12]"),
    ]
    for old, new in citation_replacements:
        text = text.replace(old, new)

    # Specific prose insertions for cohort datasets and supplementary references
    text = text.replace("from GSE172577 and", "from GSE172577 [5] and")
    text = text.replace("from GSE193766.", "from GSE193766 [6].")
    text = text.replace("GSE181919 primary oral-cavity tumours", "GSE181919 [10] primary oral-cavity tumours")
    text = text.replace("deposited in GEO.", "deposited in GEO [5,6,10].")
    text = text.replace(
        "The discovery and validation data are publicly deposited in GEO [5,6,10].",
        "The discovery and validation data are publicly deposited in GEO under accessions GSE172577 [5], GSE193766 [6], and GSE181919 [10].",
    )

    # Remove markdown links to repository tables, configs, scripts, PDFs, or GEO URLs
    text = re.sub(r"\s*(?:·\s*)?\[[^\]]+\]\([^)]+\)", "", text)
    # Remove backticks around code identifiers
    text = text.replace("`", "")
    # Clean up extra whitespace before punctuation
    text = re.sub(r"\s+([.,;:])", r"\1", text)
    text = re.sub(r"\s{2,}", " ", text).strip()
    return text


def attach_supplementary_callout(para_text: str, section_key: str) -> str:
    callouts = {
        "abstract_methods": " (Supplementary Table S1)",
        "abstract_results": " (Tables 1–2; Supplementary Tables S1–S10)",
        "2.1_p1": " (Supplementary Table S1).",
        "2.1_p2": " (Supplementary Table S1).",
        "2.2_p1": " (Supplementary Table S3).",
        "2.2_p2": " (Supplementary Tables S2–S3).",
        "2.2_p3": " (Supplementary Table S3).",
        "2.3_p1": " (Supplementary Table S4).",
        "2.3_p2": " (Supplementary Tables S4–S5).",
        "2.4_p1": " (Supplementary Table S6).",
        "2.4_p2": " (Supplementary Table S6).",
        "2.5_p1": " (Supplementary Table S7).",
        "2.5_p2": " (Supplementary Tables S8–S9).",
        "2.5_p3": " (Supplementary Table S10).",
        "3.1_p1": " (Fig. 1A–C; Supplementary Table S1; Supplementary Figs. S1–S2).",
        "3.1_p2": " (Fig. 1D; Supplementary Table S1).",
        "3.2_p1": " (Fig. 2A; Supplementary Table S2).",
        "3.2_p2": " (Fig. 2B–C; Supplementary Tables S2–S3).",
        "3.3_p1": " (Table 1; Fig. 3A; Supplementary Table S4; Supplementary Fig. S2).",
        "3.3_p2": " (Supplementary Table S4).",
        "3.3_p3": " (Fig. 3B–C; Supplementary Tables S2, S5).",
        "3.4_p1": " (Fig. 4A; Supplementary Tables S2–S3; Supplementary Fig. S3).",
        "3.4_p2": " (Fig. 4B–C; Supplementary Table S6; Supplementary Fig. S3).",
        "3.5_p1": " (Fig. 5A; Supplementary Tables S7–S8).",
        "3.5_p2": " (Table 2; Fig. 5B; Supplementary Table S8; Supplementary Fig. S4).",
        "3.5_p3": " (Supplementary Table S9).",
        "3.5_p4": " (Table 2; Fig. 5C; Supplementary Table S10).",
        "4_p1": " (Supplementary Tables S4–S5).",
        "4_p2": " (Supplementary Tables S2, S4).",
        "4_p3": " (Supplementary Tables S6–S8).",
        "4_p4": " (Supplementary Tables S1, S3).",
        "avail_p2": f" ({REPO_URL}; Supplementary Tables S1–S10).",
    }
    suffix = callouts.get(section_key)
    if not suffix:
        return para_text
    if suffix.startswith(" (") and para_text.endswith("."):
        return para_text[:-1] + suffix
    return para_text + suffix


def parse_manuscript() -> dict:
    raw = MANUSCRIPT_MD.read_text(encoding="utf-8")
    lines = raw.splitlines()
    title = lines[0].lstrip("# ").strip()

    sections: dict[str, list[str]] = {}
    current = None
    in_code = False
    for line in lines[1:]:
        if line.strip().startswith("~~~"):
            in_code = not in_code
            continue
        if in_code:
            continue
        if line.startswith("## ") or line.startswith("### "):
            current = line.lstrip("# ").strip()
            sections[current] = []
            continue
        if current is None:
            continue
        s = line.strip()
        if not s or s.startswith("![") or s.startswith("|") or s.startswith("**Figure "):
            continue
        sections[current].append(s)

    # Extract figure legends explicitly
    legends = {}
    for line in lines:
        m = re.match(r"^\*\*Figure (\d+)\.\s+([^*]+)\*\*\s*(.+)$", line.strip())
        if m:
            num = int(m.group(1))
            fig_title = m.group(2).strip()
            fig_body = strip_md_links_and_normalize(m.group(3).strip())
            legends[num] = (fig_title, fig_body)

    return {"title": title, "sections": sections, "legends": legends}


def tokenize_rich_text(text: str):
    """Parse **bold**, *italic*, genes, statistical symbols, subscripts, superscripts."""
    # First split by markdown **bold** and *italic*
    pattern = re.compile(r"(\*\*[^*]+\*\*|\*[^*]+\*)")
    tokens = []
    for part in pattern.split(text):
        if not part:
            continue
        if part.startswith("**") and part.endswith("**"):
            tokens.extend(_sub_tokenize(part[2:-2], bold=True, italic=False))
        elif part.startswith("*") and part.endswith("*"):
            tokens.extend(_sub_tokenize(part[1:-1], bold=False, italic=True))
        else:
            tokens.extend(_sub_tokenize(part, bold=False, italic=False))
    return tokens


STAT_SUP_MAP = {
    "⁰": "0", "¹": "1", "²": "2", "³": "3", "⁴": "4",
    "⁵": "5", "⁶": "6", "⁷": "7", "⁸": "8", "⁹": "9", "⁻": "-",
}
SUP_SEQ_RE = re.compile(r"[⁰¹²³⁴⁵⁶⁷⁸⁹⁻]+")
SUB_SEQ_RE = re.compile(r"log₂")
GENE_RE = re.compile(r"\b(" + "|".join(re.escape(g) for g in GENE_SYMBOLS) + r")\b")
STAT_T_RE = re.compile(r"(?<![A-Za-z0-9-])t(?=\s+(?:interval|intervals|test|tests|CI|P\b|value|values))")
STAT_P_RE = re.compile(r"(?<![A-Za-z0-9-])P(?=\s*(?:[=<>]|value|values|\())")
STAT_Q_RE = re.compile(r"(?<![A-Za-z0-9-])q(?=\s*[=<>])")
ET_AL_RE = re.compile(r"\bet al\.")


def _sub_tokenize(text: str, bold: bool, italic: bool):
    # Split sequentially for log₂, superscripts, genes, and statistical t/P/q
    combined = re.compile(
        r"(log₂|[⁰¹²³⁴⁵⁶⁷⁸⁹⁻]+|\b(?:"
        + "|".join(re.escape(g) for g in GENE_SYMBOLS)
        + r")\b|(?<![A-Za-z0-9-])t(?=\s+(?:interval|intervals|test|tests|CI|P\b|value|values))|(?<![A-Za-z0-9-])P(?=\s*(?:[=<>]|value|values|\())|(?<![A-Za-z0-9-])q(?=\s*[=<>])|\bet al\.)"
    )
    out = []
    pos = 0
    for m in combined.finditer(text):
        if m.start() > pos:
            out.append((text[pos:m.start()], bold, italic, False, False))
        val = m.group(0)
        if val == "log₂":
            out.append(("log", bold, italic, False, False))
            out.append(("2", bold, italic, False, True))
        elif SUP_SEQ_RE.fullmatch(val):
            plain = "".join(STAT_SUP_MAP[ch] for ch in val)
            out.append((plain, bold, italic, True, False))
        elif val in GENE_SET or val == "et al." or val in {"t", "P", "q"}:
            out.append((val, bold, True, False, False))
        else:
            out.append((val, bold, italic, False, False))
        pos = m.end()
    if pos < len(text):
        out.append((text[pos:], bold, italic, False, False))
    return out


def add_rich_runs(paragraph, text: str, font_size: float = 11.0, force_bold: bool = False):
    for chunk, bold, italic, sup, sub in tokenize_rich_text(text):
        run = paragraph.add_run(chunk)
        run.font.name = "Times New Roman"
        run._element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
        run.font.size = Pt(font_size)
        run.font.color.rgb = RGBColor(0x1A, 0x1A, 0x1A)
        run.bold = bold or force_bold
        run.italic = italic
        run.font.superscript = sup
        run.font.subscript = sub


def add_page_number_to_run(run):
    fld_char1 = parse_xml(r'<w:fldChar %s w:fldCharType="begin"/>' % nsdecls("w"))
    instr_text = parse_xml(r'<w:instrText %s xml:space="preserve"> PAGE </w:instrText>' % nsdecls("w"))
    fld_char2 = parse_xml(r'<w:fldChar %s w:fldCharType="separate"/>' % nsdecls("w"))
    fld_char3 = parse_xml(r'<w:fldChar %s w:fldCharType="end"/>' % nsdecls("w"))
    run._r.append(fld_char1)
    run._r.append(instr_text)
    run._r.append(fld_char2)
    run._r.append(fld_char3)


def configure_three_line_table(table, col_widths_cm: list[float], right_align_cols: set[int]):
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    table.autofit = False

    # Clear all borders at table level first
    tbl_pr = table._tbl.tblPr
    borders = parse_xml(
        f'<w:tblBorders {nsdecls("w")}>'
        '<w:top w:val="single" w:sz="12" w:space="0" w:color="111111"/>'
        '<w:left w:val="none"/>'
        '<w:bottom w:val="single" w:sz="12" w:space="0" w:color="111111"/>'
        '<w:right w:val="none"/>'
        '<w:insideH w:val="none"/>'
        '<w:insideV w:val="none"/>'
        "</w:tblBorders>"
    )
    existing = tbl_pr.find(qn("w:tblBorders"))
    if existing is not None:
        tbl_pr.remove(existing)
    tbl_pr.append(borders)

    n_rows = len(table.rows)
    for r_idx, row in enumerate(table.rows):
        tr_pr = row._tr.get_or_add_trPr()
        tr_pr.append(parse_xml(f'<w:cantSplit {nsdecls("w")}/>'))
        if r_idx == 0:
            tr_pr.append(parse_xml(f'<w:tblHeader {nsdecls("w")}/>'))
        for c_idx, cell in enumerate(row.cells):
            cell.width = Cm(col_widths_cm[c_idx])
            cell.vertical_alignment = WD_ALIGN_VERTICAL.CENTER
            tc_pr = cell._tc.get_or_add_tcPr()
            # Cell padding
            mar = parse_xml(
                f'<w:tcMar {nsdecls("w")}>'
                '<w:top w:w="70" w:type="dxa"/>'
                '<w:bottom w:w="70" w:type="dxa"/>'
                '<w:left w:w="80" w:type="dxa"/>'
                '<w:right w:w="80" w:type="dxa"/>'
                "</w:tcMar>"
            )
            tc_pr.append(mar)
            # Header bottom rule (0.75pt = sz 6) and last row bottom rule (1.5pt = sz 12)
            b_xml = []
            if r_idx == 0:
                b_xml.append('<w:top w:val="single" w:sz="12" w:space="0" w:color="111111"/>')
                b_xml.append('<w:bottom w:val="single" w:sz="6" w:space="0" w:color="111111"/>')
            elif r_idx == n_rows - 1:
                b_xml.append('<w:bottom w:val="single" w:sz="12" w:space="0" w:color="111111"/>')
            if b_xml:
                tc_borders = parse_xml(f'<w:tcBorders {nsdecls("w")}>{"".join(b_xml)}</w:tcBorders>')
                tc_pr.append(tc_borders)
            for p in cell.paragraphs:
                p.paragraph_format.space_before = Pt(1)
                p.paragraph_format.space_after = Pt(1)
                p.paragraph_format.line_spacing = 1.15
                p.alignment = WD_ALIGN_PARAGRAPH.RIGHT if c_idx in right_align_cols else WD_ALIGN_PARAGRAPH.LEFT


def build_live_tables():
    prog_rows = [r for r in read_csv_rows("11_program_summary.csv") if r["analysis"] == "primary"]
    order = ["Maturation", "NF-kB feedback", "Regulatory expression", "MHC class II"]
    prog_map = {r["program"]: r for r in prog_rows}
    t1_rows = []
    for pname in order:
        r = prog_map[pname]
        disp_name = pname.replace("NF-kB", "NF-κB")
        cov = f"{r['n_genes']}/{r['n_target_genes']}"
        mean_eff = fmt_signed(float(r["mean_effect"]), 2)
        ci = f"({fmt_signed(float(r['lower_95_t']), 2)}, {fmt_signed(float(r['upper_95_t']), 2)})"
        donors = f"{r['positive_donors']}/{r['n_donors']}"
        p_flip = fmt_p(float(r["p_sign_flip"]))
        p_t = fmt_p(float(r["p_paired_t"]))
        fdr_t = fmt_p(float(r["fdr_paired_t"]))
        fdr_fry = fmt_p(float(r["fdr_fry"]))
        t1_rows.append((disp_name, cov, f"{mean_eff} {ci}", donors, p_flip, p_t, fdr_t, fdr_fry, r))

    ext_rows = [r for r in read_csv_rows("12_cross_cohort_summary.csv") if r["status"] == "fitted"]
    comb_rows = read_csv_rows("12_combined_program_summary.csv")
    return t1_rows, ext_rows, comb_rows


def build_docx(data: dict, t1_rows, ext_rows, comb_rows):
    doc = Document()
    sec = doc.sections[0]
    sec.page_width = Cm(21.0)
    sec.page_height = Cm(29.7)
    sec.top_margin = Cm(2.54)
    sec.bottom_margin = Cm(2.54)
    sec.left_margin = Cm(2.54)
    sec.right_margin = Cm(2.54)
    sec.different_first_page_header_footer = True

    normal = doc.styles["Normal"]
    normal.font.name = "Times New Roman"
    normal._element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
    normal.font.size = Pt(11)
    normal.font.color.rgb = RGBColor(0x1A, 0x1A, 0x1A)

    # Running header and footer (pages 2+)
    hdr_p = sec.header.paragraphs[0]
    hdr_p.alignment = WD_ALIGN_PARAGRAPH.RIGHT
    hdr_r = hdr_p.add_run(RUNNING_TITLE)
    hdr_r.font.name = "Times New Roman"
    hdr_r.font.size = Pt(9)
    hdr_r.font.italic = True
    hdr_r.font.color.rgb = RGBColor(0x55, 0x55, 0x55)

    ftr_p = sec.footer.paragraphs[0]
    ftr_p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    ftr_r1 = ftr_p.add_run("Page ")
    ftr_r1.font.name = "Times New Roman"
    ftr_r1.font.size = Pt(9)
    ftr_r2 = ftr_p.add_run()
    ftr_r2.font.name = "Times New Roman"
    ftr_r2.font.size = Pt(9)
    add_page_number_to_run(ftr_r2)

    def add_body_p(text: str, space_after: float = 6.0, first_line_indent: float = 0.0, bold: bool = False):
        p = doc.add_paragraph()
        p.paragraph_format.line_spacing = 1.5
        p.paragraph_format.space_before = Pt(0)
        p.paragraph_format.space_after = Pt(space_after)
        if first_line_indent:
            p.paragraph_format.first_line_indent = Cm(first_line_indent)
        add_rich_runs(p, text, font_size=11.0, force_bold=bold)
        return p

    def add_heading_p(text: str, level: int = 1):
        p = doc.add_paragraph()
        p.paragraph_format.keep_with_next = True
        if level == 1:
            p.paragraph_format.space_before = Pt(14)
            p.paragraph_format.space_after = Pt(6)
            add_rich_runs(p, text, font_size=13.0, force_bold=True)
        else:
            p.paragraph_format.space_before = Pt(10)
            p.paragraph_format.space_after = Pt(4)
            add_rich_runs(p, text, font_size=11.5, force_bold=True)
        return p

    sections = data["sections"]

    # Prepare cleaned text blocks and word counts
    abs_paras = [strip_md_links_and_normalize(x) for x in sections["Abstract"]]
    abs_paras[1] = attach_supplementary_callout(abs_paras[1], "abstract_methods")
    abs_paras[3] = attach_supplementary_callout(abs_paras[3], "abstract_results")
    # Drop **Background** etc. header lines from abstract word count
    abstract_words = sum(
        len(re.sub(r"\*\*[^*]+\*\*", "", p).split())
        for p in abs_paras
        if not p.startswith("**Keywords**") and "Oral squamous cell carcinoma, dendritic cells" not in p
    )

    main_keys = [
        ("1. Introduction", None, ["1.1", "1.2", "1.3"]),
        ("2.1. Discovery cohorts, annotation and composition", "2.1", None),
        ("2.2. Donor-paired pseudobulk models and sensitivity analyses", "2.2", None),
        ("2.3. Held-out expression programs and competitive pathways", "2.3", None),
        ("2.4. Patient-level tumour microenvironment covariation", "2.4", None),
        ("2.5. Independent external validation and cross-cohort pooling", "2.5", None),
        ("3.1. Discovery atlas, DC heterogeneity, and matched-tissue composition", "3.1", None),
        ("3.2. Donor-paired differential expression and methodological robustness", "3.2", None),
        ("3.3. Held-out transcriptional programs and resolution of the inflammatory paradox", "3.3", None),
        ("3.4. DC lineage specificity including mucosal LC-like states and TME covariation", "3.4", None),
        ("3.5. Independent validation in GSE181919 and cross-cohort meta-analysis", "3.5", None),
        ("4. Discussion", "4", None),
    ]
    cleaned_sections: dict[str, list[str]] = {}
    main_words = 0
    for sec_name, tag, _ in main_keys:
        plist = []
        for idx, raw_p in enumerate(sections[sec_name], start=1):
            cp = strip_md_links_and_normalize(raw_p)
            if tag:
                cp = attach_supplementary_callout(cp, f"{tag}_p{idx}")
            plist.append(cp)
            main_words += len(cp.split())
        cleaned_sections[sec_name] = plist

    # ---- Title Page ----
    tp_title = doc.add_paragraph()
    tp_title.paragraph_format.space_before = Pt(24)
    tp_title.paragraph_format.space_after = Pt(14)
    tp_title.paragraph_format.line_spacing = 1.25
    add_rich_runs(tp_title, data["title"], font_size=16.0, force_bold=True)

    add_body_p(f"**Running title:** {RUNNING_TITLE}", space_after=8)
    add_body_p("**Authors:** [Author Name(s)]", space_after=4)
    add_body_p("**Affiliation:** [Department, Institution, City, Country]", space_after=4)
    add_body_p("**Corresponding author:** [Name], [Email]; ORCID: [0000-0000-0000-0000]", space_after=8)
    add_body_p(
        "**Keywords:** Oral squamous cell carcinoma; dendritic cells; pseudobulk; "
        "Langerhans-like cells; tumour microenvironment; external validation",
        space_after=8,
    )
    add_body_p(
        f"**Word counts:** Abstract: {abstract_words} words; Main text (Introduction through Discussion): "
        f"{main_words} words; Main tables: 2; Main figures: 5; References: {len(REFERENCES)}.",
        space_after=8,
    )
    add_body_p(f"**Code and reproducible archive:** {REPO_URL}", space_after=12)

    doc.add_page_break()

    # ---- Abstract ----
    add_heading_p("Abstract", level=1)
    abs_labels = ["Background", "Methods", "Results", "Conclusions"]
    abs_items = {}
    cur_lbl = None
    for item in sections["Abstract"]:
        clean_item = strip_md_links_and_normalize(item)
        if clean_item.startswith("**") and clean_item.endswith("**"):
            cur_lbl = clean_item.strip("*")
        elif cur_lbl:
            abs_items[cur_lbl] = clean_item
    abs_items["Methods"] = attach_supplementary_callout(abs_items["Methods"], "abstract_methods")
    abs_items["Results"] = attach_supplementary_callout(abs_items["Results"], "abstract_results")

    for lbl in abs_labels:
        add_body_p(f"**{lbl}:** {abs_items[lbl]}", space_after=6)
    add_body_p(f"**Keywords:** {abs_items['Keywords']}", space_after=10)

    # Recompute exact abstract word count on the 4 body sections
    abstract_words_exact = sum(len(abs_items[lbl].split()) for lbl in abs_labels)
    for p in doc.paragraphs:
        if p.text.startswith("Word counts:"):
            p.clear()
            add_rich_runs(
                p,
                f"**Word counts:** Abstract: {abstract_words_exact} words; Main text (Introduction through Discussion): "
                f"{main_words} words; Main tables: 2; Main figures: 5; References: {len(REFERENCES)}.",
                font_size=11.0,
            )

    # ---- 1. Introduction ----
    add_heading_p("1. Introduction", level=1)
    for cp in cleaned_sections["1. Introduction"]:
        add_body_p(cp)

    # ---- 2. Materials and methods ----
    add_heading_p("2. Materials and methods", level=1)
    for sec_name in [
        "2.1. Discovery cohorts, annotation and composition",
        "2.2. Donor-paired pseudobulk models and sensitivity analyses",
        "2.3. Held-out expression programs and competitive pathways",
        "2.4. Patient-level tumour microenvironment covariation",
        "2.5. Independent external validation and cross-cohort pooling",
    ]:
        add_heading_p(sec_name, level=2)
        for cp in cleaned_sections[sec_name]:
            add_body_p(cp)

    # ---- 3. Results ----
    add_heading_p("3. Results", level=1)

    add_heading_p("3.1. Discovery atlas, DC heterogeneity, and matched-tissue composition", level=2)
    for cp in cleaned_sections["3.1. Discovery atlas, DC heterogeneity, and matched-tissue composition"]:
        add_body_p(cp)

    add_heading_p("3.2. Donor-paired differential expression and methodological robustness", level=2)
    for cp in cleaned_sections["3.2. Donor-paired differential expression and methodological robustness"]:
        add_body_p(cp)

    add_heading_p("3.3. Held-out transcriptional programs and resolution of the inflammatory paradox", level=2)
    for cp in cleaned_sections["3.3. Held-out transcriptional programs and resolution of the inflammatory paradox"]:
        add_body_p(cp)

    # ---- Table 1 (inline in 3.3) ----
    t1_cap = doc.add_paragraph()
    t1_cap.paragraph_format.space_before = Pt(8)
    t1_cap.paragraph_format.space_after = Pt(4)
    t1_cap.paragraph_format.keep_with_next = True
    add_rich_runs(
        t1_cap,
        "**Table 1.** Held-out transcriptional expression programs in the discovery primary paired comparison "
        "(migDC versus cDC2; n = 5 matched tumour donors).",
        font_size=10.5,
    )
    t1_headers = [
        "Held-out program",
        "Genes",
        "Mean Δ log₂ CPM (95% t CI)",
        "Positive donors",
        "Sign-flip P",
        "Paired t P",
        "Paired t FDR",
        "Directional fry FDR",
    ]
    tbl1 = doc.add_table(rows=1 + len(t1_rows), cols=len(t1_headers))
    for j, h in enumerate(t1_headers):
        add_rich_runs(tbl1.rows[0].cells[j].paragraphs[0], h, font_size=9.0, force_bold=True)
    for i, row_data in enumerate(t1_rows, start=1):
        vals = row_data[:8]
        for j, v in enumerate(vals):
            add_rich_runs(tbl1.rows[i].cells[j].paragraphs[0], v, font_size=9.0)
    configure_three_line_table(tbl1, [3.3, 1.2, 3.6, 1.5, 1.5, 1.6, 1.6, 1.7], {1, 2, 3, 4, 5, 6, 7})

    t1_note = doc.add_paragraph()
    t1_note.paragraph_format.space_before = Pt(3)
    t1_note.paragraph_format.space_after = Pt(8)
    add_rich_runs(
        t1_note,
        "Note: Explicit annotation and gate markers are excluded from all four panels. "
        "Paired differences represent migDC minus cDC2 on TMM log₂ CPM (prior count 1).",
        font_size=9.0,
    )

    add_heading_p("3.4. DC lineage specificity including mucosal LC-like states and TME covariation", level=2)
    for cp in cleaned_sections["3.4. DC lineage specificity including mucosal LC-like states and TME covariation"]:
        add_body_p(cp)

    add_heading_p("3.5. Independent validation in GSE181919 and cross-cohort meta-analysis", level=2)
    for cp in cleaned_sections["3.5. Independent validation in GSE181919 and cross-cohort meta-analysis"]:
        add_body_p(cp)

    # ---- Table 2 (2A + 2B in 3.5) ----
    t2a_cap = doc.add_paragraph()
    t2a_cap.paragraph_format.space_before = Pt(8)
    t2a_cap.paragraph_format.space_after = Pt(4)
    t2a_cap.paragraph_format.keep_with_next = True
    add_rich_runs(
        t2a_cap,
        "**Table 2.** Independent validation in GSE181919 and combined cross-cohort program meta-analysis. "
        "**A**, Gene-level replication across fitted external comparisons (excluding assignment and gate genes). "
        "**B**, Shared-gene combined program estimates pooling discovery (n = 5) and GSE181919 patients.",
        font_size=10.5,
    )

    ext_labels = {
        "primary": "Myeloid APC pool (min 10 cells)",
        "min_5": "Myeloid APC pool (min 5 cells)",
        "cd163_neg_min10": "CD163-negative myeloid (min 10)",
        "cd163_neg_min5": "CD163-negative myeloid (min 5)",
        "hybrid_author_migDC_vs_cd163neg_cDC2": "Hybrid author migDC vs CD163-neg cDC2",
        "author_dc_migDC_vs_cDC1_min3": "Author DC migDC vs cDC1 (min 3)",
        "migDC_vs_cDC1": "Myeloid APC migDC vs cDC1 (min 5)",
    }
    t2a_headers = [
        "External model (contrast)",
        "Donors",
        "Shared genes",
        "Spearman ρ",
        "Concordant discovery hits",
        "External FDR < 0.05 hits",
    ]
    tbl2a = doc.add_table(rows=1 + len(ext_rows), cols=len(t2a_headers))
    for j, h in enumerate(t2a_headers):
        add_rich_runs(tbl2a.rows[0].cells[j].paragraphs[0], h, font_size=9.0, force_bold=True)
    for i, r in enumerate(ext_rows, start=1):
        lbl = ext_labels.get(r["analysis"], r["analysis"])
        hits = f"{r['concordant_hits']}/{r['n_discovery_hits']} ({100*float(r['concordant_hit_rate']):.1f}%)"
        vals = [lbl, r["n_donors"], f"{int(r['n_shared']):,}", f"{float(r['rho']):.3f}", hits, r["external_fdr_hits"]]
        for j, v in enumerate(vals):
            add_rich_runs(tbl2a.rows[i].cells[j].paragraphs[0], v, font_size=9.0)
    configure_three_line_table(tbl2a, [5.6, 1.3, 1.8, 1.7, 3.3, 2.3], {1, 2, 3, 4, 5})

    spacer_p = doc.add_paragraph()
    spacer_p.paragraph_format.space_before = Pt(4)
    spacer_p.paragraph_format.space_after = Pt(2)
    add_rich_runs(spacer_p, "**Table 2B.** Combined Discovery + GSE181919 shared-gene program effects.", font_size=9.5)

    comb_sel = [r for r in comb_rows if r["analysis"] in {"primary", "min_5"}]
    comb_order = {"primary": 0, "min_5": 1}
    prog_order_idx = {"Maturation": 0, "NF-kB feedback": 1, "Regulatory expression": 2, "MHC class II": 3}
    comb_sel.sort(key=lambda r: (comb_order[r["analysis"]], prog_order_idx[r["program"]]))

    t2b_headers = [
        "Model (combined n)",
        "Program",
        "Mean Δ log₂ CPM (95% t CI)",
        "Positive donors",
        "Sign-flip P",
        "Paired t P",
        "Fisher fry FDR",
    ]
    tbl2b = doc.add_table(rows=1 + len(comb_sel), cols=len(t2b_headers))
    for j, h in enumerate(t2b_headers):
        add_rich_runs(tbl2b.rows[0].cells[j].paragraphs[0], h, font_size=8.8, force_bold=True)
    for i, r in enumerate(comb_sel, start=1):
        m_lbl = f"Primary (5 + {r['n_external']} = {r['n_donors']})" if r["analysis"] == "primary" else f"Five-cell (5 + {r['n_external']} = {r['n_donors']})"
        p_lbl = r["program"].replace("NF-kB", "NF-κB")
        eff = f"{fmt_signed(float(r['mean_effect']))} ({fmt_signed(float(r['lower_95_t']))}, {fmt_signed(float(r['upper_95_t']))})"
        don = f"{r['positive_donors']}/{r['n_donors']}"
        vals = [m_lbl, p_lbl, eff, don, fmt_p(float(r["p_sign_flip"])), fmt_p(float(r["p_paired_t"])), fmt_p(float(r["fdr_fry_meta"]))]
        for j, v in enumerate(vals):
            add_rich_runs(tbl2b.rows[i].cells[j].paragraphs[0], v, font_size=8.8)
    configure_three_line_table(tbl2b, [3.2, 2.7, 3.4, 1.5, 1.7, 1.7, 1.8], {2, 3, 4, 5, 6})

    # ---- 4. Discussion ----
    add_heading_p("4. Discussion", level=1)
    for cp in cleaned_sections["4. Discussion"]:
        add_body_p(cp)

    # ---- Declarations & Data Availability ----
    add_heading_p("Declarations", level=1)
    add_body_p(
        "**Ethics approval and consent to participate:** Not applicable. This study is a secondary re-analysis "
        "of publicly deposited, de-identified human single-cell RNA sequencing datasets from the NCBI Gene "
        "Expression Omnibus (GSE172577, GSE193766, and GSE181919) [5,6,10]."
    )
    add_body_p("**Competing interests:** The authors declare no competing interests.")
    add_body_p("**Funding:** [Insert funding sources or state: This research received no external funding.]")
    add_body_p("**Authors' contributions:** [Insert CRediT author contribution statement.]")

    add_heading_p("Data and code availability", level=1)
    for idx, raw_p in enumerate(sections["Data and code availability"], start=1):
        cp = attach_supplementary_callout(strip_md_links_and_normalize(raw_p), f"avail_p{idx}")
        add_body_p(cp)

    # ---- References ----
    doc.add_page_break()
    add_heading_p("References", level=1)
    for num, body, journal, tail in REFERENCES:
        p = doc.add_paragraph()
        p.paragraph_format.line_spacing = 1.35
        p.paragraph_format.space_before = Pt(0)
        p.paragraph_format.space_after = Pt(4)
        p.paragraph_format.left_indent = Cm(0.75)
        p.paragraph_format.first_line_indent = Cm(-0.75)
        add_rich_runs(p, f"{num}. {body} *{journal}*. {tail}", font_size=10.0)

    # ---- Embedded Figures 1–5 ----
    fig_files = [
        (1, FIG_DIR / "Fig1_discovery_atlas.png"),
        (2, FIG_DIR / "Fig2_pseudobulk_de_and_robustness.png"),
        (3, FIG_DIR / "Fig3_programs_and_inflammation_paradox.png"),
        (4, FIG_DIR / "Fig4_lineage_and_tme_context.png"),
        (5, FIG_DIR / "Fig5_external_validation.png"),
    ]
    for fig_num, fig_path in fig_files:
        doc.add_page_break()
        add_heading_p(f"Figure {fig_num}", level=1)
        pic_p = doc.add_paragraph()
        pic_p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        pic_p.paragraph_format.space_before = Pt(4)
        pic_p.paragraph_format.space_after = Pt(8)
        pic_p.add_run().add_picture(str(fig_path), width=Inches(6.2))
        fig_title, fig_body = data["legends"][fig_num]
        cap_p = doc.add_paragraph()
        cap_p.paragraph_format.line_spacing = 1.35
        cap_p.paragraph_format.space_after = Pt(6)
        add_rich_runs(cap_p, f"**Figure {fig_num}. {fig_title}** {fig_body}", font_size=10.0)

    # ---- Supplementary Materials Index ----
    doc.add_page_break()
    add_heading_p("Supplementary Materials Index", level=1)
    add_body_p(
        f"All supplementary tables and figures listed below are archived in the companion repository ({REPO_URL}).",
        space_after=6,
    )
    sup_headers = ["Identifier", "Description", "Repository file(s)"]
    tbl_sup = doc.add_table(rows=1 + len(SUPPLEMENTARY_INDEX), cols=3)
    for j, h in enumerate(sup_headers):
        add_rich_runs(tbl_sup.rows[0].cells[j].paragraphs[0], h, font_size=9.0, force_bold=True)
    for i, (sid, sdesc, spath) in enumerate(SUPPLEMENTARY_INDEX, start=1):
        add_rich_runs(tbl_sup.rows[i].cells[0].paragraphs[0], sid, font_size=8.8, force_bold=True)
        add_rich_runs(tbl_sup.rows[i].cells[1].paragraphs[0], sdesc, font_size=8.8)
        add_rich_runs(tbl_sup.rows[i].cells[2].paragraphs[0], spath, font_size=8.5)
    configure_three_line_table(tbl_sup, [3.6, 6.4, 6.0], set())

    doc.save(OUT_DOCX)
    return abstract_words_exact, main_words


def tex_escape(text: str) -> str:
    rep = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    for k, v in rep.items():
        text = text.replace(k, v)
    return text


def to_tex_rich(text: str) -> str:
    parts = []
    for chunk, bold, italic, sup, sub in tokenize_rich_text(text):
        esc = tex_escape(chunk)
        esc = esc.replace("κ", r"$\kappa$").replace("ρ", r"$\rho$").replace("Δ", r"$\Delta$")
        esc = esc.replace("×", r"$\times$").replace("−", "$-$").replace("–", "--").replace("—", "---")
        if sup:
            esc = f"$^{{{esc}}}$"
        elif sub:
            esc = f"$_{{{esc}}}$"
        if bold and italic:
            esc = f"\\textbf{{\\textit{{{esc}}}}}"
        elif bold:
            esc = f"\\textbf{{{esc}}}"
        elif italic:
            esc = f"\\textit{{{esc}}}"
        parts.append(esc)
    return "".join(parts)


def build_tex(data: dict, t1_rows, ext_rows, comb_rows):
    sections = data["sections"]
    abs_items = {}
    cur_lbl = None
    for item in sections["Abstract"]:
        clean_item = strip_md_links_and_normalize(item)
        if clean_item.startswith("**") and clean_item.endswith("**"):
            cur_lbl = clean_item.strip("*")
        elif cur_lbl:
            abs_items[cur_lbl] = clean_item
    abs_items["Methods"] = attach_supplementary_callout(abs_items["Methods"], "abstract_methods")
    abs_items["Results"] = attach_supplementary_callout(abs_items["Results"], "abstract_results")

    out = [
        r"\documentclass[11pt,a4paper]{article}",
        r"\usepackage[margin=2.54cm]{geometry}",
        r"\usepackage{booktabs}",
        r"\usepackage{graphicx}",
        r"\usepackage{setspace}",
        r"\usepackage[hidelinks]{hyperref}",
        r"\usepackage{array}",
        r"\onehalfspacing",
        "",
        rf"\title{{{to_tex_rich(data['title'])}}}",
        r"\author{[Author Name(s)]\\ \small [Department, Institution, City, Country]\\ \small Corresponding author: [Email]; ORCID: [0000-0000-0000-0000]}",
        r"\date{}",
        "",
        r"\begin{document}",
        r"\maketitle",
        "",
        r"\begin{abstract}",
        rf"\noindent\textbf{{Background:}} {to_tex_rich(abs_items['Background'])}\par\medskip",
        rf"\noindent\textbf{{Methods:}} {to_tex_rich(abs_items['Methods'])}\par\medskip",
        rf"\noindent\textbf{{Results:}} {to_tex_rich(abs_items['Results'])}\par\medskip",
        rf"\noindent\textbf{{Conclusions:}} {to_tex_rich(abs_items['Conclusions'])}\par\medskip",
        rf"\noindent\textbf{{Keywords:}} {to_tex_rich(abs_items['Keywords'])}",
        r"\end{abstract}",
        "",
        r"\section{Introduction}",
    ]

    for raw_p in sections["1. Introduction"]:
        out.append(to_tex_rich(strip_md_links_and_normalize(raw_p)) + "\n")

    out.append(r"\section{Materials and methods}")
    for sec_name, tag in [
        ("2.1. Discovery cohorts, annotation and composition", "2.1"),
        ("2.2. Donor-paired pseudobulk models and sensitivity analyses", "2.2"),
        ("2.3. Held-out expression programs and competitive pathways", "2.3"),
        ("2.4. Patient-level tumour microenvironment covariation", "2.4"),
        ("2.5. Independent external validation and cross-cohort pooling", "2.5"),
    ]:
        sub_title = re.sub(r"^\d+\.\d+\.\s*", "", sec_name)
        out.append(rf"\subsection{{{to_tex_rich(sub_title)}}}")
        for idx, raw_p in enumerate(sections[sec_name], start=1):
            cp = attach_supplementary_callout(strip_md_links_and_normalize(raw_p), f"{tag}_p{idx}")
            out.append(to_tex_rich(cp) + "\n")

    out.append(r"\section{Results}")
    res_specs = [
        ("3.1. Discovery atlas, DC heterogeneity, and matched-tissue composition", "3.1", 1),
        ("3.2. Donor-paired differential expression and methodological robustness", "3.2", 2),
        ("3.3. Held-out transcriptional programs and resolution of the inflammatory paradox", "3.3", 3),
        ("3.4. DC lineage specificity including mucosal LC-like states and TME covariation", "3.4", 4),
        ("3.5. Independent validation in GSE181919 and cross-cohort meta-analysis", "3.5", 5),
    ]
    fig_pdfs = {
        1: "../results/figures/Fig1_discovery_atlas.pdf",
        2: "../results/figures/Fig2_pseudobulk_de_and_robustness.pdf",
        3: "../results/figures/Fig3_programs_and_inflammation_paradox.pdf",
        4: "../results/figures/Fig4_lineage_and_tme_context.pdf",
        5: "../results/figures/Fig5_external_validation.pdf",
    }

    for sec_name, tag, fig_num in res_specs:
        sub_title = re.sub(r"^\d+\.\d+\.\s*", "", sec_name)
        out.append(rf"\subsection{{{to_tex_rich(sub_title)}}}")
        for idx, raw_p in enumerate(sections[sec_name], start=1):
            cp = attach_supplementary_callout(strip_md_links_and_normalize(raw_p), f"{tag}_p{idx}")
            out.append(to_tex_rich(cp) + "\n")
        if fig_num == 3:
            out.extend([
                r"\begin{table}[htbp]",
                r"\centering",
                r"\caption{Held-out transcriptional expression programs in the discovery primary paired comparison (migDC versus cDC2; $n = 5$ matched tumour donors).}",
                r"\scriptsize",
                r"\begin{tabular}{lccccccc}",
                r"\toprule",
                r"Held-out program & Genes & Mean $\Delta \log_2$ CPM (95\% $t$ CI) & Positive donors & Sign-flip $P$ & Paired $t$ $P$ & Paired $t$ FDR & Directional fry FDR \\",
                r"\midrule",
            ])
            for disp_name, cov, eff_ci, donors, p_flip, p_t, fdr_t, fdr_fry, r in t1_rows:
                out.append(
                    f"{to_tex_rich(disp_name)} & {cov} & {to_tex_rich(eff_ci)} & {donors} & "
                    f"{fmt_p_tex(float(r['p_sign_flip']))} & {fmt_p_tex(float(r['p_paired_t']))} & "
                    f"{fmt_p_tex(float(r['fdr_paired_t']))} & {fmt_p_tex(float(r['fdr_fry']))} \\\\"
                )
            out.extend([
                r"\bottomrule",
                r"\end{tabular}",
                r"\end{table}",
                "",
            ])
        if fig_num == 5:
            ext_labels = {
                "primary": "Myeloid APC pool (min 10 cells)",
                "min_5": "Myeloid APC pool (min 5 cells)",
                "cd163_neg_min10": "CD163-negative myeloid (min 10)",
                "cd163_neg_min5": "CD163-negative myeloid (min 5)",
                "hybrid_author_migDC_vs_cd163neg_cDC2": "Hybrid author migDC vs CD163-neg cDC2",
                "author_dc_migDC_vs_cDC1_min3": "Author DC migDC vs cDC1 (min 3)",
                "migDC_vs_cDC1": "Myeloid APC migDC vs cDC1 (min 5)",
            }
            out.extend([
                r"\begin{table}[htbp]",
                r"\centering",
                r"\caption{Independent validation in GSE181919 and combined cross-cohort shared-gene program meta-analysis.}",
                r"\scriptsize",
                r"\begin{tabular}{lccccc}",
                r"\toprule",
                r"\multicolumn{6}{l}{\textbf{A. External gene-level replication}} \\",
                r"\midrule",
                r"External model (contrast) & Donors & Shared genes & Spearman $\rho$ & Concordant discovery hits & External FDR $<$ 0.05 hits \\",
                r"\midrule",
            ])
            for r in ext_rows:
                lbl = ext_labels.get(r["analysis"], r["analysis"])
                hits = f"{r['concordant_hits']}/{r['n_discovery_hits']} ({100*float(r['concordant_hit_rate']):.1f}\\%)"
                out.append(
                    f"{to_tex_rich(lbl)} & {r['n_donors']} & {int(r['n_shared']):,} & "
                    f"{float(r['rho']):.3f} & {hits} & {r['external_fdr_hits']} \\\\"
                )
            out.extend([
                r"\midrule",
                r"\multicolumn{6}{l}{\textbf{B. Combined Discovery ($n=5$) + GSE181919 shared-gene program effects}} \\",
                r"\midrule",
                r"Model (combined $n$) & Program & Mean $\Delta \log_2$ CPM (95\% $t$ CI) & Positive donors & Sign-flip $P$ / Paired $t$ $P$ & Fisher fry FDR \\",
                r"\midrule",
            ])
            comb_sel = [r for r in comb_rows if r["analysis"] in {"primary", "min_5"}]
            comb_order = {"primary": 0, "min_5": 1}
            prog_order_idx = {"Maturation": 0, "NF-kB feedback": 1, "Regulatory expression": 2, "MHC class II": 3}
            comb_sel.sort(key=lambda r: (comb_order[r["analysis"]], prog_order_idx[r["program"]]))
            for r in comb_sel:
                m_lbl = f"Primary (5 + {r['n_external']} = {r['n_donors']})" if r["analysis"] == "primary" else f"Five-cell (5 + {r['n_external']} = {r['n_donors']})"
                p_lbl = r["program"].replace("NF-kB", "NF-κB")
                eff = f"{fmt_signed(float(r['mean_effect']))} ({fmt_signed(float(r['lower_95_t']))}, {fmt_signed(float(r['upper_95_t']))})"
                don = f"{r['positive_donors']}/{r['n_donors']}"
                p_pair = f"{fmt_p_tex(float(r['p_sign_flip']))} / {fmt_p_tex(float(r['p_paired_t']))}"
                out.append(
                    f"{to_tex_rich(m_lbl)} & {to_tex_rich(p_lbl)} & {to_tex_rich(eff)} & {don} & "
                    f"{p_pair} & {fmt_p_tex(float(r['fdr_fry_meta']))} \\\\"
                )
            out.extend([r"\bottomrule", r"\end{tabular}", r"\end{table}", ""])

        fig_title, fig_body = data["legends"][fig_num]
        out.extend([
            r"\begin{figure}[htbp]",
            r"\centering",
            rf"\includegraphics[width=\textwidth]{{{fig_pdfs[fig_num]}}}",
            rf"\caption{{\textbf{{{to_tex_rich(fig_title)}}} {to_tex_rich(fig_body)}}}",
            r"\end{figure}",
            "",
        ])

    out.append(r"\section{Discussion}")
    for idx, raw_p in enumerate(sections["4. Discussion"], start=1):
        cp = attach_supplementary_callout(strip_md_links_and_normalize(raw_p), f"4_p{idx}")
        out.append(to_tex_rich(cp) + "\n")

    out.extend([
        r"\section*{Declarations}",
        r"\noindent\textbf{Ethics approval and consent to participate:} Not applicable. Secondary analysis of publicly deposited, de-identified GEO datasets (GSE172577, GSE193766, GSE181919) [5,6,10].\par\medskip",
        r"\noindent\textbf{Competing interests:} The authors declare no competing interests.\par\medskip",
        r"\noindent\textbf{Funding:} [Insert funding sources or state: This research received no external funding.]\par\medskip",
        r"\noindent\textbf{Authors' contributions:} [Insert CRediT author contribution statement.]",
        "",
        r"\section*{Data and code availability}",
    ])
    for idx, raw_p in enumerate(sections["Data and code availability"], start=1):
        cp = attach_supplementary_callout(strip_md_links_and_normalize(raw_p), f"avail_p{idx}")
        out.append(to_tex_rich(cp) + "\n")

    out.append(r"\begin{thebibliography}{12}")
    for num, body, journal, tail in REFERENCES:
        out.append(rf"\bibitem{{ref{num}}} {tex_escape(body)} \textit{{{tex_escape(journal)}}}. {tex_escape(tail)}")
    out.append(r"\end{thebibliography}")
    out.append("")
    out.append(r"\section*{Supplementary Materials Index}")
    out.extend([
        r"\begin{table}[htbp]",
        r"\centering",
        r"\scriptsize",
        r"\begin{tabular}{p{3.2cm}p{5.8cm}p{5.8cm}}",
        r"\toprule",
        r"Identifier & Description & Repository file(s) \\",
        r"\midrule",
    ])
    for sid, sdesc, spath in SUPPLEMENTARY_INDEX:
        out.append(f"{to_tex_rich(sid)} & {to_tex_rich(sdesc)} & \\texttt{{{tex_escape(spath)}}} \\\\")
    out.extend([r"\bottomrule", r"\end{tabular}", r"\end{table}", "", r"\end{document}"])

    OUT_TEX.write_text("\n".join(out), encoding="utf-8", newline="\n")


def main():
    data = parse_manuscript()
    t1_rows, ext_rows, comb_rows = build_live_tables()
    abs_w, main_w = build_docx(data, t1_rows, ext_rows, comb_rows)
    build_tex(data, t1_rows, ext_rows, comb_rows)
    print(f"Exported {OUT_DOCX.relative_to(ROOT)} (Abstract: {abs_w} words, Main: {main_w} words)")
    print(f"Exported {OUT_TEX.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
