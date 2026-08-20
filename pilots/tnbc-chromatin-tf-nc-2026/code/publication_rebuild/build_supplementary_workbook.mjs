#!/usr/bin/env node

import fs from "node:fs/promises";
import path from "node:path";
import zlib from "node:zlib";
import { promisify } from "node:util";
import { FileBlob, Workbook, SpreadsheetFile } from "@oai/artifact-tool";

if (process.argv[2] === "--render-existing") {
  const workbookPath = process.argv[3];
  const qaPath = process.argv[4];
  const sheetName = process.argv[5];
  if (!workbookPath || !qaPath || !sheetName) {
    throw new Error("Usage: build_supplementary_workbook.mjs --render-existing WORKBOOK QA_DIR SHEET_NAME");
  }
  const input = await FileBlob.load(workbookPath);
  const imported = await SpreadsheetFile.importXlsx(input);
  const preview = await imported.render({ sheetName, range: "A1:J35", scale: 1, format: "png" });
  await fs.mkdir(qaPath, { recursive: true });
  const bytes = new Uint8Array(await preview.arrayBuffer());
  await fs.writeFile(path.join(qaPath, `${sheetName}.png`), bytes);
  console.log(JSON.stringify({ workbookPath, sheetName, qaPath }));
  process.exit(0);
}

const gunzip = promisify(zlib.gunzip);
const inputRoot = process.argv[2];
const outputPath = process.argv[3];
const qaDir = process.argv[4];
if (!inputRoot || !outputPath || !qaDir) {
  throw new Error("Usage: build_supplementary_workbook.mjs INPUT_ROOT OUTPUT_XLSX QA_DIR");
}

const navy = "#17365D";
const blue = "#2878B5";
const paleBlue = "#DCEAF5";
const paleGold = "#FFF2CC";
const paleRed = "#FCE4D6";
const paleGreen = "#E2F0D9";
const gray = "#E7E9EB";

async function walk(dir) {
  const out = [];
  for (const entry of await fs.readdir(dir, { withFileTypes: true })) {
    const full = path.join(dir, entry.name);
    if (entry.isDirectory()) out.push(...await walk(full));
    else out.push(full);
  }
  return out;
}

function scalar(value) {
  if (value === "" || value === "NA" || value === "NaN") return null;
  if (value === "TRUE") return true;
  if (value === "FALSE") return false;
  if (/^-?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?$/.test(value)) {
    const number = Number(value);
    if (Number.isFinite(number)) return number;
  }
  return value;
}

async function readTsv(file) {
  let bytes = await fs.readFile(file);
  if (file.endsWith(".gz")) bytes = await gunzip(bytes);
  const text = bytes.toString("utf8").replace(/^\uFEFF/, "");
  const lines = text.split(/\r?\n/).filter((line) => line.length > 0);
  if (!lines.length) return [];
  return lines.map((line) => line.split("\t").map(scalar));
}

const allFiles = await walk(inputRoot);
const byBase = new Map();
for (const file of allFiles) {
  const base = path.basename(file);
  if (!byBase.has(base)) byBase.set(base, file);
}

function find(base) {
  const file = byBase.get(base);
  if (!file) throw new Error(`Missing workbook source: ${base}`);
  return file;
}

function excelColumn(n) {
  let answer = "";
  let x = n;
  while (x > 0) {
    x -= 1;
    answer = String.fromCharCode(65 + (x % 26)) + answer;
    x = Math.floor(x / 26);
  }
  return answer;
}

function styleTitle(range) {
  range.format.fill = navy;
  range.format.font = { bold: true, color: "#FFFFFF", size: 15 };
  range.format.rowHeight = 30;
  range.format.verticalAlignment = "center";
}

function styleSection(range) {
  range.format.fill = blue;
  range.format.font = { bold: true, color: "#FFFFFF", size: 10 };
  range.format.rowHeight = 22;
  range.format.verticalAlignment = "center";
}

function styleHeader(range) {
  range.format.fill = paleBlue;
  range.format.font = { bold: true, color: navy };
  range.format.borders = { preset: "all", style: "thin", color: "#B8C5D1" };
  range.format.wrapText = true;
  range.format.verticalAlignment = "center";
}

function styleBody(range) {
  range.format.borders = { preset: "all", style: "thin", color: "#D9DEE3" };
  range.format.verticalAlignment = "top";
}

function writeMatrix(sheet, startRow, sectionName, matrix, maxRows = 10000, maxCols = 60) {
  if (!matrix.length) return startRow;
  const rows = matrix.slice(0, Math.min(matrix.length, maxRows + 1));
  const cols = Math.min(Math.max(...rows.map((r) => r.length)), maxCols);
  const normalized = rows.map((row) => {
    const z = row.slice(0, cols);
    while (z.length < cols) z.push(null);
    return z;
  });
  const endCol = excelColumn(cols);
  const section = sheet.getRange(`A${startRow}:${endCol}${startRow}`);
  section.merge();
  section.values = [[`${sectionName} | ${matrix.length - 1} data rows${matrix.length > rows.length ? `; workbook preview limited to ${maxRows}` : ""}`]];
  styleSection(section);
  const headerRow = startRow + 1;
  const bodyEnd = headerRow + normalized.length - 1;
  const target = sheet.getRange(`A${headerRow}:${endCol}${bodyEnd}`);
  target.values = normalized;
  styleHeader(sheet.getRange(`A${headerRow}:${endCol}${headerRow}`));
  if (bodyEnd > headerRow) styleBody(sheet.getRange(`A${headerRow + 1}:${endCol}${bodyEnd}`));
  return bodyEnd + 2;
}

function createSheet(workbook, name, subtitle) {
  const sheet = workbook.worksheets.add(name);
  sheet.showGridLines = false;
  sheet.getRange("A1:J2").merge();
  sheet.getRange("A1:J2").values = [[name]];
  styleTitle(sheet.getRange("A1:J2"));
  sheet.getRange("A3:J4").merge();
  sheet.getRange("A3:J4").values = [[subtitle]];
  sheet.getRange("A3:J4").format.fill = "#F4F6F8";
  sheet.getRange("A3:J4").format.font = { italic: true, color: "#3C4752", size: 9 };
  sheet.getRange("A3:J4").format.wrapText = true;
  sheet.freezePanes.freezeRows(4);
  return sheet;
}

const workbook = Workbook.create();
const overview = createSheet(
  workbook,
  "Overview",
  "Navigation and frozen scientific conclusions. Full-resolution TSV/TSV.GZ files remain beside this workbook; the workbook is a readable evidence index."
);
overview.getRange("A6:D6").values = [["Evidence layer", "Frozen result", "Interpretation", "Status"]];
styleHeader(overview.getRange("A6:D6"));
const overviewRows = [
  ["Figure 1 discovery", "158 LUAD-specific TFs", "TCGA LUAD versus LUSC ARACNe3/msVIPER", "Supported"],
  ["Independent patients", "95/158 GSE81089; 70/158 GSE41271; 44 in both", "RNA-seq plus cross-platform microarray validation", "Supported, smaller than METABRIC"],
  ["Figure 2 HC gate", "31 HC-TFs", "Promoter open in >=half of each system plus anchor activity rule", "Supported"],
  ["Strict motif core", "FOXA3; NFATC4; XBP1", "HOMER q<1e-5 and >=half of samples in all three systems", "Selective, not broad"],
  ["Motif sensitivity", "ETV1; FOXA3; NFATC4; XBP1; ZNF75D", "q<=0.05 and >=half of samples", "Sensitivity only"],
  ["Figure 3 network", "31 HC-TFs; 3,011 signed edges; 2,776 targets", "All HC-TFs enter downstream analyses", "Supported"],
  ["Figure 4 prognosis", "Only TCGA OS–ZNF75D passes frozen FDR", "Most associations remain nominal/exploratory", "Weak"],
  ["Figure 5 cell lines", "Five replicated drug-TF pairs; all resistance direction", "Consistent in >=2 resources", "Replicated in vitro"],
  ["Figure 5 PDX", "0 validated pairs", "Matched public PDX response evidence is absent/insufficient", "Not validated in vivo"],
];
overview.getRange(`A7:D${6 + overviewRows.length}`).values = overviewRows;
styleBody(overview.getRange(`A7:D${6 + overviewRows.length}`));
overview.getRange(`D7:D${6 + overviewRows.length}`).conditionalFormats.add("containsText", {
  text: "Supported", format: { fill: paleGreen, font: { color: "#375623", bold: true } },
});
overview.getRange(`D7:D${6 + overviewRows.length}`).conditionalFormats.add("containsText", {
  text: "Weak", format: { fill: paleGold, font: { color: "#7F6000", bold: true } },
});
overview.getRange(`D7:D${6 + overviewRows.length}`).conditionalFormats.add("containsText", {
  text: "Not validated", format: { fill: paleRed, font: { color: "#9C0006", bold: true } },
});
overview.getRange("A:D").format.columnWidth = 26;
overview.getRange("B:B").format.columnWidth = 36;
overview.getRange("C:C").format.columnWidth = 54;
overview.getRange("D:D").format.columnWidth = 25;
overview.getRange("A6:D20").format.wrapText = true;

const sections = [
  {
    sheet: "Data01_Cohorts",
    subtitle: "Cohort manifests, histology, clinical annotations, and inclusion counts.",
    files: [
      ["SupplementaryFigure1_inclusion_counts.tsv", "Patient inclusion counts", 500],
      ["tcga_manifest.tsv", "TCGA discovery manifest", 1200],
      ["gse81089_manifest.tsv", "GSE81089 validation manifest", 500],
      ["gse41271_manifest.tsv", "GSE41271 validation manifest", 500],
      ["pdmr_manifest.tsv", "PDMR PDX manifest", 500],
      ["depmap_manifest.tsv", "DepMap cell-line manifest", 500],
    ],
  },
  {
    sheet: "Data02_TF_Discovery",
    subtitle: "Differential TF activity, independent reconstruction, and replicated effects.",
    files: [
      ["TCGA_msviper.tsv", "TCGA msVIPER TF screen", 5000],
      ["GSE81089_msviper.tsv", "GSE81089 msVIPER TF screen", 5000],
      ["GSE41271_msviper.tsv", "GSE41271 msVIPER TF screen", 5000],
      ["SupplementaryFigure2B_TF_overlap_counts.tsv", "TF overlap counts", 500],
      ["SupplementaryFigure2F_replicated_LUAD_TF_effect_zscores.tsv", "Replicated TF effects", 5000],
    ],
  },
  {
    sheet: "Data03_HC_Network",
    subtitle: "Frozen 158-to-31 gate and the signed HC-TF regulatory network.",
    files: [
      ["tf_promoter_gate_summary.tsv", "Discovery-TF promoter/activity gate", 1000],
      ["system_tf_promoter_activity_summary.tsv", "System-level gate evidence", 2000],
      ["Figure3B_target_gene_network_nodes.tsv", "Target-gene community nodes", 5000],
      ["Figure3B_target_gene_network_edges.tsv", "Target-gene community edges", 5000],
      ["Figure3B_target_gene_community_GO.tsv", "Community GO enrichment", 5000],
    ],
  },
  {
    sheet: "Data04_ATAC_QC",
    subtitle: "Independent-unit manifest, raw-QC receipts, saturation, correlations, and genomic annotation.",
    files: [
      ["figure2_atomic_sample_manifest.tsv", "Figure 2 sample manifest", 1000],
      ["figure2_raw_atac_qc.tsv", "Raw ATAC QC", 1000],
      ["saturation_summary.tsv", "Saturation summary", 2000],
      ["sample_accessibility_pearson_correlations.tsv", "Pairwise accessibility correlations", 10000],
      ["peak_genomic_annotation.tsv", "Peak genomic annotation", 10000],
    ],
  },
  {
    sheet: "Data05_Motif",
    subtitle: "Promoter evidence, JASPAR-priority/CIS-BP-fallback mapping, and strict/sensitivity motif results.",
    files: [
      ["anchor_consensus_hc_tf_triple_system_summary.tsv", "HC-TF three-system motif summary", 1000],
      ["anchor_consensus_system_tf_motif_summary.tsv", "System-level motif summary", 5000],
      ["motif_threshold_sensitivity_triple_system_tfs.tsv", "Threshold sensitivity", 1000],
    ],
  },
  {
    sheet: "Data06_Survival",
    subtitle: "All Cox models, endpoint overlaps, Kaplan-Meier receipts, and 5,000-permutation summaries.",
    files: [
      ["Figure4_all_Cox_models.tsv", "All continuous-activity Cox models", 5000],
      ["Figure4EG_endpoint_overlap.tsv", "Endpoint overlap", 1000],
      ["Figure4FH_KM_receipts.tsv", "Kaplan-Meier receipts", 1000],
      ["SupplementaryFigure7_permutation_summary.tsv", "Permutation summary", 1000],
    ],
  },
  {
    sheet: "Data07_Drugs",
    subtitle: "Dataset audit, drug identity, all cell-line screens, and the five replicated pairs.",
    files: [
      ["Figure5_dataset_audit.tsv", "Dataset audit", 1000],
      ["Figure5_replicated_pair_dataset_details.tsv", "Replicated pair details", 1000],
      ["SupplementaryFigure7B_scatter_receipts.tsv", "Scatterplot receipts", 1000],
      ["SupplementaryFigure7A_drug_membership.tsv", "Drug memberships", 5000],
      ["SupplementaryFigure7C_all_significant_associations.tsv.gz", "All significant associations", 10000],
    ],
  },
  {
    sheet: "Data08_PDX",
    subtitle: "PDX classifier state, drug coverage, all evaluable concordance tests, and the explicit zero-validation result.",
    files: [
      ["Figure5G_PDX_validation_funnel.tsv", "PDX validation funnel", 1000],
      ["Figure5_PDXE_v2_lung_state_scores.tsv", "PDX state scores", 5000],
      ["SupplementaryFigure9A_PDXE_drug_coverage.tsv", "PDX drug coverage", 5000],
      ["SupplementaryFigure9B_PDXE_all_concordance_tests.tsv", "All PDX concordance tests", 10000],
      ["Figure5GH_PDX_validated_pairs.tsv", "Validated pairs (zero rows is a result)", 1000],
    ],
  },
  {
    sheet: "Data09_Provenance",
    subtitle: "Cell-line/PDX identifiers, tissue eligibility, mapping methods, and classifier provenance.",
    files: [
      ["SupplementaryFigure8A_cell_mapping.tsv", "Cell-line mapping", 5000],
      ["SupplementaryFigure8B_drug_eligibility.tsv", "Drug eligibility", 5000],
      ["SupplementaryFigure8C_lung_classifier_validation.tsv", "Classifier validation", 5000],
      ["SupplementaryFigure8E_PDMR_classifier_scores.tsv", "PDMR classifier scores", 5000],
      ["SupplementaryFigure8F_DepMap_classifier_scores.tsv", "DepMap classifier scores", 5000],
    ],
  },
];

for (const spec of sections) {
  const sheet = createSheet(workbook, spec.sheet, spec.subtitle);
  let row = 6;
  for (const [fileName, label, maxRows] of spec.files) {
    const matrix = await readTsv(find(fileName));
    row = writeMatrix(sheet, row, label, matrix, maxRows, 60);
  }
  const used = sheet.getUsedRange();
  used.format.wrapText = false;
  used.format.autofitColumns();
  used.format.autofitRows();
  const colCount = Math.min(used.columnCount ?? 12, 60);
  for (let c = 0; c < colCount; c += 1) {
    const col = sheet.getRangeByIndexes(0, c, Math.max(1, used.rowCount ?? 1), 1);
    if ((col.format.columnWidth ?? 0) > 40) col.format.columnWidth = 40;
  }
}

const indexSheet = createSheet(workbook, "Files_Index", "Full-resolution file inventory with provenance and checksums.");
const indexMatrix = await readTsv(find("Supplementary_Data_1-9_file_index.tsv"));
writeMatrix(indexSheet, 6, "Supplementary Data file inventory", indexMatrix, 10000, 20);
indexSheet.getUsedRange().format.autofitColumns();
indexSheet.getUsedRange().format.autofitRows();

// Deterministic pre-export error scan. No formulas are used, so any Excel error
// token would necessarily have entered as source text and is reported here.
const errorTokens = ["#REF!", "#DIV/0!", "#VALUE!", "#NAME?", "#NUM!"];
for (const file of allFiles.filter((f) => /\.tsv(?:\.gz)?$/.test(f))) {
  const matrix = await readTsv(file);
  for (const row of matrix.slice(0, 10001)) {
    for (const value of row) {
      if (typeof value === "string" && errorTokens.some((token) => value.includes(token))) {
        throw new Error(`Spreadsheet error token in ${file}: ${value}`);
      }
    }
  }
}

await fs.mkdir(path.dirname(outputPath), { recursive: true });
await fs.mkdir(qaDir, { recursive: true });
const exported = await SpreadsheetFile.exportXlsx(workbook);
await exported.save(outputPath);

const inspection = await workbook.inspect({
  kind: "workbook,sheet,table",
  maxChars: 12000,
  tableMaxRows: 5,
  tableMaxCols: 8,
  tableMaxCellChars: 80,
});
await fs.writeFile(path.join(qaDir, "workbook_inspection.ndjson"), inspection.ndjson ?? String(inspection), "utf8");

const errorScan = await workbook.inspect({
  kind: "match",
  searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A",
  options: { useRegex: true, maxResults: 300 },
  summary: "final formula error scan",
});
await fs.writeFile(path.join(qaDir, "workbook_error_scan.ndjson"), errorScan.ndjson ?? String(errorScan), "utf8");

if (process.env.SKIP_WORKBOOK_RENDER !== "1") {
  for (const sheet of workbook.worksheets.items) {
    const preview = await workbook.render({ sheetName: sheet.name, range: "A1:J35", scale: 1, format: "png" });
    const bytes = new Uint8Array(await preview.arrayBuffer());
    await fs.writeFile(path.join(qaDir, `${sheet.name}.png`), bytes);
  }
}

console.log(JSON.stringify({ outputPath, sheets: workbook.worksheets.items.map((s) => s.name), qaDir }, null, 2));
