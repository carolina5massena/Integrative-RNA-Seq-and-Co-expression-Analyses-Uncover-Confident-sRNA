"""
sRNA Pipeline GUI - PySide6 desktop app for the Integrative RNA-Seq workflow.

Tabs: 01 Position Classification, 02 Co-expression (WGCNA),
03 Differential Expression (DESeq2), 04 Filters.

Run:  pip install PySide6 pandas numpy openpyxl ;  python app.py
"""

import os
import sys
import traceback

from PySide6.QtCore import QObject, QThread, Signal
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QTabWidget, QWidget, QVBoxLayout, QHBoxLayout,
    QFormLayout, QGroupBox, QLabel, QLineEdit, QPushButton, QPlainTextEdit,
    QFileDialog, QComboBox, QSpinBox, QDoubleSpinBox, QMessageBox, QScrollArea,
)

import pipeline_logic as pl


class Worker(QObject):
    log = Signal(str)
    finished = Signal(bool, str)

    def __init__(self, func, params):
        super().__init__()
        self.func = func
        self.params = params

    def run(self):
        try:
            out = self.func(self.params, lambda m: self.log.emit(str(m)))
            self.finished.emit(True, f"Done. Output: {out}")
        except Exception as e:
            self.log.emit("ERROR:\n" + traceback.format_exc())
            self.finished.emit(False, str(e))


def browse_file_row(line_edit, parent, save=False, filt="All files (*.*)", multi=False):
    def start_dir():
        existing = [x.strip() for x in line_edit.text().split(",") if x.strip()]
        return os.path.dirname(existing[-1]) if existing else line_edit.text()

    def pick():
        if save:
            path, _ = QFileDialog.getSaveFileName(parent, "Select output file", line_edit.text(), filt)
            if path:
                line_edit.setText(path)
        elif multi:
            paths, _ = QFileDialog.getOpenFileNames(parent, "Select one or more files", start_dir(), filt)
            if paths:
                existing = [x.strip() for x in line_edit.text().split(",") if x.strip()]
                for p in paths:
                    if p not in existing:
                        existing.append(p)
                line_edit.setText(", ".join(existing))
        else:
            path, _ = QFileDialog.getOpenFileName(parent, "Select file", line_edit.text(), filt)
            if path:
                line_edit.setText(path)
    btn = QPushButton("Browse...")
    btn.clicked.connect(pick)
    row = QHBoxLayout()
    row.addWidget(line_edit)
    row.addWidget(btn)
    if multi:
        clear_btn = QPushButton("Clear")
        clear_btn.clicked.connect(lambda: line_edit.clear())
        row.addWidget(clear_btn)
    w = QWidget()
    w.setLayout(row)
    return w


def browse_dir_row(line_edit, parent):
    def pick():
        path = QFileDialog.getExistingDirectory(parent, "Select folder", line_edit.text())
        if path:
            line_edit.setText(path)
    btn = QPushButton("Browse...")
    btn.clicked.connect(pick)
    row = QHBoxLayout()
    row.addWidget(line_edit)
    row.addWidget(btn)
    w = QWidget()
    w.setLayout(row)
    return w


def scrollable(inner_widget):
    sa = QScrollArea()
    sa.setWidgetResizable(True)
    sa.setWidget(inner_widget)
    return sa


class AnalysisTab(QWidget):
    def __init__(self, main):
        super().__init__()
        self.main = main
        self.thread = None
        self.worker = None
        self.run_btn = None

    def build_log(self):
        self.log_view = QPlainTextEdit()
        self.log_view.setReadOnly(True)
        self.log_view.setPlaceholderText("Run output appears here...")
        box = QGroupBox("Log")
        lay = QVBoxLayout()
        lay.addWidget(self.log_view)
        box.setLayout(lay)
        return box

    def append_log(self, msg):
        self.log_view.appendPlainText(msg)

    def run_async(self, func, params):
        if self.thread and self.thread.isRunning():
            QMessageBox.warning(self, "Busy", "A task is already running.")
            return
        self.log_view.clear()
        if self.run_btn:
            self.run_btn.setEnabled(False)
        self.thread = QThread()
        self.worker = Worker(func, params)
        self.worker.moveToThread(self.thread)
        self.thread.started.connect(self.worker.run)
        self.worker.log.connect(self.append_log)
        self.worker.finished.connect(self._on_finished)
        self.worker.finished.connect(self.thread.quit)
        self.thread.start()

    def _on_finished(self, success, message):
        if self.run_btn:
            self.run_btn.setEnabled(True)
        if success:
            QMessageBox.information(self, "Finished", message)
        else:
            QMessageBox.critical(self, "Failed", message)


class PositionTab(AnalysisTab):
    def __init__(self, main):
        super().__init__(main)
        inner = QWidget()
        form = QFormLayout(inner)
        self.gtf_file = QLineEdit()
        form.addRow("GTF file:", browse_file_row(self.gtf_file, self,
                    filt="GTF (*.gtf *.gff *.txt);;All (*.*)"))
        self.gtf_feature = QLineEdit("exon")
        form.addRow("GTF feature type:", self.gtf_feature)
        self.gtf_id_field = QLineEdit("exon_id")
        form.addRow("GTF attribute key (ID):", self.gtf_id_field)
        self.srna_prefix = QLineEdit("srn_")
        form.addRow("sRNA ID prefix:", self.srna_prefix)
        self.gene_prefix = QLineEdit("gene")
        form.addRow("Gene ID prefix:", self.gene_prefix)
        self.srna_rename = QLineEdit("sRNA")
        form.addRow("sRNA rename prefix:", self.srna_rename)
        self.gene_rename = QLineEdit("Gene")
        form.addRow("Gene rename prefix:", self.gene_rename)
        self.overlap = QDoubleSpinBox()
        self.overlap.setRange(0.0, 1.0)
        self.overlap.setSingleStep(0.05)
        self.overlap.setValue(0.90)
        form.addRow("Overlap threshold (fraction):", self.overlap)
        self.utr = QSpinBox()
        self.utr.setRange(0, 100000)
        self.utr.setValue(150)
        form.addRow("UTR window distance (bp):", self.utr)
        self.output_file = QLineEdit(os.path.join(os.getcwd(), "sRNA_annotation.xlsx"))
        form.addRow("Output file (.xlsx/.csv):", browse_file_row(self.output_file, self,
                    save=True, filt="Excel (*.xlsx);;CSV (*.csv)"))
        self.run_btn = QPushButton("Run Position Classification")
        self.run_btn.clicked.connect(self.run)
        form.addRow(self.run_btn)
        lay = QVBoxLayout(self)
        lay.addWidget(scrollable(inner))
        lay.addWidget(self.build_log())

    def run(self):
        params = {
            "gtf_file": self.gtf_file.text().strip(),
            "gtf_feature": self.gtf_feature.text().strip(),
            "gtf_id_field": self.gtf_id_field.text().strip(),
            "srna_prefix": self.srna_prefix.text(),
            "gene_prefix": self.gene_prefix.text(),
            "srna_rename_prefix": self.srna_rename.text(),
            "gene_rename_prefix": self.gene_rename.text(),
            "overlap_threshold": self.overlap.value(),
            "utr_distance": self.utr.value(),
            "output_file": self.output_file.text().strip(),
        }
        self.run_async(pl.run_position_classification, params)


class WGCNATab(AnalysisTab):
    def __init__(self, main):
        super().__init__(main)
        inner = QWidget()
        form = QFormLayout(inner)
        self.expr = QLineEdit()
        form.addRow("Expression matrix TSV:", browse_file_row(self.expr, self,
                    filt="TSV (*.tsv *.txt);;All (*.*)"))
        self.traits = QLineEdit()
        form.addRow("Trait table CSV:", browse_file_row(self.traits, self,
                    filt="CSV (*.csv);;All (*.*)"))
        self.gene_id_col = QLineEdit("Geneid")
        form.addRow("Gene ID column:", self.gene_id_col)
        self.trait_sample_col = QLineEdit("rownames")
        form.addRow("Trait sample-ID column:", self.trait_sample_col)
        self.treatment_col = QLineEdit("treatment")
        form.addRow("Treatment column:", self.treatment_col)
        self.out_dir = QLineEdit(os.getcwd())
        form.addRow("Output directory:", browse_dir_row(self.out_dir, self))
        self.network_type = QComboBox()
        self.network_type.addItems(["unsigned", "signed"])
        form.addRow("Network type:", self.network_type)
        self.corfnc = QComboBox()
        self.corfnc.addItems(["bicor", "cor"])
        form.addRow("Correlation function:", self.corfnc)
        self.soft_power = QSpinBox()
        self.soft_power.setRange(1, 30)
        self.soft_power.setValue(10)
        form.addRow("Soft-threshold power:", self.soft_power)
        self.min_module = QSpinBox()
        self.min_module.setRange(2, 1000)
        self.min_module.setValue(30)
        form.addRow("Min module size:", self.min_module)
        self.deep_split = QSpinBox()
        self.deep_split.setRange(0, 4)
        self.deep_split.setValue(2)
        form.addRow("Deep split:", self.deep_split)
        self.medissthres = QDoubleSpinBox()
        self.medissthres.setRange(0.0, 1.0)
        self.medissthres.setSingleStep(0.05)
        self.medissthres.setValue(0.25)
        form.addRow("Module merge cut height (MEDissThres):", self.medissthres)
        self.cyt_thr = QDoubleSpinBox()
        self.cyt_thr.setRange(0.0, 1.0)
        self.cyt_thr.setDecimals(3)
        self.cyt_thr.setSingleStep(0.005)
        self.cyt_thr.setValue(0.02)
        form.addRow("Cytoscape edge threshold:", self.cyt_thr)
        self.cut_height = QSpinBox()
        self.cut_height.setRange(0, 100000)
        self.cut_height.setValue(140)
        form.addRow("Sample tree cut height (line only):", self.cut_height)
        self.threads = QLineEdit()
        self.threads.setPlaceholderText("blank = auto")
        form.addRow("Threads:", self.threads)
        self.rscript = QLineEdit(pl.find_rscript())
        form.addRow("Rscript executable:", browse_file_row(self.rscript, self))
        self.run_btn = QPushButton("Run WGCNA")
        self.run_btn.clicked.connect(self.run)
        form.addRow(self.run_btn)
        lay = QVBoxLayout(self)
        lay.addWidget(scrollable(inner))
        lay.addWidget(self.build_log())

    def run(self):
        params = {
            "expr_counts_tsv": self.expr.text().strip(),
            "traits_csv": self.traits.text().strip(),
            "gene_id_col": self.gene_id_col.text().strip(),
            "trait_sample_col": self.trait_sample_col.text().strip(),
            "treatment_col": self.treatment_col.text().strip(),
            "out_dir": self.out_dir.text().strip(),
            "networkType": self.network_type.currentText(),
            "corFnc": self.corfnc.currentText(),
            "softPower": self.soft_power.value(),
            "minModuleSize": self.min_module.value(),
            "deepSplit": self.deep_split.value(),
            "MEDissThres": self.medissthres.value(),
            "cytoscape_edge_threshold": self.cyt_thr.value(),
            "sample_tree_cut_height": self.cut_height.value(),
            "n_threads": self.threads.text().strip(),
            "rscript_path": self.rscript.text().strip(),
        }
        self.run_async(pl.run_wgcna, params)


class DETab(AnalysisTab):
    def __init__(self, main):
        super().__init__(main)
        inner = QWidget()
        form = QFormLayout(inner)
        self.counts = QLineEdit()
        form.addRow("Count matrix TSV:", browse_file_row(self.counts, self,
                    filt="TSV (*.tsv *.txt);;All (*.*)"))
        self.meta = QLineEdit()
        form.addRow("Sample metadata CSV:", browse_file_row(self.meta, self,
                    filt="CSV (*.csv);;All (*.*)"))
        self.gene_id_col = QLineEdit("Geneid")
        form.addRow("Gene ID column (counts):", self.gene_id_col)
        self.sample_col = QLineEdit("sample")
        form.addRow("Sample-ID column (metadata):", self.sample_col)
        self.condition_col = QLineEdit("condition")
        form.addRow("Condition column (metadata):", self.condition_col)
        self.comparison_col = QLineEdit()
        self.comparison_col.setPlaceholderText("blank = single run; e.g. 'comparison' for per-group runs")
        self.comparison_col.setToolTip(
            "Optional. Name of a metadata column that groups samples into separate "
            "comparisons (e.g. values A, B, C, D). DESeq2 runs once per group, "
            "writing one DEG file per value, suffixed with the group name."
        )
        form.addRow("Comparison column (optional):", self.comparison_col)
        self.ref_level = QLineEdit("control")
        form.addRow("Reference level:", self.ref_level)
        self.test_level = QLineEdit("treated")
        form.addRow("Test level:", self.test_level)
        self.min_count = QSpinBox()
        self.min_count.setRange(0, 1000000)
        self.min_count.setValue(10)
        form.addRow("Min total count per gene:", self.min_count)
        self.output = QLineEdit(os.path.join(os.getcwd(), "DEG_result.tsv"))
        form.addRow("Output TSV:", browse_file_row(self.output, self,
                    save=True, filt="TSV (*.tsv)"))
        self.rscript = QLineEdit(pl.find_rscript())
        form.addRow("Rscript executable:", browse_file_row(self.rscript, self))
        self.run_btn = QPushButton("Run DESeq2")
        self.run_btn.clicked.connect(self.run)
        form.addRow(self.run_btn)
        lay = QVBoxLayout(self)
        lay.addWidget(scrollable(inner))
        lay.addWidget(self.build_log())

    def run(self):
        params = {
            "counts_tsv": self.counts.text().strip(),
            "meta_csv": self.meta.text().strip(),
            "gene_id_col": self.gene_id_col.text().strip(),
            "sample_col": self.sample_col.text().strip(),
            "condition_col": self.condition_col.text().strip(),
            "comparison_col": self.comparison_col.text().strip(),
            "ref_level": self.ref_level.text().strip(),
            "test_level": self.test_level.text().strip(),
            "min_count": self.min_count.value(),
            "output_tsv": self.output.text().strip(),
            "rscript_path": self.rscript.text().strip(),
        }
        self.run_async(pl.run_deseq2, params)


class FiltersTab(AnalysisTab):
    def __init__(self, main):
        super().__init__(main)
        from PySide6.QtWidgets import QCheckBox
        self.QCheckBox = QCheckBox
        inner = QWidget()
        form = QFormLayout(inner)

        # ---- Input files ----
        # Every file is optional: an empty path skips the step that needs it.
        note = QLabel("Leave any file blank to skip that step. If the Prediction CSV "
                      "is blank, sRNA-target pairs are built from the edges file.")
        note.setWordWrap(True)
        form.addRow(note)
        self.pred = QLineEdit()
        self.pred.setPlaceholderText("optional - blank = build pairs from the edges file")
        self.pred.setToolTip(
            "Optional. The scored sRNA-target prediction table. Only 'sRNA' and "
            "'Target' columns are required; any program score column that is missing "
            "is simply skipped. Leave blank to build the pairs from the WGCNA edges "
            "file instead (energy/probability filtering is then skipped)."
        )
        form.addRow("Prediction CSV:", browse_file_row(self.pred, self,
                    filt="CSV (*.csv);;All (*.*)"))
        self.deg = QLineEdit()
        self.deg.setPlaceholderText("One or more DEG TSV paths (comma-separated) - used for consensus")
        self.deg.setToolTip(
            "Add multiple DEG files to build a consensus. Browse... lets you select "
            "several at once (Ctrl/Shift-click) and adds them to any already listed. "
            "A gene is called consistent when at least 'Min strains consistent' files agree."
        )
        form.addRow("DEG TSV file(s):", browse_file_row(self.deg, self,
                    filt="TSV (*.tsv *.txt);;All (*.*)", multi=True))
        self.nodes = QLineEdit()
        form.addRow("WGCNA nodes file:", browse_file_row(self.nodes, self))
        self.edges = QLineEdit()
        form.addRow("WGCNA edges file:", browse_file_row(self.edges, self))
        self.srna_prefix = QLineEdit("sRNA")
        self.srna_prefix.setToolTip(
            "Used only when the Prediction CSV is blank and pairs are built from the "
            "edges file: a node whose name starts with this prefix is treated as the "
            "sRNA, the other endpoint as the target."
        )
        form.addRow("sRNA ID prefix (edges mode):", self.srna_prefix)

        # ---- Energy / probability filters (each toggleable) ----
        ebox = QGroupBox("Energy / probability filters (untick to disable; "
                         "missing columns are skipped automatically)")
        ef = QFormLayout(ebox)
        self.e_inta = self._dspin(-2.44, -1e6, 1e6, 2); self.e_inta_on = QCheckBox(); self.e_inta_on.setChecked(True)
        ef.addRow("E_intaRNA max (<=):", self._with_toggle(self.e_inta, self.e_inta_on))
        self.e_rnaplex = self._dspin(-32.6, -1e6, 1e6, 2); self.e_rnaplex_on = QCheckBox(); self.e_rnaplex_on.setChecked(True)
        ef.addRow("E_Rnaplex max (<=):", self._with_toggle(self.e_rnaplex, self.e_rnaplex_on))
        self.e_target = self._dspin(-5.13, -1e6, 1e6, 2); self.e_target_on = QCheckBox(); self.e_target_on.setChecked(True)
        ef.addRow("E_TargetRNA3 max (<=):", self._with_toggle(self.e_target, self.e_target_on))
        self.p_target = self._dspin(0.06, 0, 1, 3); self.p_target_on = QCheckBox(); self.p_target_on.setChecked(True)
        ef.addRow("Probability_TargetRNA3 min (>=):", self._with_toggle(self.p_target, self.p_target_on))
        self.p_rf = self._dspin(0.40, 0, 1, 3); self.p_rf_on = QCheckBox(); self.p_rf_on.setChecked(True)
        ef.addRow("Probability_sRNARFTarget min (>=):", self._with_toggle(self.p_rf, self.p_rf_on))
        self.energy_min_pass = QSpinBox(); self.energy_min_pass.setRange(0, 5); self.energy_min_pass.setValue(0)
        self.energy_min_pass.setSpecialValueText("All enabled (AND)")
        ef.addRow("Min metrics that must pass:", self.energy_min_pass)
        form.addRow(ebox)

        # ---- DEG consistency + direction ----
        dbox = QGroupBox("Differential expression (DEG) filters")
        df2 = QFormLayout(dbox)
        self.padj = QDoubleSpinBox(); self.padj.setRange(0.0, 1.0); self.padj.setDecimals(4)
        self.padj.setSingleStep(0.01); self.padj.setValue(0.05)
        df2.addRow("padj cutoff:", self.padj)
        self.min_strains = QSpinBox(); self.min_strains.setRange(1, 1000); self.min_strains.setValue(4)
        df2.addRow("Min strains consistent:", self.min_strains)
        self.deg_min_basemean = QDoubleSpinBox(); self.deg_min_basemean.setRange(0.0, 1e9)
        self.deg_min_basemean.setDecimals(1); self.deg_min_basemean.setSingleStep(10.0)
        self.deg_min_basemean.setValue(0.0)
        self.deg_min_basemean.setToolTip(
            "Minimum DESeq2 baseMean (mean normalized count) required for a gene to "
            "count as a DEG. 0 = off. Files without a baseMean column are unaffected."
        )
        df2.addRow("Min baseMean (0 = off):", self.deg_min_basemean)
        self.req_srna = self.QCheckBox("sRNA must be a DEG"); self.req_srna.setChecked(True)
        df2.addRow(self.req_srna)
        self.req_target = self.QCheckBox("Target must be a DEG"); self.req_target.setChecked(True)
        df2.addRow(self.req_target)
        self.srna_dir = QComboBox(); self.srna_dir.addItems(["any", "upregulated", "downregulated"])
        df2.addRow("Required sRNA direction:", self.srna_dir)
        self.target_dir = QComboBox(); self.target_dir.addItems(["any", "upregulated", "downregulated"])
        df2.addRow("Required target direction:", self.target_dir)
        self.pair_rel = QComboBox(); self.pair_rel.addItems(["any", "same", "opposite"])
        df2.addRow("sRNA vs target direction:", self.pair_rel)
        form.addRow(dbox)

        # ---- WGCNA module + network weight ----
        mbox = QGroupBox("WGCNA module & network filters")
        mf = QFormLayout(mbox)
        self.module_mode = QComboBox(); self.module_mode.addItems(["same", "different", "any"])
        mf.addRow("sRNA / target module:", self.module_mode)
        self.require_net = self.QCheckBox("Interaction must have a network edge"); self.require_net.setChecked(True)
        mf.addRow(self.require_net)
        self.min_weight = QDoubleSpinBox(); self.min_weight.setRange(0.0, 1.0); self.min_weight.setDecimals(4)
        self.min_weight.setSingleStep(0.01); self.min_weight.setValue(0.0)
        mf.addRow("Min edge weight (0 = off):", self.min_weight)
        self.selection = QComboBox(); self.selection.addItems(["best_per_target", "best_per_srna", "all"])
        mf.addRow("Keep:", self.selection)
        form.addRow(mbox)

        self.output = QLineEdit(os.path.join(os.getcwd(), "filtered_weight.csv"))
        form.addRow("Output CSV:", browse_file_row(self.output, self,
                    save=True, filt="CSV (*.csv)"))
        self.run_btn = QPushButton("Run Filters")
        self.run_btn.clicked.connect(self.run)
        form.addRow(self.run_btn)

        lay = QVBoxLayout(self)
        lay.addWidget(scrollable(inner))
        lay.addWidget(self.build_log())

    @staticmethod
    def _dspin(val, lo, hi, decimals):
        s = QDoubleSpinBox(); s.setDecimals(decimals); s.setRange(lo, hi); s.setValue(val)
        return s

    @staticmethod
    def _with_toggle(spin, check):
        row = QHBoxLayout(); row.addWidget(check); row.addWidget(spin)
        w = QWidget(); w.setLayout(row)
        return w

    def run(self):
        deg_files = [x.strip() for x in self.deg.text().split(",") if x.strip()]
        params = {
            "predictions_csv": self.pred.text().strip(),
            "deg_files": deg_files,
            "module_nodes_file": self.nodes.text().strip(),
            "module_edges_file": self.edges.text().strip(),
            "srna_prefix": self.srna_prefix.text().strip(),
            "filters_outlier": {
                "E_intaRNA_max": self.e_inta.value(),
                "E_Rnaplex_max": self.e_rnaplex.value(),
                "E_TargetRNA3_max": self.e_target.value(),
                "Probability_TargetRNA3_min": self.p_target.value(),
                "Probability_sRNARFTarget_min": self.p_rf.value(),
            },
            "energy_enabled": {
                "E_intaRNA": self.e_inta_on.isChecked(),
                "E_Rnaplex": self.e_rnaplex_on.isChecked(),
                "E_TargetRNA3": self.e_target_on.isChecked(),
                "Probability_TargetRNA3": self.p_target_on.isChecked(),
                "Probability_sRNARFTarget": self.p_rf_on.isChecked(),
            },
            "energy_min_pass": self.energy_min_pass.value(),
            "padj_cutoff": self.padj.value(),
            "min_strains_consistent": self.min_strains.value(),
            "deg_min_basemean": self.deg_min_basemean.value(),
            "require_srna_deg": self.req_srna.isChecked(),
            "require_target_deg": self.req_target.isChecked(),
            "srna_direction": self.srna_dir.currentText(),
            "target_direction": self.target_dir.currentText(),
            "pair_relationship": self.pair_rel.currentText(),
            "module_mode": self.module_mode.currentText(),
            "require_in_network": self.require_net.isChecked(),
            "min_edge_weight": self.min_weight.value(),
            "selection_mode": self.selection.currentText(),
            "output_csv": self.output.text().strip(),
        }
        self.run_async(pl.run_filters, params)


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("sRNA Pipeline - Integrative RNA-Seq & Co-expression")
        self.resize(860, 820)
        tabs = QTabWidget()
        self.position_tab = PositionTab(self)
        self.wgcna_tab = WGCNATab(self)
        self.de_tab = DETab(self)
        self.filters_tab = FiltersTab(self)
        tabs.addTab(self.position_tab, "01 Position Classification")
        tabs.addTab(self.wgcna_tab, "02 Co-expression analysis")
        tabs.addTab(self.de_tab, "03 Differential Expression")
        tabs.addTab(self.filters_tab, "04 Filters")
        self.setCentralWidget(tabs)


def main():
    app = QApplication(sys.argv)
    win = MainWindow()
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
