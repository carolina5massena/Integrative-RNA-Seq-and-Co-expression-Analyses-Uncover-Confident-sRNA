"""
sRNA Pipeline GUI - PySide6 desktop app for the Integrative RNA-Seq workflow.

Tabs: 01 Position Classification, 02 Co-expression (WGCNA),
03 Differential Expression (DESeq2), 04 Combine Predictions, 05 Filters.

Filter order (05 Filters tab): Energy -> DEG -> WGCNA -> Position (always last).

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


def results_dir():
    """Default output location: a 'Results' folder under the current directory.

    Created on demand so every tab's default output path points inside it.
    """
    d = os.path.join(os.getcwd(), "Results")
    os.makedirs(d, exist_ok=True)
    return d


class Worker(QObject):
    log = Signal(str)
    finished = Signal(bool, str)

    def __init__(self, func, params):
        super().__init__()
        self.func = func
        self.params = params

    def run(self):
        # Catch everything so a bad input or unexpected data never takes the app
        # down: the full traceback goes to the log and the run ends cleanly with a
        # failure status instead of crashing.
        try:
            out = self.func(self.params, lambda m: self.log.emit(str(m)))
            self.finished.emit(True, f"Done. Output: {out}")
        except BaseException as e:
            self.log.emit(
                "\n" + "=" * 64 +
                "\nERROR - the run stopped. Nothing was crashed; details below.\n" +
                "=" * 64 + "\n" + traceback.format_exc()
            )
            self.finished.emit(False, f"{type(e).__name__}: {e}")


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
        try:
            self.log_view.appendPlainText(msg)
        except Exception:
            pass  # never let a logging hiccup take the app down

    def run_async(self, func, params):
        if self.thread is not None and self.thread.isRunning():
            QMessageBox.warning(self, "Busy", "A task is already running.")
            return
        # Make sure any thread from a previous run has fully stopped and been
        # released before we start a new one. Reusing/recreating a QThread
        # without this is the usual cause of the app closing unexpectedly after
        # several runs ("QThread: Destroyed while thread is still running").
        self._teardown_thread()
        self.log_view.clear()
        if self.run_btn:
            self.run_btn.setEnabled(False)

        self.thread = QThread()
        self.worker = Worker(func, params)
        self.worker.moveToThread(self.thread)
        self.thread.started.connect(self.worker.run)
        self.worker.log.connect(self.append_log)
        self.worker.finished.connect(self._on_finished)
        # Ordered, safe teardown. Stop the thread's event loop when the work is
        # done, then let Qt delete the worker and the thread ONLY after the loop
        # has fully stopped (on thread.finished). We must NOT delete the worker on
        # its own `finished` signal: _on_finished shows a modal dialog that spins
        # the event loop, and freeing the signal's sender mid-emission would
        # hard-crash the app. Deferring deletion to thread.finished avoids that.
        self.worker.finished.connect(self.thread.quit)
        self.thread.finished.connect(self.worker.deleteLater)
        self.thread.finished.connect(self.thread.deleteLater)
        self.thread.finished.connect(self._on_thread_finished)
        self.thread.start()

    def _on_thread_finished(self):
        # Drop our Python references only after Qt has finished with the objects.
        self.worker = None
        self.thread = None

    def _teardown_thread(self):
        """Stop and wait for a lingering thread from a previous run, if any."""
        t = self.thread
        if t is not None:
            try:
                if t.isRunning():
                    t.quit()
                    t.wait()
            except RuntimeError:
                pass  # underlying C++ object was already deleted
        self.thread = None
        self.worker = None

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
        self.output_file = QLineEdit(os.path.join(results_dir(), "sRNA_annotation.xlsx"))
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
        self.out_dir = QLineEdit(results_dir())
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
        self.output = QLineEdit(os.path.join(results_dir(), "DEG_result.tsv"))
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


class CombineTab(AnalysisTab):
    """Optional: merge the outputs of the four sRNA-target prediction programs
    into one Prediction table (same columns as Prediction_test.csv)."""

    def __init__(self, main):
        super().__init__(main)
        from PySide6.QtWidgets import QCheckBox
        inner = QWidget()
        form = QFormLayout(inner)

        note = QLabel(
            "Combine the native output of up to four prediction programs into a "
            "single table like Prediction_test.csv. Every file is optional - leave "
            "a program blank to skip it. Files are merged on each (Target, sRNA) "
            "pair. intaRNA and RNAplex report only an energy, so their p-value "
            "column is computed as an empirical p-value from the energy distribution."
        )
        note.setWordWrap(True)
        form.addRow(note)

        self.intarna = QLineEdit()
        self.intarna.setPlaceholderText("intaRNA output (';'-separated: id1;id2;...;E)")
        form.addRow("intaRNA file:", browse_file_row(self.intarna, self,
                    filt="intaRNA (*.txt *.csv);;All (*.*)"))
        self.rnaplex = QLineEdit()
        self.rnaplex.setPlaceholderText("RNAplex output (>target / >sRNA / structure (energy))")
        form.addRow("RNAplex file:", browse_file_row(self.rnaplex, self,
                    filt="RNAplex (*.txt);;All (*.*)"))
        self.targetrna3 = QLineEdit()
        self.targetrna3.setPlaceholderText("TargetRNA3 output (tab; needs an sRNA column)")
        self.targetrna3.setToolTip(
            "Tab-separated TargetRNA3 table. Because TargetRNA3 runs one sRNA at a "
            "time, the file must include an 'sRNA' column identifying the sRNA, plus "
            "'Target', 'Energy', 'P-value' and 'Probability'.")
        form.addRow("TargetRNA3 file:", browse_file_row(self.targetrna3, self,
                    filt="TargetRNA3 (*.txt *.tsv);;All (*.*)"))
        self.srnarf = QLineEdit()
        self.srnarf.setPlaceholderText("sRNARFTarget output (sRNA_ID, mRNA_ID, Prediction_Probability)")
        form.addRow("sRNARFTarget file:", browse_file_row(self.srnarf, self,
                    filt="sRNARFTarget (*.csv *.tsv *.txt);;All (*.*)"))

        obox = QGroupBox("Options")
        of = QFormLayout(obox)
        self.empirical_p = QCheckBox("Compute empirical p-values for intaRNA / RNAplex")
        self.empirical_p.setChecked(True)
        of.addRow(self.empirical_p)
        self.require_all = QCheckBox("Keep only pairs predicted by every provided program")
        self.require_all.setChecked(True)
        self.require_all.setToolTip(
            "Inner join: a (Target, sRNA) pair is kept only if it appears in all the "
            "files you provided. Untick to keep every pair (outer join), leaving "
            "blanks where a program has no value.")
        of.addRow(self.require_all)
        self.relabel = QCheckBox("Relabel IDs to Gene# / sRNA#")
        self.relabel.setChecked(False)
        self.relabel.setToolTip(
            "Replace the original Target/sRNA names with tidy Gene1, Gene2... / "
            "sRNA1, sRNA2... labels and save a mapping file next to the output.")
        of.addRow(self.relabel)
        form.addRow(obox)

        self.output = QLineEdit(os.path.join(results_dir(), "Prediction_combined.csv"))
        form.addRow("Output CSV:", browse_file_row(self.output, self,
                    save=True, filt="CSV (*.csv)"))
        self.run_btn = QPushButton("Combine Predictions")
        self.run_btn.clicked.connect(self.run)
        form.addRow(self.run_btn)

        lay = QVBoxLayout(self)
        lay.addWidget(scrollable(inner))
        lay.addWidget(self.build_log())

    def run(self):
        out = self.output.text().strip()
        mapping = ""
        if self.relabel.isChecked() and out:
            base, _ = os.path.splitext(out)
            mapping = base + "_id_mapping.csv"
        params = {
            "intarna_file": self.intarna.text().strip(),
            "rnaplex_file": self.rnaplex.text().strip(),
            "targetrna3_file": self.targetrna3.text().strip(),
            "srnarftarget_file": self.srnarf.text().strip(),
            "compute_empirical_p": self.empirical_p.isChecked(),
            "require_all_programs": self.require_all.isChecked(),
            "relabel_ids": self.relabel.isChecked(),
            "mapping_csv": mapping,
            "output_csv": out,
        }
        self.run_async(pl.run_combine_predictions, params)


class FiltersTab(AnalysisTab):
    def __init__(self, main):
        super().__init__(main)
        from PySide6.QtWidgets import QCheckBox
        self.QCheckBox = QCheckBox
        inner = QWidget()
        form = QFormLayout(inner)

        # ---- Input files ----
        # Each filter is switched on simply by providing its file here. When a file
        # is left blank the matching filter is skipped AND its detail options stay
        # hidden below (no separate enable checkbox is needed).
        note = QLabel(
            "Each filter is turned on by giving it a file below. Leave a file blank "
            "to skip that analysis - its options stay hidden. Pairs come from the "
            "Prediction CSV if given, otherwise from the WGCNA edges file. The "
            "position filter always runs last."
        )
        note.setWordWrap(True)
        form.addRow(note)
        self.pred = QLineEdit()
        self.pred.setPlaceholderText("blank = build pairs from the edges file; given = enables energy filters")
        self.pred.setToolTip(
            "The scored sRNA-target prediction table. Only 'sRNA' and 'Target' "
            "columns are required; any program score column that is missing is "
            "simply skipped. Providing it shows the Energy / probability filters. "
            "Leave blank to build the pairs from the WGCNA edges file instead."
        )
        form.addRow("Prediction CSV:", browse_file_row(self.pred, self,
                    filt="CSV (*.csv);;All (*.*)"))
        self.deg = QLineEdit()
        self.deg.setPlaceholderText("blank = no DEG filter; one or more TSV paths (comma-separated) = enables it")
        self.deg.setToolTip(
            "Add one or more DEG files to switch on the DEG filter and its consensus. "
            "Browse... lets you select several at once (Ctrl/Shift-click) and adds "
            "them to any already listed. Leave blank to skip DEG filtering entirely."
        )
        form.addRow("DEG TSV file(s):", browse_file_row(self.deg, self,
                    filt="TSV (*.tsv *.txt);;All (*.*)", multi=True))
        self.nodes = QLineEdit()
        self.nodes.setPlaceholderText("blank = no module filter; given = enables the module-colour filter")
        self.nodes.setToolTip(
            "WGCNA CytoscapeInput-nodes file. Providing it shows the module-colour "
            "filter (keep pairs in the same or different modules)."
        )
        form.addRow("WGCNA nodes file:", browse_file_row(self.nodes, self))
        self.edges = QLineEdit()
        self.edges.setPlaceholderText("blank = no network filter; given = enables edge/weight/best-pair filters")
        self.edges.setToolTip(
            "WGCNA CytoscapeInput-edges file. Providing it shows the network filters "
            "(require an edge, minimum weight, keep best pair per group). It is also "
            "used as the pair source when no Prediction CSV is given."
        )
        form.addRow("WGCNA edges file:", browse_file_row(self.edges, self))
        self.loc_file = QLineEdit()
        self.loc_file.setPlaceholderText("blank = no position filter; given = enables it (runs LAST)")
        self.loc_file.setToolTip(
            "The sRNA_annotation table from the 01 Position Classification tab "
            "(.xlsx/.csv). Providing it switches on the position filter, which always "
            "runs last. Must contain the sRNA ID column and a 'location' column."
        )
        form.addRow("sRNA annotation (position):", browse_file_row(self.loc_file, self,
                    filt="Annotation (*.xlsx *.csv);;Excel (*.xlsx);;CSV (*.csv);;All (*.*)"))
        self.srna_prefix = QLineEdit("sRNA")
        self.srna_prefix.setToolTip(
            "Used only when the Prediction CSV is blank and pairs are built from the "
            "edges file: a node whose name starts with this prefix is treated as the "
            "sRNA, the other endpoint as the target."
        )
        form.addRow("sRNA ID prefix (edges mode):", self.srna_prefix)

        # ---- Energy / probability filters (shown only with a Prediction CSV) ----
        self.ebox = QGroupBox("Energy / probability filters (untick a metric to skip it; "
                              "missing columns are skipped automatically)")
        ef = QFormLayout(self.ebox)
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
        self.pv_inta = self._dspin(0.41, 0, 1, 3); self.pv_inta_on = QCheckBox(); self.pv_inta_on.setChecked(True)
        ef.addRow("p_intaRNA max (<=):", self._with_toggle(self.pv_inta, self.pv_inta_on))
        self.pv_rnaplex = self._dspin(0.76, 0, 1, 3); self.pv_rnaplex_on = QCheckBox(); self.pv_rnaplex_on.setChecked(True)
        ef.addRow("p_Rnaplex max (<=):", self._with_toggle(self.pv_rnaplex, self.pv_rnaplex_on))
        self.pv_target = self._dspin(0.89, 0, 1, 3); self.pv_target_on = QCheckBox(); self.pv_target_on.setChecked(True)
        ef.addRow("p_TargetRNA3 max (<=):", self._with_toggle(self.pv_target, self.pv_target_on))
        self.energy_min_pass = QSpinBox(); self.energy_min_pass.setRange(0, 8); self.energy_min_pass.setValue(0)
        self.energy_min_pass.setSpecialValueText("All enabled (AND)")
        ef.addRow("Min metrics that must pass:", self.energy_min_pass)
        form.addRow(self.ebox)

        # ---- DEG consistency + direction (shown only with DEG file(s)) ----
        self.dbox = QGroupBox("Differential expression (DEG) filters "
                              "(shown when DEG file(s) are provided above)")
        df2 = QFormLayout(self.dbox)
        self.padj = QDoubleSpinBox(); self.padj.setRange(0.0, 1.0); self.padj.setDecimals(4)
        self.padj.setSingleStep(0.01); self.padj.setValue(0.05)
        df2.addRow("padj cutoff (<=):", self.padj)
        # baseMean, with its own enable toggle
        self.deg_basemean_on = QCheckBox()
        self.deg_min_basemean = QDoubleSpinBox(); self.deg_min_basemean.setRange(0.0, 1e9)
        self.deg_min_basemean.setDecimals(1); self.deg_min_basemean.setSingleStep(10.0)
        self.deg_min_basemean.setValue(0.0)
        self.deg_min_basemean.setToolTip(
            "Tick to enable. Minimum DESeq2 baseMean (mean normalized count) required "
            "for a gene to count as a DEG within a file. Files without a baseMean "
            "column are unaffected."
        )
        df2.addRow("Min baseMean (>=):", self._with_toggle(self.deg_min_basemean, self.deg_basemean_on))
        # log2FC inclusion band, applied per DEG file. A gene is called when its
        # log2FC is >= the up threshold OR <= the down threshold (the band between
        # them is discarded). Untick a side to disable it; untick both to keep any.
        self.lfc_up_on = QCheckBox(); self.lfc_up_on.setChecked(True)
        self.lfc_up = self._dspin(1.0, -1e6, 1e6, 2); self.lfc_up.setSingleStep(0.5)
        df2.addRow("Keep log2FC >= :", self._with_toggle(self.lfc_up, self.lfc_up_on))
        self.lfc_down_on = QCheckBox(); self.lfc_down_on.setChecked(True)
        self.lfc_down = self._dspin(-1.0, -1e6, 1e6, 2); self.lfc_down.setSingleStep(0.5)
        df2.addRow("Keep log2FC <= :", self._with_toggle(self.lfc_down, self.lfc_down_on))
        for w in (self.lfc_up, self.lfc_down):
            w.setToolTip(
                "log2FC inclusion band, applied per DEG file when calling a gene a "
                "DEG. The gene passes when log2FC >= the up value OR <= the down "
                "value (the band between them is discarded). Untick a side to disable "
                "it; untick both to keep any log2FC."
            )
        # Consensus: type (sign vs presence) + how many files (X)
        self.consensus_mode = QComboBox()
        self.consensus_mode.addItems(["sign of log2FC", "presence (ignore log2FC)"])
        self.consensus_mode.setToolTip(
            "How consensus across the DEG files is decided.\n"
            "'sign of log2FC': the gene must be up (log2FC>0) in >= X files, or down "
            "(log2FC<0) in >= X files.\n"
            "'presence (ignore log2FC)': the gene only has to be called (pass padj / "
            "baseMean / log2FC band) in >= X files, regardless of sign."
        )
        df2.addRow("Consensus type:", self.consensus_mode)
        self.min_strains = QSpinBox(); self.min_strains.setRange(1, 1000); self.min_strains.setValue(2)
        self.min_strains.setToolTip(
            "X: the minimum number of DEG files a gene must be called in to reach "
            "consensus. A gene absent from a file is treated as failing that file."
        )
        df2.addRow("Consensus: min files (X):", self.min_strains)
        # Direction / relationship sub-filters, each with an enable toggle.
        self.srna_dir_on = QCheckBox()
        self.srna_dir = QComboBox(); self.srna_dir.addItems(["upregulated", "downregulated"])
        df2.addRow("Required sRNA direction:", self._with_toggle(self.srna_dir, self.srna_dir_on))
        self.target_dir_on = QCheckBox()
        self.target_dir = QComboBox(); self.target_dir.addItems(["upregulated", "downregulated"])
        df2.addRow("Required target direction:", self._with_toggle(self.target_dir, self.target_dir_on))
        self.pair_rel_on = QCheckBox()
        self.pair_rel = QComboBox(); self.pair_rel.addItems(["same", "opposite"])
        df2.addRow("sRNA vs target direction:", self._with_toggle(self.pair_rel, self.pair_rel_on))
        form.addRow(self.dbox)

        # ---- WGCNA module filter (shown only with a nodes file) ----
        self.module_box = QGroupBox("WGCNA module filter "
                                    "(shown when a nodes file is provided above)")
        modf = QFormLayout(self.module_box)
        self.module_mode = QComboBox(); self.module_mode.addItems(["same", "different"])
        self.module_mode.setToolTip(
            "Keep only pairs whose sRNA and target modules are the same colour, or "
            "the different colours. Providing the nodes file switches this filter on; "
            "leave it blank to keep any module pairing."
        )
        modf.addRow("Keep module colours:", self.module_mode)
        form.addRow(self.module_box)

        # ---- WGCNA network filter (shown only with an edges file) ----
        self.network_box = QGroupBox("WGCNA network filter "
                                     "(shown when an edges file is provided above)")
        netf = QFormLayout(self.network_box)
        self.require_net = self.QCheckBox("Interaction must have a network edge"); self.require_net.setChecked(True)
        netf.addRow(self.require_net)
        self.min_weight_on = self.QCheckBox()
        self.min_weight = QDoubleSpinBox(); self.min_weight.setRange(0.0, 1.0); self.min_weight.setDecimals(4)
        self.min_weight.setSingleStep(0.01); self.min_weight.setValue(0.0)
        self.min_weight.setToolTip("Tick to enable. Keep only pairs whose network "
                                   "edge weight is >= this value.")
        netf.addRow("Min edge weight (>=):", self._with_toggle(self.min_weight, self.min_weight_on))
        self.selection_on = self.QCheckBox()
        self.selection = QComboBox(); self.selection.addItems(["best_per_target", "best_per_srna"])
        self.selection.setToolTip(
            "Tick to enable. Keep only the best pair per target (unique targets) or "
            "per sRNA (unique sRNAs), ranked by edge weight. Untick to keep every "
            "pair that survived the earlier filters."
        )
        netf.addRow("Keep best per:", self._with_toggle(self.selection, self.selection_on))
        self.top_n = QSpinBox(); self.top_n.setRange(0, 1000000); self.top_n.setValue(0)
        self.top_n.setSpecialValueText("Best only (1)")
        self.top_n.setToolTip(
            "With 'Keep best per' enabled, keep the top N partners per group, ranked "
            "by network edge weight (highest first). best_per_srna + Top N = 10 keeps "
            "the top 10 targets for each sRNA; best_per_target + Top N = 10 keeps the "
            "top 10 sRNAs for each target. 0 = keep only the single best."
        )
        netf.addRow("Top N per group:", self.top_n)
        form.addRow(self.network_box)

        # ---- Position filter (ALWAYS LAST; shown only with an annotation file) ----
        self.lbox = QGroupBox("Position filter (runs LAST; shown when the sRNA "
                              "annotation file is provided above)")
        lf = QFormLayout(self.lbox)
        self.loc_id_col = QLineEdit("exon_id")
        self.loc_id_col.setToolTip(
            "Column in the annotation file holding the sRNA name that matches the "
            "'sRNA' column of the pairs (default 'exon_id')."
        )
        lf.addRow("sRNA ID column:", self.loc_id_col)
        self.loc_checks = {}
        loc_row = QHBoxLayout()
        for name in ["Intragenic", "Antisense", "5' UTR", "3' UTR", "Intergenic"]:
            cb = self.QCheckBox(name); cb.setChecked(True)
            self.loc_checks[name] = cb
            loc_row.addWidget(cb)
        loc_w = QWidget(); loc_w.setLayout(loc_row)
        lf.addRow("Keep locations:", loc_w)
        form.addRow(self.lbox)

        # ---- File-driven visibility -------------------------------------------
        # A filter's option box is shown only when its input file is provided.
        for le in (self.pred, self.deg, self.nodes, self.edges, self.loc_file):
            le.textChanged.connect(self._update_filter_visibility)

        self.output = QLineEdit(os.path.join(results_dir(), "filtered_weight.csv"))
        form.addRow("Output CSV:", browse_file_row(self.output, self,
                    save=True, filt="CSV (*.csv)"))
        # When ticked, a separate table is written after every filtering step
        # (base pairs, energy, DEG, WGCNA, position), named after the output CSV
        # with a step suffix, so each step's result can be inspected on its own.
        self.save_intermediate = self.QCheckBox("Generate all intermediate tables")
        self.save_intermediate.setToolTip(
            "Write one table per filtering step. After each step the current "
            "result is saved next to the Output CSV, using the output name plus a "
            "step suffix: _step1_base, _step2_energy, _step3_deg, _step4_wgcna "
            "(same-module + network-edge requirements), _step5_bestpair (after "
            "best-pair-per-group prioritization) and _step6_position. A step that "
            "is skipped still produces a table (identical to the previous step) so "
            "the full funnel is captured."
        )
        form.addRow(self.save_intermediate)
        self.run_btn = QPushButton("Run Filters")
        self.run_btn.clicked.connect(self.run)
        form.addRow(self.run_btn)

        lay = QVBoxLayout(self)
        lay.addWidget(scrollable(inner))
        lay.addWidget(self.build_log())

        # Set the initial show/hide state to match the (empty) file fields.
        self._update_filter_visibility()

    def _update_filter_visibility(self):
        """Show each filter's option box only when its input file is provided."""
        self.ebox.setVisible(bool(self.pred.text().strip()))
        self.dbox.setVisible(bool(self.deg.text().strip()))
        self.module_box.setVisible(bool(self.nodes.text().strip()))
        self.network_box.setVisible(bool(self.edges.text().strip()))
        self.lbox.setVisible(bool(self.loc_file.text().strip()))

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
        # A filter is enabled purely by the presence of its input file.
        has_pred = bool(self.pred.text().strip())
        has_deg = bool(deg_files)
        has_nodes = bool(self.nodes.text().strip())
        has_edges = bool(self.edges.text().strip())
        has_loc = bool(self.loc_file.text().strip())
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
                "p_intaRNA_max": self.pv_inta.value(),
                "p_Rnaplex_max": self.pv_rnaplex.value(),
                "p_TargetRNA3_max": self.pv_target.value(),
            },
            "energy_enabled": {
                "E_intaRNA": self.e_inta_on.isChecked(),
                "E_Rnaplex": self.e_rnaplex_on.isChecked(),
                "E_TargetRNA3": self.e_target_on.isChecked(),
                "Probability_TargetRNA3": self.p_target_on.isChecked(),
                "Probability_sRNARFTarget": self.p_rf_on.isChecked(),
                "p_intaRNA": self.pv_inta_on.isChecked(),
                "p_Rnaplex": self.pv_rnaplex_on.isChecked(),
                "p_TargetRNA3": self.pv_target_on.isChecked(),
            },
            "energy_min_pass": self.energy_min_pass.value(),
            "deg_enabled": has_deg,
            "padj_cutoff": self.padj.value(),
            "min_strains_consistent": self.min_strains.value(),
            "consensus_mode": ("presence" if self.consensus_mode.currentIndex() == 1
                               else "sign"),
            "deg_basemean_enabled": self.deg_basemean_on.isChecked(),
            "deg_min_basemean": self.deg_min_basemean.value(),
            "lfc_up_enabled": self.lfc_up_on.isChecked(),
            "lfc_up": self.lfc_up.value(),
            "lfc_down_enabled": self.lfc_down_on.isChecked(),
            "lfc_down": self.lfc_down.value(),
            "srna_dir_enabled": self.srna_dir_on.isChecked(),
            "srna_direction": self.srna_dir.currentText(),
            "target_dir_enabled": self.target_dir_on.isChecked(),
            "target_direction": self.target_dir.currentText(),
            "pair_rel_enabled": self.pair_rel_on.isChecked(),
            "pair_relationship": self.pair_rel.currentText(),
            "module_filter_enabled": has_nodes,
            "module_mode": self.module_mode.currentText(),
            "require_in_network": self.require_net.isChecked() and has_edges,
            "min_weight_enabled": self.min_weight_on.isChecked() and has_edges,
            "min_edge_weight": self.min_weight.value(),
            "selection_enabled": self.selection_on.isChecked() and has_edges,
            "selection_mode": self.selection.currentText(),
            "top_n_per_group": self.top_n.value(),
            "location_filter_enabled": has_loc,
            "srna_annotation_file": self.loc_file.text().strip(),
            "location_id_col": self.loc_id_col.text().strip() or "exon_id",
            "locations_keep": [name for name, cb in self.loc_checks.items() if cb.isChecked()],
            "output_csv": self.output.text().strip(),
            "save_intermediate": self.save_intermediate.isChecked(),
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
        self.combine_tab = CombineTab(self)
        self.filters_tab = FiltersTab(self)
        tabs.addTab(self.position_tab, "01 Position Classification")
        tabs.addTab(self.wgcna_tab, "02 Co-expression analysis")
        tabs.addTab(self.de_tab, "03 Differential Expression")
        tabs.addTab(self.combine_tab, "04 Combine Predictions")
        tabs.addTab(self.filters_tab, "05 Filters")
        self.setCentralWidget(tabs)


def _install_global_excepthook():
    """Report uncaught exceptions instead of letting Qt abort the process.

    Without this, an exception raised inside a Qt slot (e.g. a button click or a
    signal handler on the main thread) can terminate the whole application. Here
    we print the full traceback and show a dialog, then let the app keep running.
    (Errors inside a running analysis are already caught by Worker.run and sent
    to that tab's log.)
    """
    def hook(exc_type, exc, tb):
        if issubclass(exc_type, KeyboardInterrupt):
            sys.__excepthook__(exc_type, exc, tb)
            return
        sys.stderr.write("".join(traceback.format_exception(exc_type, exc, tb)))
        try:
            QMessageBox.critical(
                None, "Unexpected error",
                "An unexpected error occurred, but the application is still "
                f"running:\n\n{exc_type.__name__}: {exc}"
            )
        except Exception:
            pass
    sys.excepthook = hook


def main():
    app = QApplication(sys.argv)
    _install_global_excepthook()
    win = MainWindow()
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
