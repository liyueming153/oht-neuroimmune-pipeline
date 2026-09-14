#!/bin/bash
# Self-healing pipeline runner v2. Per-sample + slim checkpoints on /mnt.
M=/mnt/agents/work/pipe
W=/tmp/ohtrun
mkdir -p $W/out $M/h5ad/samples $M/h5ad/files
export GSE293358_DIR=/mnt/agents/work/gse293358
export OUT_DIR=$W/out
export CKPT_DIR=$M/h5ad/samples
export SKIP_MERGED=1
export PYTHONPATH=$W/pylibs
cd /mnt/agents/work/oht-neuroimmune-pipeline/scripts

echo "[runner] start $(date)"
if ! python3 -c "import scanpy,scrublet" 2>/dev/null; then
  echo "[runner] installing libs from wheelhouse"
  mkdir -p $W/pylibs
  pip install -q --no-index --no-deps --target=$W/pylibs /mnt/agents/work/wheelhouse/*.whl >> $M/out/runner.log 2>&1
  tar xzf /mnt/agents/work/wheelhouse/annoy_built.tgz -C $W/pylibs
fi
python3 -c "import scanpy,scrublet" || { echo "[runner] LIBS FAIL"; exit 1; }
echo "[runner] libs ok"

push_file(){ # $1=local $2=name  (single file <100MB)
  sz=$(stat -c%s "$1"); [ "$sz" -lt 100000000 ] || { echo "[runner] TOO BIG $2 $sz"; return 1; }
  cp "$1" "$M/h5ad/files/$2" && sync && echo "[runner] pushed $2 ($sz)"; }
pull_file(){ [ -f "$M/h5ad/files/$1" ] && cp "$M/h5ad/files/$1" "$2" && echo "[runner] pulled $1"; }
push_small(){ cp "$1" "$M/out/$(basename $1)" && sync; }

# 01: per-sample QC (resumable per sample)
n=$(ls $M/h5ad/samples/*.h5ad 2>/dev/null | wc -l)
if [ "$n" -ge 13 ]; then echo "[runner] 01 cached ($n samples)"; else
  python3 01_load_qc.py && cp $W/out/qc_table.csv $M/out/ && sync || echo "[runner] 01 incomplete"
fi
n=$(ls $M/h5ad/samples/*.h5ad 2>/dev/null | wc -l)
[ "$n" -ge 13 ] || { echo "[runner] 01 partial ($n/13), exit for next window"; exit 0; }
echo "[runner] 01 done"

# 02a: cluster (vstack from sample cache), slim ckpt
if pull_file clustered_slim.h5ad $W/out/clustered_slim.h5ad; then :; else
  python3 02a_cluster.py && push_file $W/out/clustered_slim.h5ad clustered_slim.h5ad || { echo "[runner] 02a FAIL/RETRY"; exit 0; }
fi
echo "[runner] 02a done $(date)"

# 02: annotate
if pull_file annotated.h5ad $W/out/annotated.h5ad; then :; else
  python3 02_annotate.py && push_file $W/out/annotated.h5ad annotated.h5ad && push_small $W/out/cell_class_zscores.csv.gz || { echo "[runner] 02 FAIL/RETRY"; exit 0; }
fi
echo "[runner] 02 done $(date)"

# 08: pseudotime
if pull_file microglia_dpt.h5ad $W/out/microglia_dpt.h5ad; then :; else
  python3 08_pseudotime.py && push_file $W/out/microglia_dpt.h5ad microglia_dpt.h5ad && push_small $W/out/dpt_trajectory_genes.csv || { echo "[runner] 08 FAIL/RETRY"; exit 0; }
fi
echo "[runner] 08 done $(date)"

# single-cell figures
if [ -f $M/out/fig3d_new.png ]; then echo "[runner] figs cached"; else
  cd /mnt/agents/work/pipe/figs && python3 fig_sc_all.py || { echo "[runner] figs FAIL/RETRY"; exit 0; }
  for f in fig1a_new.png fig1b_new.png fig1c_new.png fig1d_new.png fig2a_new.png fig2b_new.png fig2c_new.png fig3a_new.png fig3b_new.png fig3c_new.png fig3d_new.png fig5b_new.png composition_stats.csv mgstate_stats.csv dpt_group_tests.csv; do push_small $W/out/$f 2>/dev/null || push_small ./$f; done; cd /mnt/agents/work/oht-neuroimmune-pipeline/scripts
fi
echo "[runner] figs done $(date)"

# root sensitivity
if [ -f $M/out/root_sensitivity.csv ]; then echo "[runner] rootsens cached"; else
  cd /mnt/agents/work/pipe/figs && python3 root_sens.py && push_small root_sensitivity.csv; cd /mnt/agents/work/oht-neuroimmune-pipeline/scripts
fi
echo "[runner] ALL DONE $(date)"
